# Detox
Accurate determination of node and arc multiplicities in a de Bruijn graph using Conditional Random Fields.

## Prerequisites
In order to build Detox you need to the following software or library packages:

  * CMake version 2.6.3 or higher
  * GCC version 4.7 or higher
  * Google's [sparsehash](https://github.com/sparsehash/sparsehash)
  * ZLIB (optional)

In order to run Detox you also need to install and compile [BCALM 2](https://github.com/GATB/bcalm).

## Compilation

Clone the Detox Github repository:

```bash
git clone https://github.com/biointec/detox
```

Next, compile the C++ code as follows:

```bash
cd detox
mkdir build
cd build
cmake ..
make
```

After compilation, the `detox` binary can be found in:

```bash
detox/build/src/detox
```

## Usage

We use BCALM 2 to build a compacted de Bruijn graph from sequencing data. Next, Detox is used to infer its the node and arc multiplicities.
### 1. Generating the de Bruijn graph using BCALM 2
The first step is to prepare a manifest file, e.g. named `reads.mf`, in which the input FASTQ files are listed. For example, suppose the reads are stored in two FASTQ files named `reads_1.fastq` and `reads_2.fastq`, then the `reads.mf` file would simply list these input files, one file per line:

```
reads_1.fastq
reads_2.fastq
```
The set of reads can be in FASTA (`.fa`/`.fasta`) or FASTQ (`.fq`/`.fastq`) format or a gzipped version (`.gz`) of these formats.

Using BCALM 2, one can then create a de Bruijn graph as follows:

```bash
bcalm -in reads.mf -kmer-size 21 -abundance-min 2
```

Depending on the size of the genome, and the number of CPU cores of your machine, thus process may take minutes (for bacterial genomes) or many hours (for human genome scale genomes). The command line option `-abundance-min 2` instructs BCALM 2 to only retain k-mers that occur at least twice in the input reads, thus removing a large fraction of the sequencing errors. Remark that BCALM 2's default value for `-abundance-min` is 2. If you want a de Bruijn graph containing all k-mers in your read set to run through Detox, you should set this value to 1.

When BCALM 2 is finished, it produces the file `reads.unitigs.fa`, in which the de Bruijn graph's unitig node sequences as well as their connectivity are encoded in a FASTA-format. This file serves as input for Detox.

### 2. Inferring the node/arc multiplicities using Detox

The pipeline takes as input a de Bruijn graph in the FASTA-format of BCALM 2 (`reads.unitigs.fa`) as well as the original read set (`reads.mf`). Detox can then be run as follows:

```bash
detox -use-qual -abundance-min 2 -crf-nb-size 3 reads.unitigs.fa reads.mf
```

The two mandatory input arguments are `reads.unitigs.fa` and `reads.mf`, in that order. All other flags or command line options should be listed prior to these two input arguments.

The flag `-use-qual` instructs Detox to weigh the k-mer counts using the Phred quality scores. When BCALM 2 was run using the `-abundance-min`  option, it is best to provide this information also to Detox. When Detox is fitting the error model to the data, it then knows that all k-mers that occur fewer times than the value provided by `abundance-min` are missing. It takes this into account when fitting the model to the data. If you use Detox with unweighted k-mer counts, the value of `-abundance-min` can be automatically inferred.

Finally, the `-crf-nb-size` option specifies the size of the neighborhood subgraph around the node/arc whose multiplicity is inferred. For smaller, bacterial genomes, a value of 5 can be chosen in order to get the highest possible accuracy. For human genome-scale genomes, a value of 1 is recommended to prevent very long runtimes. In most cases, using a small neighborhood size is enough, however, if regions of extremely high or low coverage span large regions of the de Bruijn graph, it might take a large neighborhood size to obtain a CRF that can infer the correct multiplicities.

Upon completion, the estimated node and arc multiplicities are stored in files `estmult.node` and `estmult.arc`, respectively. The file `estmult.node` contains the following fields (tab-separated):

* The node identifier (integer value from 1 to #nodes)
* The estimated multiplicity m (non-negative integer value)
* The log-probability _log[ P(mult. is m) ]_
* The average node coverage (weighed by the Phred quality use if the `-use-qual` flag was used)
* The length of the node expressed as the number of k-mers.

The file `estmult.arc` contains the following fields (tab-separated):

* The node identifier of the node from which the arc is originating
* The node identifier of the destination node
* The estimated multiplicity m (non-negative integer value)
* The log-probability _log[ P(mult. is m) ]_

Note that the node identifiers in the `estmult.*` files are equal to (the BCALM node identifier + 1). In other words, we number nodes 1, 2, 3, ... whereas BCALM start numbering from 0. A negative node identifier refers to the reverse-complement node.

### 3. Advanced user information 
The Detox pipeline consists of 3 stages. Each stage outputs a number of intermediate files. When re-running Detox (e.g. using different parameter settings), the existence of these intermediate files is checked. If they exist, stages 1 or 2 may be skipped. In order to force re-running stages 1 or 2, simply remove the appropriate files: `rm *.st1` for stage 1 or `rm *.st2` for stage 2. We will now explain the stages in more detail below.

#### Stage 1
The input file `reads.unitigs.fa`, produced by BCALM 2 is read. Even though BCALM 2 provides average k-mer counts for each node (unitig), it does not provide (k+1)-mer counts for the arcs. Also, when using the `-use-qual` flag, Detox needs to weigh the k-mer counts using the Phred quality scores. Note that the Phred score ASCII base value can be set using `-phred-base <value>` (default=33).

So in either case, after reading the de Bruijn graph into memory, Detox needs to stream through all reads in order to establish both the appropriate node and arc counts. The results are written to `reads.unitigs.fa.st1`, a binary file containing the de Bruijn graph annotated with q- or k-mer coverage. This file is a prerequisite for stage 2.

If the file `genome.fasta` is present in current working directory, Detox will assume it contains the reference genome and it will also compute the true node and edge multiplicities. These will be stored in files `truemult.node` and `truemult.edge`, respectively.

In principle, one should only re-run stage 1 when switching between k-mer or q-mer counts.

#### Stage 2
In stage 2, a k-mer and (k+1)-mer mixture model is fitted to the data. By default, Detox selects a number of random nodes and random arcs from the de Bruijn graph on which the model is then trained. The number of nodes and arcs can be specified using the `-em-train-size` (default = 10.000) option. Larger values lead to larger runtimes, but model parameters will be estimated more accurately.

Both models are fitted using the expectation-maximization (EM) algorithm. During the E-step, the multiplicity is inferred for selected nodes and arcs based on the current model estimation. During the M-step, the model parameters are updated based on the estimated multiplicities. This process continues until either the model has converged (the relative convergence tolerance can be set using `-em-conv-eps` (default = 0.001) or until the maximum number of EM iterations is reached (this value can be set using the `-em-max-iter` (default = 25) option.

Stage 2 generates following files:

  * `nodemodel.st2` and `edgemodel.st2` that contain the parameters of the fitted models to the k-mer or q-mer data for the nodes and arcs respectively. These files are prerequisites for stage 3.
  * `nodes.gnuplot`, `nodes.dat`, `edges.gnuplot` and `edges.dat`. These are auxiliary files that can be used to plot the model fit to the histogram data. It can be used to manually check if the fit is good.  You will need the `gnuplot` for that. Simply run `gnuplot nodes.gnuplot` and `gnuplot edges.gnuplot` to generate the files `nodesspectrum.pdf` and `edgesspectrum.pdf`, respectively.

The following command line options are relevant to set the **initial** parameter estimates to the q-mer or k-mer models. It is generally not necessary to provide them. Only when the EM algorithm does not converge to the correct solution, manually setting these these initial estimates might help.

* `-mm-coverage` (default = auto) to provide the average coverage of the non-repeated nodes/arcs of the dataset.
* `-mm-err-cov` (default = 1) to provide the average coverage of erroneous nodes or arcs. 
* `-mm-odf` (default = 1.5) to provide the overdispersion factor of the Negative Binomial distribution. This is the ratio of the variance and the mean.

The number of mixture model components is controlled using the `-mm-components` option. By default, we use 6 components: one to model the erroneous nodes (resp. arcs) and 5 components to model multiplicities 1 to 5.

It is important to note that during the E-step, the multiplicities of nodes and arcs are inferred using the CRF methodology. Therefore, the `-crf-nb-size` option to specify the neighborhood size is relevant. In general, a value larger than 0 (i.e., using the CRF methodology) reduces the required number of EM iterations and proves more robust against poor choices of initial model parameter values.

One should only re-run stage 2 to re-train the mixture models. It suffices to remove model files: `rm *.st2`.

#### Stage 3
By default, in stage 3, the node and arc multiplicities are computed for **all** nodes and arcs in de Bruijn graph. One can list as subset of nodes and edges in files `nodes.list` and `edges.list`, respectively. If those files exist, the node and edge multiplicities are only inferred for the entries listed in those files.

The most important command line option is the `-crf-nb-size` to specify the neighborhood size. Note that is possible use a different neighborhood size in stage 2 (to train to the models) and stage 3 (to infer the multiplicities).

#### Other parameters that are relevant for CRF inference:

  * `-crf-margin` (default = 2) specifies the number of alternative multiplicities (one-sided) used within the probabilistic factors. For example, for the default value of 2, if the estimated multiplicity of a node is 3 based on its coverage, then the probabilities will be computed for multiplicities {1,2,3,4,5}. If the estimated multiplicity of a node is 0, then probabilities will be computed for multiplicities {0,1,2}.
  * `-crf-flow` (default 1e7) defines the probability ratio for the scenario in which the conservation of flow holds and the scenario in which it does not hold. The higher this value, the stronger the CRF model will impose the conservation of flow rule. If you have a dataset with a very high coverage (much larger than _50x_) you might set this parameter to a higher value, as multiplicity distributions that are far apart require a stronger conservation of flow rule.
  * `-crf-max-fact` (default 1e6) defines the maximum factor size that is allowed during the variable elimination algorithm. If this threshold is exceeded, the algorithm with reduce the `-crf-nb-size` (only for the node/arc for which this threshold was exceeded).

 ### 4. Other options

 * `-num-threads` to specify the number of threads to use. By default, the number of threads equals the number of CPU cores or twice that number for CPUs that support hyperthreading.
 * `-cyt-graph <node ID>`  to produce files that can be imported in Cytoscape to visualize part of the de Bruijn graph (i.e. a neighborhood centered around node specified by `<node ID>`. The size of the neighborhood is specified by the `-crf-nb-size` option.
 * `-help` to display the help page.
 * `-version` to display the version number.

 ### 5. Checking the fit of the model produced in stage 2

If you want to check the fit to the histogram of the model that Detox fits in stage 2, this can be done visually using [Gnuplot](http://www.gnuplot.info/).
Intermediate files `nodes.gnuplot`, `nodes.dat`, `edges.gnuplot` and `edges.dat` can be used to plot the histogram overlaid with the estimated model. Just run
```
gnuplot nodes.gnuplot
```
This will produce a pdf output of the plot.
Detox auto-determines a good scope for the plot, but you can manually change this in the `.gnuplot` files if need be.

#### What to do when the fit seems bad

We have determined what we think are the best default settings for the parameters Detox uses, based on tests on several datasets. However, if your visual inspection of the model fit seems bad, you can try to manually alter some of these parameter values:

##### Check if convergence was obtained in stage 2

By default, the maximum number of EM-iterations is 25. You will be notified by Detox when stage 2 finished with this many iterations. If you see that there was still a relatively large change in parameter values between iteration 24 and 25, it can be beneficial to rerun stage 2 with a higher number of EM-iterations. Option `-em-max-iter [iters]` sets the maximum number of EM-iterations.

##### Adjusting the initial estimates for the parameters of the error distribution and of the multiplicity 1 distribution

If you see that peaks of the estimated distributions do not correspond with the correct peaks in the histogram you can initialize the means of the error- and the multiplicity 1 distribution yourself (`-mm-coverage`, resp. `-mm-err-cov`). Another option is to initialize the multiplicity 1 mean to the average read coverage of your dataset (if this is known to you). If the fit of the distributions seems too wide or too tight, try a different initialization of the overdispersion factor (`-mm-odf`).

 ### 6. Visualizing parts of the de Bruijn graph

You can visualize a de Bruijn graph neighborhood around a central node, as is used to construct the CRF that determines the multiplicity of that node, using [Cytoscape](https://cytoscape.org/). Just run:
```
/Path/To/detox/src/detox -crf-nb-size [nb-size] -cyt-graph [nodeID] reads.unitigs.fa reads.mf
```
This will run an alternative **stage 3** that extracts a subgraph (of size `nb-size`) of the de Bruijn graph that is used in the CRF model that determines the multiplicity of the central node (with identifier `nodeID`). A CRF for a neighborhood of size `nb-size` is then used to determine multiplicities for all nodes and arcs in the subgraph.
If there are no `.stX` files present in the working directory, stage X will be rerun as well.

The output is provided as two files `cytgraph.<nodeID>.nodes` and `cytgraph.<nodeID>.arcs`. Within Cytoscape you first import the `.arcs` file with `File > Import > Network from file...` and then you import the `.nodes` file with `File > Import > Table from file...`
