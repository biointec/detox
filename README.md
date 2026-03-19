# Detox
An experimental, but functional de novo genome assembler for Illumina data. Detox uses Conditional Random Fields (CRFs) to estimate node and edge multiplicities in the de Bruijn graph. On some datasets, it is able to outperform the state-of-the-art SPAdes assembler in terms of NGA50 (higher contig length) and accuracy (fewer misassemblies, higher per-base accuracy).

For the versions of Detox related to our earlier papers that focussed on multiplicity estimation alone, please refer to the `approx_inference` branch and the releases.

## Citing
If you make use of Detox, please cite us appropriatly.
A publication describing the assembly pipeline and the use of MPNNs is currently under review. 

> A. Steyaert, S. Heemeryck, and J. Fostier, "Deep learning based de Bruijn graph multiplicity estimation improves short read genome assembly", UNDER REVIEW

Until this is published, please refer to: 

> Steyaert, A., Audenaert, P. & Fostier, J. Accurate determination of node and arc multiplicities in de bruijn graphs using conditional random fields. BMC Bioinformatics 21, 402 (2020). https://doi.org/10.1186/s12859-020-03740-x

> A. Steyaert, P. Audenaert and J. Fostier, "Improved Node and Arc Multiplicity Estimation in De Bruijn Graphs Using Approximate Inference in Conditional Random Fields," in IEEE/ACM Transactions on Computational Biology and Bioinformatics, vol. 20, no. 3, pp. 1995-2006, 1 May-June 2023, doi: 10.1109/TCBB.2022.3229085.

## Prerequisites
In order to build Detox you need to the following software or library packages:

  * CMake version 2.8.12 or higher
  * GCC version 9 or higher (or any other compiler with C++17 support)
  * Google's [sparsehash](https://github.com/sparsehash/sparsehash)
  * ZLIB (optional)

Note that Google's sparsehash can typically be installed used the Linux package manager (`apt`, `aptitude` or `dnf`). For example, in Ubuntu you can run `sudo apt install libsparsehash-dev`.

To build the source code:

 1. `git clone https://github.com/biointec/detox detox`
 2. `cd detox`
 3. `mkdir build`
 4. `cd build`
 5. `cmake ..`
 6. `make`

The binary to run detox is in `build/src/detox`.

## Basic usage

### Input
The pipeline takes as input a set of reads. The set of reads can be in fasta, fastq, sam or raw format or a gzipped version of any of these formats. These files are passed via a manifest file containing a list of multiple reads-files (e.g. `reads.mf`) (one filename/line). 

To make use of paired-end reads, a manifest file containing the names of both files on the same line (`/path/to/reads.1.fastq    /path/to/reads.2.fastq`), separated by a tab should be passed. If both filenames are not on the same line, both read sets will be used, but will not be interpreted as paired-end reads.

### Running Detox

You can simply run detox without passing any options values. Only pass the reads file(s) (`reads.mf`):
```
/Path/To/detox/src/detox reads.mf
```
All options are described below.

### Output 

Detox produces some intermediate files (see further for what their purpose is). The final contigs are written in fasta format to a file `detox.fasta`, while the final graph - including multiplicity estimates for its nodes and arcs is written in cystoscape format to `graph.arcs` and `graph.nodes`.
      
## Using quality scores

By default, Detox estimated a coverage spectrum based on the q-mer coverage. Contributions of each k-mer are  weighed by the quality scores of the nucleotides. This behaviour can be switched off with the `-no-qual` option:

**IMPORTANT**: Detox saves only either k-mer coverage or q-mer coverage in stage 1. It is thus necessary to rerun Detox from stage 1. To ensure this, please remove the intermediate files containing `.stX` extension (`X` being either 1-4) or move them to another folder.

## Running Stage 4 with MPNN-based multiplicity inference

As an alternative to CRF-based multiplicity computation in stage 4, it is possible to pass a pre-trained ANN that predicts multiplicities based on coverage (and de Bruijn graph topology if using a GNN).
To this end we make use of the [libtorch library](https://docs.pytorch.org/cppdocs/installing.html)
Make sure this library is located in your path, or alternatively pass the appropriate location to cmake when compiling detox:
```
cmake -DCMAKE_PREFIX_PATH=/path/to/libtorch ..
```

Detox can then be run with the option `-nn-model path/to/model_GNN.pt` .  
Python code with the MPNN model, example training flow and command to save the model to the correct format can be found in the `gnn_training` folder of this repository.
 
## Intermediate files
 
  * Stage 1 outputs a set of files with suffix `*.st1`. These are binary files containing the constructed de Bruijn graph annotated with q- or k-mer coverage. This file is a prerequisite for stage 2.
  * Stage 2 outputs the initial mixture model representing the relationship between coverage and multiplicity. The model parameters are written to  `model.node.st2` and `model.edge.st2`. These files are prerequisites for Stage 3. Files to visualise these models and their fit to the histogram, with gnuplot are output as well (`histogram.*`)
  * Stage 3 outputs a binairy file with the cleaned de Bruijn graph (`dbgraph.bin.st3`), as well as the re-estimated mixture model. These files are prerequisites for Stage 4.
  * Stage 4 outputs, besides the output files described above, the (paired-end) read alignments used during contig extraction, as well as a binary version of the final de Bruijn graph. If you want to re-run stage 4 (e.g., to compare CRF-based with GNN-based output) you should remove this file (`rm *.st4`). 
  
 
 ## Options
 
 ### Important options/parameters and their influence
 
 * do not use quality scores `-no-qual`   
    Use quality score weighted contributions to the coverage of a node/arc. 
    The phred score ASCII base can be set with `-phred-base <value>` (default=33).
 * neighourhood size `-crf-nb-size <size>` (default = 5)   
    This is the size of the neighbourhood used to create CRF for the determination of the multiplicity of a node or arc.
    A larger value of this parameter will result in a slight increase in accuracy, however it also strongly influences the runtime.
 * flow conservation strength `-crf-flow <strength>` (default = 1.0e7)   
    This parameter determines the difference in value that the flow-conservation factors assign to multiplicity combinations for which conservation of flow holds, vs. combinations for which it is violated.
    A large value often results in a higher accuracy, but the difference in accuracy became smaller for values of 1.0e5 and higher in our experiments.
 * EM training set size `-em-train-size <value>` (default = 1e4)   
    Higher values will result in less variability in parameter estimates and accuracy between runs, but it will also increase runtimes. Experiments showed that a value of 1e4 resulted in small enough variability and low runtimes such that the runtime is dominated by Stage 3.
 * Parameters relating to path resolution during assembly, both parameters will influence how much (consistent) support is required to keep paths supported by paired-end reads: 
      * `max-confl-degree` (default=0.25): maximum conflict degree (between paths)
      * `min-count` (default=2): minimum support from paired-end reads required
    
### Other options
Note that the usage info might give even more options, however, if they are not listed here, they are not yet fully supported!

  * margin around prior multiplicity estimate `-crf-margin <value>` (default = 2)   
    For value `x` and prior multiplicity estimate `m` the model will consider all values in `[m-x,m+x]`. The size of the interval will influence exact inference runtime, as this increases/decreases the number of computations needed.
  * maximum intermediate factor size `-crf-max-fact <value>` (default = 1.0e6)   
    If the number of computations during variable elimination would exceed this value, computations are restarted with a lower neighbourhood size.
  * Set initial coverage estimate `-mm-coverage <value>` (default = auto-infer)
  * Set initial error coverage estimate `-mm-err-cov <value>` (default = 1.0)
  * Set initial overdispersion factor estimate `-mm-odf <value>` (default = 1.5)
  * Set number of components in the mixture model `-mm-components <value>` (default = 6)
  * Set maximum number of EM-iterations `-em-max-iter <value>` (default = 25)
  * relative EM convergence criterion `-emconv-eps <value>` (default = 1.0e-3)   
    This is the minimum fraction of multiplicity estimates that should change to determine convergence has not yet been reached.


## [ADVANCED USER INFO] Checking the fit of the model produced in stage 2

If you want to check the fit to the histogram of the model that Detox fits in stage 2, this can be done visually using [Gnuplot](http://www.gnuplot.info/).
Intermediate files `histogram.node.gnuplot`, `histogram.node.dat`, `histogram.edge.gnuplot` and `histogram.edge.dat` can be used to plot the histogram overlaid with the estimated model. Just run
```
gnuplot histogram.node.gnuplot
```
This will produce a pdf output of the plot.
Detox auto-determines a good scope for the plot, but you can manually change this in the `.gnuplot` files if need be.

### What to do when the fit seems bad

We have determined what we think are the best default settings for the parameters Detox uses, based on tests on several datasets. However, if your visual inspection of the model fit seems bad, you can try to manually alter some of these parameter values:

#### Check if convergence was obtained in stage 2

By default, the maximum number of EM-iterations is 25. You will be notified by Detox when stage 2 finished with this many iterations. If you see that there was still a relatively large change in parameter values between iteration 24 and 25, it can be beneficial to rerun stage 2 with a higher number of EM-iterations. Option `-em-max-iter [iters]` sets the maximum number of EM-iterations.

#### Adjusting the initial estimates for the mean of the error distribution and of the multiplicity 1 distribution

By default the error distribution mean is initialised to 1.0, while the multiplicity 1 mean is initialised to an automatically inferred estimate of the average k-mer coverage. If you see that peaks of the estimated distributions do not correspond with the correct peaks in the histogram you can initialise these two distribution means yourself to where you see the center of the histogram peaks are. Another option is to initialise the multiplicity 1 mean to the average read coverage of your dataset (if this is known to you).

The intial estimate for the mean of the error distribution can be set with the option `-mm-err-cov [value]`. The initial estimate for the mean of the multiplicity 1 distribution can be set with the option `-mm-coverage [value].`

### What to try when your fit seems OK, but you think your downstream results could be better

If your read set is large, the subset size used for EM-training in stage 2 might be too small. In this case you Detox probably choose a subset that did not represent the full dataset well. If you re-run stage 2 and you see a big change in model estimates, this was probably the case. You can change the subset size with option `-em-train-size <size>` to a larger value (default value is `1e4`).

If you think there is little to no difference in results when using Detox as compared to just using a cut-off in the coverage-spectrum, There are two things your can try. 

First, if you are okay with a larger runtime, try using a larger neighbourhood size (option `-crf-nb-size [nb-size]`). If regions of extremely high or low coverage span large regions of the de Bruijn graph, it might take a large neighbourhood size to obtain a CRF that can detect the correct multiplicities.

Secondly, you can alter the strength of the flow conservation factors. This could  be beneficial for example when you have a high coverage dataset and the strength of beliefs based on coverage alone is high. Flow conservation strength can be set with `-crf-flow <strength>`, where strength is a positive number (default value is `1e7`).