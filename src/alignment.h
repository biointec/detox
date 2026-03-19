/******************************************************************************
 *   Copyright (C) 2014 - 2020 Jan Fostier (jan.fostier@ugent.be)             *
 *   This file is part of Detox                                               *
 *                                                                            *
 *   This program is free software; you can redistribute it and/or modify     *
 *   it under the terms of the GNU Affero General Public License as published *
 *   by the Free Software Foundation; either version 3 of the License, or     *
 *   (at your option) any later version.                                      *
 *                                                                            *
 *   This program is distributed in the hope that it will be useful,          *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of           *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            *
 *   GNU General Public License for more details.                             *
 *                                                                            *
 *   You should have received a copy of the GNU Affero General Public License *
 *   along with this program; if not, see <https://www.gnu.org/licenses/>.    *
 ******************************************************************************/

#ifndef ALIGNMENT_H
#define ALIGNMENT_H

#include "string.h"
#include "global.h"

#include <iostream>
#include <vector>
#include <limits>

// ============================================================================
// CLASS PROTOTYPES
// ============================================================================

class NodeChain;

// ============================================================================
// CIGAR STATES
// ============================================================================

enum CIGARstate { M, I, D, NOTHING };

// ============================================================================
// ALIGNMENT RESULT
// ============================================================================

class AlnRes {

//private:
public:
        /**
         * s1len and s2len can be shorter than the corresponding sequence length
         * if the alignment contains trailing gaps in either of the sequences.
         * In that case s1len/s2len contain the sequence length that were
         * aligned to characters or internal gaps in the other sequence.
         */

        int s1len;      // number of characters in s1 that were actually aligned
        int s2len;      // number of characters in s2 that were actually aligned
        int ib, ie;     // part of s1 that is aligned, i.e., s1[ib, ie[
        int jb, je;     // part of s2 that is aligned, i.e., s2[jb, je[
        int score;      // alignment score
public:
        AlnRes() : s1len(0), s2len(0), score(0) {}

        AlnRes(int s1len, int s2len, int score) :
                s1len(s1len), s2len(s2len), score(score) {}

        AlnRes(int ib, int ie, int jb, int je, int score) :
                ib(ib), ie(ie), jb(jb), je(je), score(score) {}
};

// ============================================================================
// ALIGNMENT RESULT 2
// ============================================================================

class AlnRes2 {

public:
        int lenX;       // number of characters in X that were actually aligned
        int lenY;       // number of characters in Y that were actually aligned
        int score;      // alignment score
        int maxAtt;     // maximum attainable score

public:
        AlnRes2() : lenX(0), lenY(0), score(0), maxAtt(0) {}

        AlnRes2(int lenX, int lenY, int score, int maxAtt) :
                lenX(lenX), lenY(lenY), score(score), maxAtt(maxAtt) {}
};

// ============================================================================
// MATRIX CLASS
// ============================================================================

template<class T>
class Matrix {

private:
        size_t numRow;          // number of actively used rows
        size_t numCol;          // number of actively used columns
        std::vector<T> data;    // actual data

public:
        /**
         * Default constructor (initializes to 100x100 matrix)
         */
        Matrix() : numRow(0), numCol(0) {}

        /**
         * Default constructor
         * @param numRow Number of rows
         * @param numCol Number of columns
         */
        Matrix(size_t numRow, size_t numCol) : numRow(numRow), numCol(numCol) {
                data.resize(numRow * numCol);
        }

        /**
         * Resize the current matrix, allocate memory only if necessary
         * @param numRow Number of rows
         * @param numCol Number of columns
         */
        void resize(size_t numRow, size_t numCol) {
                this->numRow = numRow;
                this->numCol = numCol;
                if (data.size() < (numRow * numCol))
                        data.resize(numRow * numCol);
        }

        /**
         * Set all elements of a matrix to a target value
         * @param target Target value
         */
        void fill(const T& target) {
                for (size_t i = 0; i < numRow * numCol; i++)
                        data[i] = target;
        }

        /**
         * Operator () overloading to access element (i, j) of the matrix
         * @param i Row index
         * @param j Column index
         * @return Copy of the element at position (i, j)
         */
        T operator() (int i, int j) const {
                return data[i * numCol + j];    // row-major
        }

        /**
         * Operator () overloading to access element (i, j) of the matrix
         * @param i Row index
         * @param j Column index
         * @return Reference to element at position (i, j)
         */
        T& operator() (int i, int j) {
                return data[i * numCol + j];
        }

        /**
         * Compute the determinant of the matrix
         * @return The determinant
         */
        T determinant() {
                Matrix<T>& A = *this;

                // compute the LU decomposition
                for (int i = 0; i < numRow-1; i++) {
                        for (int j = i+1; j < numRow; j++) {
                                double m = A(j,i) / A(i,i);
                                for (int k = i; k < numRow; k++)
                                        A(j,k) = A(j,k) - m * A(i,k);
                        }
                }

                T retVal = 1;
                for (int i = 0; i < numRow; i++)
                        retVal *= A(i,i);

                return retVal;
        }
};

// ============================================================================
// ALIGNMENT CLASS
// ============================================================================

class NWAligner {

private:
        size_t currSize;// size of matrix (as a linear array)
        int numRow;

        int maxIndel;   // maximum number of indels
        int M;          // match score
        int I;          // mismatch penalty
        int G;          // gap score
        int *matrix;    // alignment matrix

        Matrix<int> D;  // dense score matrix

        /**
         * Allocate memory to align sequences with lengths l1 and l2
         * @param l1 sequence length 1
         * @param l2 sequence length 2
         */
        void reserveBanded(size_t l1, size_t l2);

        /**
         * Return the alignment score when aligning two characters
         * @param a character a in ACTG+N alphabet, N = wildcard
         * @param b character b in ACTG+N alphabet, N = wildcard
         * @return alignment score (diagonal)
         */
        int S(char a, char b) const;

public:
        /**
         * Default constructor
         * @param maxIndel Maximum number of insertions or deletions
         * @param M Match score
         * @param I Mismatch penalty
         * @param G Gap score
         */
        NWAligner(int maxIndel = 3, int M = 1, int I = -1, int G = -3);

        /**
         * Allocate memory to align a pattern up to length m
         * @param m Pattern length
         * @param startScore Start score
         */
        void reserveBanded(size_t m, int startScore);

        /**
         * Destructor
         */
        ~NWAligner() {
                delete [] matrix;
        }

        /**
         * Perform the alignment between two sequences
         * @param s1 First string to align
         * @param s2 Second string to align
         * @return The alignment score (higher is better)
         */
        AlnRes align(const std::string &s1, const std::string &s2);

        /**
         * Get the match score (positive number)
         * @return The match score
         */
        int getMatchScore() const {
                return M;
        }

        /**
         * Get the gap score (negative number)
         * @return The gap score
         */
        int getGapScore() const {
                return G;
        }

        /**
         * Print matrix to stdout
         */
        void printAlignment(const AlnRes& alnRes,
                            const std::string &s1,
                            const std::string &s2) const;

        /**
         * Remove trailing gaps from the alignment
         * @param alnRes Alignment result input
         * @return Alignment result with gaps removed (if any)
         */
        AlnRes trimTrailingGaps(const AlnRes& alnRes) const;

        /**
         * Print matrix to stdout
         */
        void printMatrix(const AlnRes& alnRes) const;

        int operator() (int i, int j) const {
                int k = std::max(i, j);
                int l = maxIndel - i + j;
                return matrix[k * (2 * maxIndel + 1) + l];
        }

        int& operator() (int i, int j) {
                int k = std::max(i, j);
                int l = maxIndel - i + j;
                return matrix[k * (2 * maxIndel + 1) + l];
        }

        /**
         * Perform the banded alignment between two sequences
         * @param s1 First string
         * @param s2 Second string
         * @return The alignment score (higher is better)
         */
        int alignBanded(const std::string &s1, const std::string &s2);

        /**
         * Continue a banded alignment between two sequences
         * @param X Sequence X (horizontal)
         * @param Y Sequence Y (vertical)
         * @param offsetY Number of characters were already matched to X
         * in previous calls, alignment will be appended
         * @return The alignment result
         */
        AlnRes2 alignBandedContd(const std::string& X, const std::string& Y,
                                 int offsetY);

        /**
         * Get the maximum score / score position at a specific row
         * @param row Desired row index
         * @param sizeX Size of sequence X (horizontal)
         * @return pair <maxScore, column index>
         */
        std::pair<int, int> getMaxScoreAtRow(int row, int sizeX)
        {
                NWAligner& F = *this;   // shorthand notation
                std::pair<int, int> res(std::numeric_limits<int>::min(), 0);

                const int jBegin = std::max<int>(0, row - maxIndel);
                const int jMid   = std::min<int>(sizeX, row);
                const int jLast  = std::min<int>(sizeX, row + maxIndel);
                // in case of ex-aequo, report position closest to the middle
                for (int j = jBegin; j <= jMid; j++)
                        if (res.first <= F(row, j))
                                res = std::make_pair(F(row, j), j);
                for (int j = jMid + 1; j <= jLast; j++)
                        if (res.first < F(row, j))
                                res = std::make_pair(F(row, j), j);
                return res;
        }

        AlnRes2 alignBandedContdFullY(const std::string& X, const std::string& Y,
                                      int offsetY);

        /**
         * Global alignment but don't penalize trailing gaps in either s1 OR s2
         * @param s1 First string
         * @param s2 Second string
         * @return The alignment result
         */
        AlnRes alnGlobFreeEndGap(const std::string &s1, const std::string &s2);

        /**
         * Global alignment but don't penalize trailing gaps in s1
         * @param s1 First string
         * @param s2 Second string
         * @return The alignment result
         */
        AlnRes alnGlobFreeEndGapS2(const std::string &s1,
                                   const std::string &s2);

        /**
         * Get the maximal attainable score
         * @param l length of the sequence
         * @return maximal attainable score
         */
        int getMaxScore(size_t l) const {
                return M * l;
        }

        /**
         * Print matrix to stdout
         */
        void printAlignmentBanded(const std::string &s1,
                                  const std::string &s2) const;

        /**
         * Get the number of matches, substitutions and indels in the alignment
         * @param s1 First string
         * @param s2 Second string
         * @param nMatch Number of matches (output)
         * @param nSubst Number of substitutions (output)
         * @param nIndel Number of indels (output)
         */
        void getAlnStats(const std::string& s1, const std::string& s2,
                         size_t& nMatch, size_t& nSubst, size_t& nIndel) const;

        /**
         * Overlap alignment
         * @param s1 First string
         * @param s2 Second string
         * @return The alignment result
         */
        AlnRes overlapAln(const std::string &s1, const std::string &s2);

        /**
         * Overlap alignment between two node chains (score computation)
         * @param s1 Leftmost node chain
         * @param s2 Rightmost node chain
         * @return The alignment result (only ie, je and score are computed)
         */
        AlnRes overlapAln(const NodeChain& s1, const NodeChain& s2);

        /**
         * Overlap alignment between two node chains (backtracking)
         * @param s1 Leftmost node chain
         * @param s2 Rightmost node chain
         * @param AlnRes Alignment result (input/output: ib and jb are computed)
         * @param CIGAR CIGAR string (output)

         */
        void overlapBT(const NodeChain& s1, const NodeChain& s2, AlnRes& alnRes,
                       std::vector<std::pair<char, size_t>>& CIGAR) const;
};

#endif
