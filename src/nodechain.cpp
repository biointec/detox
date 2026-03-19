/******************************************************************************
 *   Copyright (C) 2014 - 2023 Jan Fostier (jan.fostier@ugent.be)             *
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

#include <cassert>
#include <climits>

#include "nodechain.h"

using namespace std;

// ============================================================================
// NODE CHAIN CLASS
// ============================================================================

ostream &operator<<(ostream &out, const NodeChain& nc)
{
        if (nc.empty())         // nothing to do...
                return out;

        for (size_t i = 0; i < nc.size(); i++)
                out << nc[i] << "\t";
        out << "(" << nc.getCount() << ")";
        return out;
}

// ============================================================================
// CONSENSUS NODE CHAIN CLASS
// ============================================================================

ostream &operator<<(ostream &out, const ConsNodeChain& nc)
{
        if (nc.empty())         // nothing to do...
                return out;

        out << nc[0].first << " (" << nc[0].second << ")";
        for (size_t i = 1; i < nc.size(); i++)
                out << "\t" << nc[i].first << " (" << nc[i].second << ")";
        return out;
}

ConsNodeChain ConsNodeChain::merge(const ConsNodeChain& X,
                                   const ConsNodeChain& Y,
                                   float MCD)
{
        const int MIN = std::numeric_limits<int>::min();
        int m = X.size(), n = Y.size();

        int bestOC, bestScore = 0;
        vector<int> bestScoreX, bestScoreY;
        for (int OC = -n + 1; OC < m; OC++)
        {
                int i = (OC < 0) ? 0 : OC;
                int j = (OC < 0) ? -OC : 0;

                vector<int> scoreX, scoreY;
                scoreX.reserve(min<int>(m - i, n - j) + 1);
                scoreY.reserve(min<int>(m - i, n - j) + 1);
                scoreX.push_back(OC < 0 ? MIN : 0);     // Y is leading
                scoreY.push_back(OC > 0 ? MIN : 0);     // X is leading

                for ( ; (i < m) && (j < n); i++, j++)
                {
                        const int& vX = X[i].second;
                        const int& vY = Y[j].second;

                        // in case of a match
                        if (X[i].first == Y[j].first) {
                                int prev = max(scoreX.back(), scoreY.back());
                                scoreX.push_back(prev + min(vX, vY));
                                scoreY.push_back(prev + min(vX, vY));
                                continue;
                        }

                        // in case of mismatch
                        float CDX = float(vY) / float(vX + vY);
                        scoreX.push_back((CDX > MCD) || (scoreX.back() == MIN) ?
                                MIN : scoreX.back() - vY);
                        float CDY = float(vX) / float(vX + vY);
                        scoreY.push_back((CDY > MCD) || (scoreY.back() == MIN) ?
                                MIN : scoreY.back() - vX);

                        // stop early
                        if ((scoreX.back() == MIN) && (scoreY.back() == MIN))
                                break;
                }

                int score = (i < m) ? scoreX.back() :   // X is trailing
                           ((j < n) ? scoreY.back() :   // Y is trailing
                            max(scoreX.back(), scoreY.back()));

                // keep track of the best overlap
                if (score > bestScore) {
                        bestScore = score;
                        bestOC = OC;
                        bestScoreX = scoreX;
                        bestScoreY = scoreY;
                }
        }

        if (bestScore == 0)     // no alignment with positive score...
                return ConsNodeChain();

        int i = min(m, n + bestOC);
        int j = min(n, m - bestOC);
        assert( (i == m) || (j == n) );

        // trailing part of the merged node chain
        ConsNodeChain res = (i < m) ?
                ConsNodeChain(X.begin() + i, X.end()) :
                ConsNodeChain(Y.begin() + j, Y.end());

        // overlapping part of the merged node chain
        bool atX = (i < m);     // trailing part of X?
        bool atY = (j < n);     // trailing part of Y?
        while ((i > 0) && (j > 0)) {
                i--; j--;
                if (!atX && !atY) {     // choose path based on score
                        atX = bestScoreX.back() >= bestScoreY.back();
                        atY = !atX;
                }
                bestScoreX.pop_back(); bestScoreY.pop_back();
                res.push_front(atX ? X[i] : Y[j]);

                // at a match, we can switch chains
                if (X[i].first == Y[j].first) {
                        atX = atY = false;
                        res.front().second = X[i].second + Y[j].second;
                }
        }

        // leading part of the merged node chain
        if (i > 0)
                res.insert(res.begin(), X.begin(), X.begin() + i);
        else
                res.insert(res.begin(), Y.begin(), Y.begin() + j);

        return res;
}
