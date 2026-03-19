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

#include "dsnode.h"

using namespace std;

bool DSNode::deleteLeftArc(NodeID targetID)
{
        for (int i = 0; i < arcInfo.p.numLeft; i++) {
                if (lArc[i].getNodeID() != targetID)
                        continue;

                // shift remaing arcs one position down
                for (i++; i < arcInfo.p.numLeft; i++)
                        lArc[i - 1] = lArc[i];

                arcInfo.p.numLeft--;
                lArc[arcInfo.p.numLeft].deleteArc();
                return true;
        }

        return false;   // arc is not found, get out
}

bool DSNode::deleteRightArc(NodeID targetID)
{
        for (int i = 0; i < arcInfo.p.numRight; i++) {
                if (rArc[i].getNodeID() != targetID)
                        continue;

                // shift remaing arcs one position down
                for (i++; i < arcInfo.p.numRight; i++)
                        rArc[i - 1] = rArc[i];

                arcInfo.p.numRight--;
                rArc[arcInfo.p.numRight].deleteArc();
                return true;
        }

        return false;   // arc is not found, get out
}

// ============================================================================
// HIGHER-ORDER NODE IDENTIFIER
// ============================================================================

ostream& operator<<(ostream& out, const HoNodeID &id)
{
        for (size_t i = 0; i < id.size(); i++)
                out << id[i] << " ";
        return out;
}

// ============================================================================
// HIGHER-ORDER EDGE IDENTIFIER
// ============================================================================

ostream& operator<<(ostream& out, const HoEdgeID &id)
{
        cout << id.first << " --> " << id.second;
        return out;
}
