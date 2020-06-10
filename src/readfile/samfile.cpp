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

#include "samfile.h"

using namespace std;

// ============================================================================
// SAM FILE
// ============================================================================

bool SamFile::getNextRead(string &read)
{
        // empty output strings
        read.clear();

        // read past header files
        while (rfHandler->good() && rfHandler->peekCharacter() == '@')
                rfHandler->getLine();

        // read sequence identifier
        const char *result = rfHandler->getLine();
        istringstream oss(result);

        string qname, flag, rname, position, dummy;
        oss >> dummy >> flag >> rname >> position >> dummy >> dummy >>
                        dummy >> dummy >> dummy >> read;

        if (read.empty())       // end of file might be reached
                return false;

        return true;
}

bool SamFile::getNextRecord(ReadRecord& record)
{
        // empty output strings
        record.clear();

        // read past header files
        while (rfHandler->good() && rfHandler->peekCharacter() == '@')
                record.preRead.append(rfHandler->getLine());

        // read sequence identifier
        string result = rfHandler->getLine();

        if (result.empty())
                return false;

        istringstream oss(result);

        string qname, flag, rname, position, dummy;
        oss >> dummy >> flag >> rname >> position >> dummy >> dummy >>
                        dummy >> dummy >> dummy;

        record.preRead.append(result.substr(0, oss.tellg()));
        record.preRead.push_back('\t');
        oss >> record.read;
        record.postRead.append(result.substr(oss.tellg()));

        return !record.read.empty();
}
