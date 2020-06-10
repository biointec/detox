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

#include "rawfile.h"

using namespace std;

// ============================================================================
// FASTA FILE
// ============================================================================

bool RawFile::getNextRead(std::string &read)
{
        // empty output strings
        read.clear();

        // read until something non-empty is encountered
        if (!rfHandler->good())
                return false;

        read = rfHandler->getLine();
        if (!read.empty() && read.back() == '\n')
                read.pop_back();

        return !read.empty();
}

bool RawFile::getNextRecord(ReadRecord& record)
{
        // empty output strings
        record.clear();

        // read until something non-empty is encountered
        if (!rfHandler->good())
                return false;

        record.read = rfHandler->getLine();
        if (!record.read.empty() && record.read.back() == '\n') {
                record.postRead = '\n';
                record.read.pop_back();
        }

        return !record.read.empty();
}
