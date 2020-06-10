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

#ifndef FASTAFILE_H
#define FASTAFILE_H

#include "readfile.h"

class FastAFile : public ReadFile
{

public:
        /**
         * Default constructor
         * @param gzipped True if the file is gzipped
         */
        FastAFile(bool gzipped) : ReadFile(gzipped) {};

        /**
         * Get the next read from a file
         * @param read String containing the read (output)
         */
        bool getNextRead(std::string &read);

        /**
         * Get the next record from a file
         * @param record Record to store output
         */
        bool getNextRecord(ReadRecord& record);

        /**
         * Write a record to file
         * @param record Record to write
         */
        void writeRecord(const ReadRecord& record);
};

#endif
