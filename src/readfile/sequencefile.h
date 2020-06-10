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

#ifndef SEQUENCEFILE_H
#define SEQUENCEFILE_H

#include "readfile.h"

class SequenceFile : public ReadFile
{

public:
        /**
         * Default constructor
         * @param gzipped True if the file is gzipped
         */
        SequenceFile(bool gzipped) : ReadFile(gzipped) {};

        /**
         * Get the next read from a file
         * @param read String containing the read (output)
         * @param description String containing the read description (output)
         */
        bool getNextRead(std::string &read, std::string &description);

        /**
         * Write a read to file
         * @param read Read to write (input)
         * @param description Description of the read (input)
         */
        void writeRead(const std::string &read,
                       const std::string &description);
};

#endif
