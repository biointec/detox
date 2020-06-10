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

#ifndef READFILE_H
#define READFILE_H

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdio>
#include <cstdlib>
#include <vector>

#ifdef HAVE_ZLIB
        #include "zlib.h"
#endif

// ============================================================================
// ENUM TYPES
// ============================================================================

typedef enum { READ, WRITE } ReadFileMode;

// ============================================================================
// READ RECORD
// ============================================================================

class ReadRecord {

public:
        /**
         * Default constructor
         */
        ReadRecord() {}

        void clear() {
                preRead.clear();
                read.clear();
                qual.clear();
                postRead.clear();
        }

        std::string getRead() const {
                return read;
        }

        std::string& getRead() {
                return read;
        }

        std::vector<int> getNodeChain() const {
                return nodeChain;
        }

        std::vector<int>& getNodeChain() {
                return nodeChain;
        }

        size_t getReadLength() const {
                return read.length();
        }

        std::string getQStr() const {
                return qual;
        }

        std::string& getQStr() {
                return qual;
        }

//private:
        std::string preRead;    // everything in the record that precedes the read
        std::string read;       // read itself
        std::string qual;       // quality string
        std::string postRead;   // everything in the record that procedes the read

        std::vector<int> nodeChain; // node chain to which the read aligns
};

// ============================================================================
// READFILE HANDLER
// ============================================================================

class ReadFileHandler {

protected:
        ReadFileMode mode;

public:
        /**
         * Virtual destructor for your convenience
         */
        virtual ~ReadFileHandler() {};

        /**
         * Open a file with a given filename
         * @param filename File to open
         * @param mode Read (default) or write
         * @return True upon succes, false otherwise
         */
        virtual bool open(const std::string& filename,
                          ReadFileMode mode = READ) = 0;

        /**
         * Check if the input is still valid
         * @return True of false
         */
        virtual bool good() = 0;

        /**
         * Read a line from file
         * @return Pointer to an internal buffer containing the line
         */
        virtual const char *getLine() = 0;

        /**
         * Write a line to file
         * @param line String to write
         */
        virtual void writeLine(const std::string &line) = 0;

        /**
         * Write a character to file
         * @param c Character to write
         */
        virtual void writeChar(char c) = 0;

        /**
         * Peek at the next character in the stream
         * @return The character
         */
        virtual char peekCharacter() = 0;

        /**
         * Get one character from the stream
         * @return The character
         */
        virtual char getCharacter() = 0;

        /**
         * Close the file
         */
        virtual void close() = 0;

        /**
         * Reset the inputfile
         */
        virtual void reset() = 0;
};

// ============================================================================
// REGULAR READFILE HANDLER
// ============================================================================

class RegularReadFileHandler : public ReadFileHandler {

protected:
        FILE *fh;                               // file handler
        int bufSize;
        char *buffer;

public:
        RegularReadFileHandler() {
                bufSize = 1024;
                buffer = (char*)malloc(bufSize * sizeof(char));
        }

        /**
         * Virtual destructor for you convenience
         */
        ~RegularReadFileHandler() {
                free(buffer);
        }

        /**
         * Open a file with a given filename
         * @param filename File to open
         * @param mode Read (default) or write
         * @return True upon succes, false otherwise
         */
        bool open(const std::string& filename, ReadFileMode mode = READ);

        /**
         * Check if the input is still valid
         * @return True of false
         */
        bool good() {
                return !feof(fh);
        }

        /**
         * Read a line from file
         */
        const char *getLine() {
                size_t readSize = bufSize;
                size_t bufOffset = 0;

                buffer[0] = '\0';
                buffer[bufSize-1] = '\n'; // make the final charachter non-terminating

                do {
                        // the if-statement avoids a compiler warning
                        if (fgets(buffer + bufOffset, readSize, fh) != NULL) ;

                        // check if we're done reading the string
                        if (buffer[bufSize-1] != '\0' || buffer[bufSize-2] == '\n')
                                return buffer;

                        readSize = bufSize + 1;
                        bufOffset = bufSize - 1;
                        bufSize *= 2;

                        buffer = (char*)realloc(buffer, bufSize);
                        buffer[bufSize-1] = '\n';

                } while (good());

                return buffer;
        }

        /**
         * Write a line to file
         * @param line String to write
         */
        void writeLine(const std::string &line) {
                fputs(line.c_str(), fh);
        }

        /**
         * Write a character to file
         * @param c Character to write
         */
        void writeChar(char c) {
               fputc(c, fh);
        }

        /**
         * Peek at the next charachter in the stream
         * @return The character
         */
        char peekCharacter() {
                char c = fgetc(fh);
                ungetc(c, fh);
                return c;
        }

        /**
         * Get one character from the stream
         * @return The character
         */
        char getCharacter() {
                char c = fgetc(fh);
                return c;
        }

        /**
         * Close the file
         */
        void close();

        /**
         * Reset the inputfile
         */
        void reset();
};

// ============================================================================
// GZIPPED READFILE HANDLER
// ============================================================================

#ifdef HAVE_ZLIB

class GZipReadFileHandler : public ReadFileHandler {

protected:
        gzFile ifs;     // gzipped input file stream
        int bufSize;
        char *buffer;

public:
        /**
         * Default constructor
         */
        GZipReadFileHandler() {
                bufSize = 1024;
                buffer = (char*)malloc(bufSize * sizeof(char));
        }

        /**
         * Destructor
         */
        ~GZipReadFileHandler() {
                free(buffer);
        }

        /**
         * Open a file with a given filename
         * @param filename File to open
         * @param mode Read (default) of write
         * @return True upon succes, false otherwise
         */
        bool open(const std::string& filename, ReadFileMode mode);

        /**
         * Check if the input is still valid
         * @return True of false
         */
        bool good() {
                return !gzeof(ifs);
        }

        const char *getLine()  {
                size_t readSize = bufSize;
                size_t bufOffset = 0;

                buffer[0] = '\0';
                buffer[bufSize-1] = '\n'; // make the final charachter non-terminating

                do {
                        gzgets(ifs, buffer + bufOffset, readSize);

                        // check if we're done reading the string
                        if (buffer[bufSize-1] != '\0' || buffer[bufSize-2] == '\n')
                                return buffer;

                        readSize = bufSize + 1;
                        bufOffset = bufSize - 1;
                        bufSize *= 2;

                        buffer = (char*)realloc(buffer, bufSize);
                        buffer[bufSize-1] = '\n';

                } while (good());

                return buffer;
        }

        /**
         * Write a line to file
         * @param line String to write
         */
        void writeLine(const std::string &line) {
                gzwrite(ifs, line.c_str(), line.length());
        }

        /**
         * Write a character to file
         * @param c Character to write
         */
        void writeChar(char c) {
                gzputc(ifs, c);
        }

        /**
         * Peek at the next charachter in the stream
         * @return The character
         */
        char peekCharacter() {
                char c;
                c = gzgetc(ifs);
                gzungetc(c, ifs);
                return c;
        }

        /**
         * Get one character from the stream
         * @return The character
         */
        char getCharacter() {
                return gzgetc(ifs);
        }

        /**
         * Close the file
         */
        void close();

        /**
         * Reset the inputfile
         */
        void reset();
};

#endif

// ============================================================================
// READFILE
// ============================================================================

class ReadFile {

protected:

        ReadFileHandler *rfHandler;     // read file handler

public:
        /**
         * Default constructor
         * @param gzipped True if the file is gzipped
         */
        ReadFile(bool gzipped);

        /**
         * Destructor
         */
        virtual ~ReadFile() {
                delete rfHandler;
        }

        /**
         * Open a file
         * @param filename File to open
         * @param mode READ or WRITE
         */
        void open(const std::string& filename, ReadFileMode mode = READ) {
                rfHandler->open(filename, mode);
        }

        /**
         * Check if the input is still valid
         * @return True of false
         */
        bool good() {
                return rfHandler->good();
        }

        /**
         * Reset file to starting position
         */
        void reset() {
                rfHandler->reset();
        }

        /**
         * Close a read file
         */
        void close() {
                rfHandler->close();
        }

        /**
         * Get the next read from a file
         * @param read String containing the read (output)
         */
        virtual bool getNextRead(std::string &read) = 0;

        /**
         * Get the next record from a file
         * @param record Record to store output
         * @return true upon succesful reading
         */
        virtual bool getNextRecord(ReadRecord& record) = 0;

        /**
         * Write a record to file
         * @param record Record to write
         */
        virtual void writeRecord(const ReadRecord& record);
};

#endif
