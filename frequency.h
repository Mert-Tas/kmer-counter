/****************************************************************
DESCRIPTION: Class to calculate the top n frequencies of k-mers of
             size k from a given FASTQ file.
AUTHOR:      Mert TAS

******************************************************************/

#ifndef FREQUENCY_H
#define FREQUENCY_H

// Output of the program (top frequent kmers) is written to a file
#define OUTPUT_TO_FILE  1

#include <iostream>
#include <fstream>      // file stream
#include <vector>
#include <map>          // map & multimap
#include <unordered_map>
#include <stdlib.h>     // atoi

// Mask to decode bits to char
#define SHIFT_MASK 0x3	// binary 11
// The maximum kmer length that could be encoded
#define KMER_LIMIT  32

using namespace std;

// unsorted map with uint64_t encoded k-mers
typedef unordered_map<uint64_t, uint32_t> encodedMap;
typedef multimap<uint32_t, uint64_t> encodedMultiMap;
typedef vector<pair<string, uint32_t> > vecPair;

// unsorted map with no encoding
typedef unordered_map<string, uint32_t> frqMap;
typedef multimap<uint32_t, string> frqMultiMap;

enum BinaryBases
{
	SEQ_A = 0x0,	// binary: 00
	SEQ_C = 0x1,	// binary: 01
	SEQ_G = 0x2,	// binary: 10
	SEQ_T = 0x3,	// binary: 11
};

class Frequency
{
    public:
         Frequency(const char* fname, const size_t &ksize, const int &knum,
                   const char &includeLast = 'N');
        ~Frequency();

        inline uint64_t encoder(string kmer);
        inline string decoder(uint64_t code);

        int findKmersEncoded(ifstream &in);
        int findKmers(ifstream &in);

        vecPair sortKmersEncoded();
        vecPair sortKmers();

        int run();

    private:
        // Manually set buffer size for the ifstream
        // A larger buffer size (as far as the system can handle) is expected to
        // speed up file reading by limiting interaction with the hard disk drive
        // and the number of system calls.
        enum
        {
            BufferSize = 1024 //1048576 // 2^20  // 1024 // 16384
        };
        char _buffer[BufferSize];

        const char* filename;
        size_t kmer_size;
        int top_kmer_num;
        // Counts the frequency of k-mers
        encodedMap mapFrq;
        frqMap kmerFrq;

        // If user selects to view the kmers that are not in the top list
        // but have the same frequency with the last kmer in the list
        bool includeLastKmer;
        int  cntLastFrq;

        // Timer start value
        uint64_t start;

    #if (OUTPUT_TO_FILE)
        ofstream outFile;
    #endif
};

#endif // FREQUENCY_H
