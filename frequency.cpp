/**
  Frequency class implementation
**/

#include "frequency.h"
#include "evaltime.h"

Frequency::Frequency(const char* fname, const size_t &ksize, const int &knum,
                     const char &includeLast):
                     filename(fname), kmer_size(ksize), top_kmer_num(knum), cntLastFrq(0)
{
    includeLastKmer = (toupper(includeLast) == 'Y') ? true : false;

#if (OUTPUT_TO_FILE)
    outFile.open("out.txt");
#endif
}

Frequency::~Frequency()
{
#if (OUTPUT_TO_FILE)
    outFile.close();
#endif
}

vecPair Frequency::sortKmers()
{
    start = GetTimeMs64();
    frqMultiMap kmer_swapped;

    try
    {
        // Swaps the keys and values of kmer map, sorting the elements according to their values. O(nlogn)
        for (auto it = kmerFrq.begin(); it != kmerFrq.end(); ++it)
        {
            kmer_swapped.insert(pair<uint32_t, string>(it -> second, it -> first));
        }
        // Clear the containers that are not used anymore.
        kmerFrq.clear();
    }
    catch (const std::bad_alloc& e)
    {
        cerr << "Swapping key-value pairs. Allocation failed: " << e.what() << "\n";
        kmerFrq.clear();
        return vecPair();
    }

    auto rev_end = kmer_swapped.rbegin();
    try
    {
        // Takes the most frequent n elements from the end of the map which is sorted in ascending order.
        advance(rev_end, top_kmer_num);
    }
    catch (const std::bad_alloc& e)
    {
        cerr << "Getting most frequent elements. Allocation failed: " << e.what() << "\n";
        return vecPair();
    }

    vecPair kmerVector;
    try
    {
        // Copies the top n kmers from the end of the multimap to a vector in descending order.
        for (auto rev = kmer_swapped.rbegin(); rev != rev_end; ++rev)
        {
            kmerVector.push_back(make_pair(rev->second, rev->first));
        }
    }
    catch (const std::bad_alloc& e)
    {
        cerr << "Copying to vector. Allocation failed: " << e.what() << "\n";
        kmerVector.clear();
        return vecPair();
    }

    // If user selected to view the kmers that are not on the top list but have
    // enough frequency to be in the top list.
    if (includeLastKmer)
    {
        for (auto rev = rev_end; rev != kmer_swapped.rend(); ++rev, ++cntLastFrq)
        {
            // Add these kmers to the top frequency list
            if (rev->first == kmerVector.back().second)
            {
                kmerVector.push_back(make_pair(rev->second, rev->first));
            }
            else
            {
                break;
            }
        }
    }

    // Clear the containers that are not used anymore.
    kmer_swapped.clear();
    cout << "kmer sort time: " << GetTimeMs64() - start << " ms. \n";

    // Write the result to the screen and file.
    cout << "Top " << top_kmer_num <<  " frequent k-mers of size " << kmer_size << " and their frequencies:\n\n";
#if (OUTPUT_TO_FILE)
    outFile << "\nTop " << top_kmer_num <<  " frequent k-mers of size " << kmer_size << " and their frequencies:\n\n";
#endif

    for (size_t i = 0; i < kmerVector.size(); ++i)
    {
        cout.fill(' ');
        cout.width(2);
        cout << (i+1) << ".  " << kmerVector[i].first << " : " << kmerVector[i].second <<  "\n";
#if (OUTPUT_TO_FILE)
        // Print to file
        outFile.fill(' ');
        outFile.width(2);
        outFile << (i+1) << ".  " << kmerVector[i].first << " : " << kmerVector[i].second <<  "\n";
#endif
    }

    if (includeLastKmer)
    {
        cout << "\n* " << cntLastFrq << " more kmer(s) are also listed with the same frequencies of the last kmer(s)\n";
#if (OUTPUT_TO_FILE)
        outFile << "\n* " << cntLastFrq << " more kmer(s) are also listed with the same frequencies of the last kmer(s)\n";
#endif
    }

    return kmerVector;
}

vecPair Frequency::sortKmersEncoded()
{
    start = GetTimeMs64();
    encodedMultiMap kmer_swapped;

    try
    {
        // Swaps the keys and values of kmer map, sorting the elements according to their values. O(nlogn)
        for (auto it = mapFrq.begin(); it != mapFrq.end(); ++it)
        {
            kmer_swapped.insert(pair<uint32_t, uint64_t>(it -> second, it -> first));
        }
        // Clear the containers that are not used anymore.
        mapFrq.clear();
    }
    catch (const std::bad_alloc& e)
    {
        cerr << "Swapping key-value pairs. Allocation failed: " << e.what() << "\n";
        mapFrq.clear();
        return vecPair();
    }

    auto rev_end = kmer_swapped.rbegin();
    try
    {
        // Takes the most frequent n elements from the end of the map which is sorted in ascending order.
        advance(rev_end, top_kmer_num);
    }
    catch (const std::bad_alloc& e)
    {
        cerr << "Getting most frequent elements. Allocation failed: " << e.what() << "\n";
        return vecPair();
    }

    vecPair kmerVector;
    try
    {
        // Copies the top n kmers from the end of the multimap to a vector in descending order.
        for (auto rev = kmer_swapped.rbegin(); rev != rev_end; ++rev)
        {
            // The decoded kmers are encoded while copying them.
            kmerVector.push_back(make_pair(decoder(rev->second), rev->first));
        }
    }
    catch (const std::bad_alloc& e)
    {
        cerr << "Copying to vector. Allocation failed: " << e.what() << "\n";
        kmerVector.clear();
        return vecPair();
    }

    // If user selected to view the kmers that are not on the top list but have
    // enough frequency to be in the top list.
    if (includeLastKmer)
    {
        for (auto rev = rev_end; rev != kmer_swapped.rend(); ++rev, ++cntLastFrq)
        {
            // Decode these kmers and add them to the top frequency list
            if (rev->first == kmerVector.back().second)
            {
                kmerVector.push_back(make_pair(decoder(rev->second), rev->first));
            }
            else
            {
                break;
            }
        }
    }

    // Clear the containers that are not used anymore.
    kmer_swapped.clear();
    cout << "kmer sort time: " << GetTimeMs64() - start << " ms. \n";

    // Write the result to the screen and file.
    cout << "Top " << top_kmer_num <<  " frequent k-mers of size " << kmer_size << " and their frequencies:\n\n";
#if (OUTPUT_TO_FILE)
    outFile << "\nTop " << top_kmer_num <<  " frequent k-mers of size " << kmer_size << " and their frequencies:\n\n";
#endif

    for (size_t i = 0; i < kmerVector.size(); ++i)
    {
        cout.fill(' ');
        cout.width(2);
        cout << (i+1) << ".  " << kmerVector[i].first << " : " << kmerVector[i].second <<  "\n";
#if (OUTPUT_TO_FILE)
        // Print to file
        outFile.fill(' ');
        outFile.width(2);
        outFile << (i+1) << ".  " << kmerVector[i].first << " : " << kmerVector[i].second <<  "\n";
#endif
    }

    if (includeLastKmer)
    {
        cout << "\n* " << cntLastFrq << " more kmer(s) are also listed with the same frequencies of the last kmer(s)\n";
#if (OUTPUT_TO_FILE)
        outFile << "\n* " << cntLastFrq << " more kmer(s) are also listed with the same frequencies of the last kmer(s)\n";
#endif
    }

    return kmerVector;
}

int Frequency::findKmersEncoded(ifstream &in)
{
    start = GetTimeMs64();
    string line;
    // Get the second line which is the first dna sequence
    getline(in, line);
    getline(in, line);
    // In this case, line counter starts from the second line
    uintmax_t cnt_lines = 2;

    // Checks if the kmer size is smaller than the sequence length.
    if (kmer_size > line.length())
    {
        cerr << "Error! k-mer size is larger than the sequence length. Exiting ...\n";
        exit(0);
    }

    // Calculate the kmer number in a single dna sequence
    int kmer_num = line.length() - kmer_size + 1;
    // Finds the kmers of the first dna sequence and count their frequencies  O(kmer_num + kmer_size)
    for (int i = 0; i < kmer_num; ++i)
    {
        ++mapFrq[encoder(line.substr(i, kmer_size))];
    }

    while (getline(in, line))
    {
        // Get the second line (dna sequence) of each  nucleotide
        if (++cnt_lines % 4 == 2)
        {
            // Find the kmers of the remaining dna sequences and count their frequencies
            for (int i = 0; i < kmer_num; ++i)
            {
                ++mapFrq[encoder(line.substr(i, kmer_size))];
            }
        }
    }
    cout << "kmer find time: " << GetTimeMs64() - start << " ms. \n";
    return 1;
}

int Frequency::findKmers(ifstream &in)
{
    start = GetTimeMs64();
    string line;
    // Get the second line which is the first dna sequence
    getline(in, line);
    getline(in, line);
    // In this case, line counter starts from the second line
    uintmax_t cnt_lines = 2;

    // Checks if the kmer size is smaller than the sequence length.
    if (kmer_size > line.length())
    {
        cerr << "Error! k-mer size is larger than the sequence length. Exiting ...\n";
        exit(0);
    }
    cout << "k-mer size is larger than " << KMER_LIMIT << ". Encoding will not be used.\n";

    // Calculate the kmer number in a single dna sequence
    int kmer_num = line.length() - kmer_size + 1;
    // Finds the kmers of the first dna sequence and count their frequencies  O(kmer_num + kmer_size)
    for (int i = 0; i < kmer_num; ++i)
    {
        ++kmerFrq[line.substr(i, kmer_size)];
    }

    while (getline(in, line))
    {
        // Get the second line (dna sequence) of each  nucleotide
        if (++cnt_lines % 4 == 2)
        {
            // Find the kmers of the remaining dna sequences and count their frequencies
            for (int i = 0; i < kmer_num; ++i)
            {
                ++kmerFrq[line.substr(i, kmer_size)];
            }
        }
    }
    cout << "kmer find time: " << GetTimeMs64() - start << " ms. \n";
    return 1;
}

int Frequency::run()
{
    ifstream fastq_file(filename);
    // Set the buffer size manually
    fastq_file.rdbuf()->pubsetbuf(_buffer, BufferSize);

    if (fastq_file.is_open())
    {
        cout << "Input file " << filename << " opened for extracting k-mers.\n";

        if (kmer_size <= KMER_LIMIT)
            findKmersEncoded(fastq_file);
        else
            findKmers(fastq_file);

        fastq_file.close();
    }
    else
    {
        cerr << "Unable to open file" << endl;
        return 0;
    }

    // Print the top kmers
    if (kmer_size <= KMER_LIMIT)
        sortKmersEncoded();
    else
        sortKmers();

    return 1;
}

// Encoder function that returns the binary representation of
// the k-mers in unsigned integer format
inline uint64_t Frequency::encoder(string kmer)
{
    uint64_t code = 0;
    for (size_t i = 0; i < kmer.size(); ++i)
    {
        code *= 4;
        switch(kmer[i])
        {
            case 'A':
                code += SEQ_A;
                break;
            case 'C':
                code += SEQ_C;
                break;
            case 'G':
                code += SEQ_G;
                break;
            case 'T':
                code += SEQ_T;
                break;
        }
    }
    return code;
}

// Decoder function that returns the string represenations of
// the binary encoded k-mers
inline string Frequency::decoder(uint64_t code)
{
    string seq;
    for (size_t i = 0; i < kmer_size; ++i)
    {
        switch(code & SHIFT_MASK)
        {
            case 0:
                seq = "A" + seq;
                break;
            case 1:
                seq = "C" + seq;
                break;
            case 2:
                seq = "G" + seq;
                break;
            case 3:
                seq = "T" + seq;
                break;
        }
        code >>= 2;
    }
    return seq;
}
