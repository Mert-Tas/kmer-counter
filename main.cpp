#include "frequency.h"

int main(int argc, char* argv[])
{
    if (argc < 4)
    {
        cerr << "Usage: " << argv[0] << " <FILENAME> <K-MER SIZE> <NUMBER OF K-MERS>\n";
        cerr << "(*optional) <INCLUDE LAST K-MERS>\n\n";
        cerr << "(*) NOTE: Last argument is optional. If you want to list the kmers\n";
        cerr << "that are not in the top n frequent k-mers, but have the same frequency\n";
        cerr << "of the last kmers in the top list select this option.\n";
        cerr << "\t(Y)es : Include last k-mers\n\t(N)o : Do not include last k-mers\n";
        cerr << "Please note that even if this option is selected there might be no\n";
        cerr << "more k-mers to be listed.\n";
        return 0;
    }

    // Basic input controls
    if (atoi(argv[2]) < 1)
    {
        cerr << "Invalid k-mer size.\n";
        return 0;
    }
    if (atoi(argv[3]) < 1)
    {
        cerr << "Invalid number of top k-mers.\n";
        return 0;
    }

    Frequency freq(argv[1], atoi(argv[2]), atoi(argv[3]), (argv[4] == NULL) ? 'N' : argv[4][0]);
    freq.run();

    return 1;
}
