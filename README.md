# kmer-counter
Finds the k-mers of each DNA sequence in a given FASTQ file and counts their frequencies.  It then selects the most frequent n k-mers and displays them.

# Executing the program
The program accepts three (compulsory *) inputs via command line:
(1) filename, (2) k-mer size, (3) top frequency count

and it also accepts an optional input that is used to list the k-mers that are not in 
the top frequency list but have the same frequency as the last k-mer of the top frequency
list (so they deserve to be in the list). By selecting this last input they can also be listed.

(4) include last k-mers = (Y)es / (N)o  To include the last kmers to the top frequency list.


"kmercnt.exe <FILE_NAME> <KMER_SIZE> <TOP_FREQUENCY_CNT> <Y/N>"

# Performance

This program is tested on a laptop with a 3GB of RAM. 

The approximate time to count and list the top frequent kmers are as follows:
- fastq file Size:  45,1 Mb   Number of lines in file: 799,440     DNA sequence length: 90   elapsed time (appr.):  1 min 50 sec. 
- fastq file Size:  1,53 Gb   Number of lines in file: 40 million  DNA sequence length: 90   elapsed time (appr.):  27 min.

# Other notes

- This program stores the output in a vector and also prints it to the screen. To print the 
output to a file please set the value of the define as following:

	#define OUTPUT_TO_FILE  1 => Write the output to file
	#define OUTPUT_TO_FILE  0 => Do not write the output to file
 
- To measure the operation time of the code snippets a function by Andreas Bonini is used. <Thanks!>
stackoverflow.com/questions/1861294/how-to-calculate-execution-time-of-a-code-snippet-in-c

- This program is tested with the files from the 1000 Genomes Project.
http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG01595/sequence_read/
