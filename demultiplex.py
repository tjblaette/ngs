import sys

"""
T.J.Blätte
2019
"""

"""
Demultiplex FASTQ files by extracting all those reads whose
index matches a given sequence and saving these to separate
FASTQ files. At most one mismatch is allowed within each
extracted index sequence.

Args:
    - Paired FASTQ files containing actual read sequences (_R1, _R2).
    - Paired FASTQ files containing index sequences to be used
        for demultiplexing (_I1, _I2)
    - The index sequences to extract. These must be contained within
        the I1 and I2 index files, respectively. For Illumina, the i7
        index will be contained within I1 and should be provided before
        the i5 index, which will be contained within I2.

Output:
    - Paired read and index FASTQ files (_R1, _R2, _I1, _I2)
        containing all those reads from the input FASTQ files
        that could be mapped to the input index sequence.
        Output files names correspond to a concatenation of
        the input FASTQ file prefix and the index sequences
        whose reads were extracted.
"""


def extract(read_files, index_files, correct_indices):
    """
    Read in read and index FASTQ files as well as the index combination
    that is to be extracted and save all those reads with matching indices
    to separate FASTQ files. Matching indices are allowed to harbor one
    mismatch each, corresponding to a Hamming distance of at most 1.

    Args:
        read_files ([str, str]): List of FASTQ filenames containing
            the read sequences to be demultiplexed (_R1, _R2).
        index_files ([str, str]): List of FASTQ filenames containing
            the index sequences to be demultiplexed (_I1, _I2).
        correct_indices ([str, str]): List of indices to be extracted
            where correct_indices[0] has to match the index in I1 FASTQ
            and correct_indices[1] has to match the one in I2.
    """
    correct_index1, correct_index2 = correct_indices
    prnt_msg = "\nIndex combination that will be extracted: {}"
    print(prnt_msg.format(correct_indices))

    index_file1, index_file2 = index_files
    read_file1, read_file2 = read_files
    print("Processing read files: {}".format(read_files))
    print("Processing index files: {}".format(index_files))
    out_file_prefix = (
            read_file1.split("R")[0] 
            + correct_index1 
            + "-" 
            + correct_index2
            + "_")
    n_demultiplexed = 0
    n_total = 0

    with open(index_file1, "r") as (index1
            ), open(index_file2, "r") as (index2
            ), open(read_file1, "r") as (read1
            ), open(read_file2, "r") as (read2
            ), open(out_file_prefix + "R1.fastq", "w") as (out_read1
            ), open(out_file_prefix + "R2.fastq", "w") as (out_read2
            ), open(out_file_prefix + "I1.fastq", "w") as (out_index1
            ), open(out_file_prefix + "I2.fastq", "w") as (out_index2
            ):
        next_index1 = [index1.readline() for i in range(4)]
        while(next_index1[0]):
            next_index2 = [index2.readline() for i in range(4)]
            next_read1 = [read1.readline() for i in range(4)]
            next_read2 = [read2.readline() for i in range(4)]
            n_total += 1

            if (hamming_distance(next_index1[1].strip("\n"), correct_index1) <= 1 
                    and hamming_distance(next_index2[1].strip("\n"), correct_index2) <= 1):
                n_demultiplexed += 1
                for out_file, this_read in zip(
                        [out_read1, out_read2, out_index1, out_index2], 
                        [next_read1, next_read2, next_index1, next_index2]
                        ):
                    for i in range(4):
                        out_file.write(this_read[i])

            next_index1 = [index1.readline() for i in range(4)]

    prnt_msg = "---\nTotal number of reads processed: {}"
    print(prnt_msg.format(n_total))
    prnt_msg = "Number of reads successfully demultiplexed: {}"
    print(prnt_msg.format(n_demultiplexed))


def hamming_distance(seq1, seq2):
    """
    Calculate the Hamming distance between two sequences 
    of the same length. 
    The Hamming distance corresponds to the number of characters
    that differ between these two sequences.

    Args:
        seq1 (str), seq2 (str): Sequences to compare.

    Returns:
        Hamming distance (int) of seq1 and seq2.
    """
    dist = sum([char1 != char2 for char1, char2 in zip(seq1, seq2)])
    return dist



# read in commandline arguments
read_fastqs = sys.argv[1:3]
index_fastqs = sys.argv[3:5]
correct_indices = sys.argv[5:7]

# demultiplex
extract(read_fastqs, index_fastqs, correct_indices)




