import numpy as np
import argparse
import pprint
import os
import copy

"""
T.J.Blätte
2018
"""

"""
To extract molecular barcodes sequenced as part of the main reads
into separate FASTQ files.

Args:
    fastq1: FASTQ file of forward reads (REQUIRED)
    fastq2: FASTQ file of reverse reads (REQUIRED)
    sampleID: sample ID used as output file prefix (REQUIRED)
    primer1: Gene-specific primer sequence of forward reads (REQUIRED)
    primer2: Gene-specific primer sequence of reverse reads (REQUIRED)

Output:
    sampleID_indexed_R1.fastq: fastq1 minus the extracted sequence
    sampleID_indexed_R2.fastq: fastq2 minus the extracted sequence
    sampleID_indexed_I1.fastq: the sequences extracted from fastq1
    sampleID_indexed_I2.fastq: the sequences extracted from fastq2
    sampleID_indexed_I.fastq: concatenation of extracted sequences
"""

class Read(object):
    """
    Sequencing read.
    """
    def __init__(
                self,
                seq,
                name=None,
                desc="+",
                index=None,
                sense=1,
                bqs=None,
                index_bqs=None,
                counts=1,
                al_score=None,
                al_seq=None,
                al_ref=None,
                al_file=None):
        """
        Initialize Read.

        Args:
            seq (str): Base pair sequence of the read.
            name (str): Read ID = 1st line of FASTQ file entry.
            desc (str): Read description = 3rd line of FASTQ file entry.
            index (int): Read index to identify paired reads
                for paired-end sequencing
                --> will have the same index
            sense (int): 1 for forward reads, -1 for reverse.
            bqs (str): Base quality scores of the read.
            index_bqs (str): Base quality scores of the
                corresponding index read.
            counts (int): Total number of reads with these
                exact attributes.
            al_score (int): Score of the read-to-reference
                alignment.
            al_seq (str): Read sequence aligned to the reference.
            al_ref (str): Reference sequence aligned to the read.
            al_file (str): Name of the file containing the
                read-to-reference alignment.
        """
        self.seq = seq
        self.name = name
        self.desc = desc
        self.index = index
        self.bqs = bqs
        self.index_bqs = index_bqs
        self.length = len(seq)
        self.counts = counts
        self.sense = sense
        self.al_score = al_score
        self.al_seq = al_seq
        self.al_ref = al_ref
        self.al_file = al_file
        assert self.counts > 0
        if self.bqs is not None:
            assert len(self.seq) == len(self.bqs)

    def print(self):
        """
        Pretty print Read.
        """
        pprint.pprint(vars(self))

    def get_barcode(self, primer):
        """
        Extract molecular barcode from Read,
        change Read accordingly in-place.

        Returns:
            Read() object for barcode.
        """
        substr_index = self.seq.find(primer)
        if substr_index != -1:
            barcode = self.seq[:substr_index]
            barcode_bqs = self.bqs[:substr_index]
            self.seq = self.seq[substr_index:]
            self.bqs = self.bqs[substr_index:]
            self.length = len(self.seq)
            return Read(seq=barcode, bqs=barcode_bqs, name=self.name, desc=self.desc, index=self.index, sense=self.sense)



def read_fastq(fastq_file):
    """
    Read sequence fastq file and extract sequences and BQS.

    Args:
        fastq_file: Name of the fastq file to read, R1 or R2.

    Returns:
        List of Read() objects.
    """
    reads = []
    read_index = 0
    try:
        with open(fastq_file,'r') as f:
            line = f.readline().rstrip('\n')
            while line:
                read_id = line
                read_seq = f.readline().rstrip('\n')
                read_desc = f.readline().rstrip('\n')
                read_bqs = f.readline().rstrip('\n')
                assert len(read_seq) == len(read_bqs)
                reads.append(Read(seq=read_seq, name=read_id, desc=read_desc, index=read_index, bqs=read_bqs, sense=1))
                line = f.readline().rstrip('\n')
                read_index += 1
    # catch missing file or permissions
    except IOError as e:
        print("---\nCould not read fastq file {}!\n---".format(fastq_file))
    return reads


def write_fastq(reads, fastq_file):
    """
    Write list of reads of FASTQ file.

    Args:
        reads (list): List of reads to write.
        fastq_file (str): Filename of file to write to.
    """
    try:
        with open(fastq_file,'w') as f:
            for read in reads:
                f.write(read.name + "\n")
                f.write(read.seq + "\n")
                f.write(read.desc + "\n")
                f.write(read.bqs + "\n")
    # catch missing file or permissions
    except IOError as e:
        print("---\nCould not write fastq file {}!\n---".format(fastq_file))




if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("fastq1", help="FASTQ file of forward reads (REQUIRED)")
    parser.add_argument("fastq2", help="FASTQ file of reverse reads (REQUIRED)")
    parser.add_argument("sampleID", help="sample ID used as output file prefix (REQUIRED)")
    parser.add_argument("primer1", help="Gene-specific primer sequence of forward reads (REQUIRED)")
    parser.add_argument("primer2", help="Gene-specific primer sequence of reverse reads (REQUIRED)")
    cmd_args = parser.parse_args()

    R1 = cmd_args.fastq1
    R2 = cmd_args.fastq2
    SAMPLE = cmd_args.sampleID
    PRIMER1 = cmd_args.primer1
    PRIMER2 = cmd_args.primer2

    ##################################
    ## Get Reads

    reads = read_fastq(R1)
    reverse_reads = read_fastq(R2)

    for read in reverse_reads:
        read.sense = -1

    print("Total read pairs: {}".format(len(reads)))


    ##################################
    ## Get Barcodes

    barcodes = [read.get_barcode(PRIMER1) for read in reads]
    reverse_barcodes = [read.get_barcode(PRIMER2) for read in reverse_reads]

    barcodes = [barcode for barcode in barcodes if barcode is not None]
    reverse_barcodes = [barcode for barcode in reverse_barcodes if barcode is not None]

    print("Forward reads with primer: {}".format(len(barcodes)))
    print("Reverse reads with primer: {}".format(len(reverse_barcodes)))


    ##################################
    ## Filter Barcodes

    barcodes = [barcode for barcode in barcodes if barcode.length == 16]
    reverse_barcodes = [barcode for barcode in reverse_barcodes if barcode.length == 16]

    print("Forward reads with primer and 16 bp barcode: {}".format(len(barcodes)))
    print("Reverse reads with primer and 16 bp barcode: {}".format(len(reverse_barcodes)))


    ##################################
    ## COLLECT READ INDEXES FOR WHICH EITHER R1 OR R2 PASS BARCODE AND PRIMER FILTERS
    passing_forward_indices = [barcode.index for barcode in barcodes]
    passing_reverse_indices = [barcode.index for barcode in reverse_barcodes]
    passing_indices = set(passing_forward_indices).intersection(set(passing_reverse_indices))
    print("Read pairs where forward and reverse read contain primer and 16 bp barcode: {}".format(len(passing_indices)))

    passing_reads = [read for read in reads if read.index in passing_indices]
    passing_reverse_reads = [read for read in reverse_reads if read.index in passing_indices]
    passing_barcodes = [barcode for barcode in barcodes if barcode.index in passing_indices]
    passing_reverse_barcodes = [barcode for barcode in reverse_barcodes if barcode.index in passing_indices]


    passing_full_barcodes = []
    for pb, prb in zip(passing_barcodes, passing_reverse_barcodes):
        full_barcode = copy.deepcopy(pb)
        full_barcode.seq = pb.seq + prb.seq
        full_barcode.bqs = pb.bqs + prb.bqs
        full_barcode.length = len(full_barcode.seq)
        passing_full_barcodes.append(full_barcode)


    ##################################
    ## FOR READS AND BARCODES WITH SUCH AN INDEX, WRITE FASTQ FILES

    write_fastq(passing_reads, SAMPLE + "_indexed_R1.fastq")
    write_fastq(passing_reverse_reads, SAMPLE + "_indexed_R2.fastq")
    write_fastq(passing_barcodes, SAMPLE + "_indexed_I1.fastq")
    write_fastq(passing_reverse_barcodes, SAMPLE + "_indexed_I2.fastq")
    write_fastq(passing_full_barcodes, SAMPLE + "_indexed_I.fastq")

