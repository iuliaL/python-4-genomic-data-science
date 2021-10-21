import read_fasta
import sys
import seq_lengths

"""
In molecular biology, a reading frame is a way of dividing the DNA sequence of nucleotides into a set of consecutive,
non-overlapping triplets (or codons).
Depending on where we start, there are six possible reading frames:
three in the forward (5' to 3') direction and three in the reverse (3' to 5').
For instance, the three possible forward reading frames for the sequence AGGTGACACCGCAAGCCTTATATTAGC are:

AGG TGA CAC CGC AAG CCT TAT ATT AGC

A GGT GAC ACC GCA AGC CTT ATA TTA GC

AG GTG ACA CCG CAA GCC TTA TAT TAG C

These are called reading frames 1, 2, and 3 respectively.

An open reading frame (ORF) is the part of a reading frame that has the potential to encode a protein.
Has the potential underlined, does not mean the aminoacid sequence will translate to proteins necessarily.
It starts with a start codon (ATG), and ends with a stop codon (TAA, TAG or TGA).
For instance, ATGAAATAG is an ORF of length 9, which in this case would translate to one aminoacid only. AAA encodes for Lysine.
Given an input reading frame on the forward strand (1, 2, or 3) your program should be able to identify
all ORFs present in each sequence of the FASTA file, and answer the following questions:"""


def seq_as_reading_frames(seq, start=0):
    """Start can be index 0, 1, 2. Default is 0.
        Returns list of 3-mers"""
    result = []
    for i in range(start, len(seq), 3):
        curr_3mer = seq[i:i + 3]
        if len(curr_3mer) == 3:  # account for the possible 1 or 2mers at the end of the sequence if start is > 0
            result.append(curr_3mer)
    return result


def find_orfs_in_seq(seq):
    """ The input is a list of 3 lists corresponding to the reading frames corresponding to a sequence for each start 0,1,2
        Returns a list the ORFs for a sequence """
    reading_frames = [
        seq_as_reading_frames(seq),
        seq_as_reading_frames(seq, 1),
        seq_as_reading_frames(seq, 2)
    ]
    orfs_found = []
    for i in range(0, 3):
        # print('Start position *****************************', i)
        orf = []
        for frame in reading_frames[i]:
            if is_START_codon(frame):
                # print(frame, "is start codon and orf so far", orf)
                if (not orf):  # if this frame would be the first in the ORF
                    orf.append(frame)
                    continue
                else:  # reset
                    orf = []
            if is_STOP_codon(frame):
                if (not orf):  # if the stop codon frame is the first in the ORF
                    continue
                elif len(orf) < 2:  # if there are less than 2 frames before
                    orf = []
                    continue
                else:
                    # print(frame, "is STOP codon and orf so far", orf)
                    orf.append(frame)
                    # reset orf
                    orfs_found.append("".join(orf))
                    orf = []
            else:  # if neither START nor STOP codons
                # print(frame, "not start not stop codon and orf so far", orf)
                orf.append(frame)
    # do not forget to filter out those who dont start by the START codon and end with any of the STOP codons
    filtered = [o for o in orfs_found if o.startswith('ATG') and (len(o) >= 9) and is_STOP_codon(o[len(o) - 3: len(o)])]
    return filtered


def longest_orf_for_seq(orfs):
    return max(orfs, key=len)


def is_START_codon(reading_frame):
    """Returns True is ATG"""
    return reading_frame == "ATG"


def is_STOP_codon(reading_frame):
    """Returns True if is TAA, TAG or TGA"""
    return reading_frame in ('TAA', 'TAG', 'TGA')


def longest_ORF_in_file(f):
    """ What is the identifier of the sequence containing the longest ORF?
        Returns List of Tuple(identifier, length) because maybe i would have 2 longest ORF with the same length"""
    records = read_fasta.read(f)
    records_of_longest_orfs = {identifier: longest_orf_for_seq(find_orfs_in_seq(seq)) for identifier, seq in records.items()}
    records_of_orfs_lengths = seq_lengths.compute_lengths(records_of_longest_orfs)
    longest_Orf_in_file = max(records_of_orfs_lengths.values())
    return [(identifier,  longest_Orf_in_file) for identifier in records_of_orfs_lengths if records_of_orfs_lengths[identifier] == longest_Orf_in_file]



if __name__ == "__main__":
    print(longest_ORF_in_file(sys.argv[1]))
