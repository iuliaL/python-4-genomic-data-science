import read_fasta
import sys


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


def find_orfs_in_reading_frames(rfs):
    """The input is a list of 3 lists corresponding to the reading frames for each start 0,1,2 """
    orfs_found = []
    for i in range(0, 3):
        # print('Start position *****************************', i)
        orf = []
        for frame in rfs[i]:
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
    filtered = [o for o in orfs_found if (o.startswith('ATG') and is_STOP_codon(o[len(o) - 3: len(o)]))]
    return filtered


def is_START_codon(reading_frame):
    """Returns True is ATG"""
    return reading_frame == "ATG"


def is_STOP_codon(reading_frame):
    """Returns True if is TAA, TAG or TGA"""
    return reading_frame in ('TAA', 'TAG', 'TGA')


def longest_ORF_in_file(f):
    """ What is the identifier of the sequence containing the longest ORF?
        Returns Tuple(identifier, length) """
    records = read_fasta.read(f)
    for identifier, seq in records.items():
        reading_frames = [
            seq_as_reading_frames(seq),
            seq_as_reading_frames(seq, 1),
            seq_as_reading_frames(seq, 2)
        ]
        orfs_for_seq = find_orfs_in_reading_frames(reading_frames)
        print("Sequence identifier:", identifier, orfs_for_seq)

    return None
    # return (identifier, length)


# def longest_ORF_length_in_file(f):
#     """ What is the length of the longest ORF in the file?
#         Returns integer """
#     return longest_ORF_in_file(f)[1]

# def longest_ORF_id_in_file(f):
#     """ What is the length of the longest ORF in the file?
#         Returns string """
#     return longest_ORF_in_file(f)[0]


if __name__ == "__main__":
    print(longest_ORF_in_file(sys.argv[1]))
