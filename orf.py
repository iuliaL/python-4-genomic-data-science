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


def is_START_codon(reading_frame):
    """Returns True if ATG"""
    return reading_frame == "ATG"


def is_STOP_codon(reading_frame):
    """Returns True if is TAA, TAG or TGA"""
    return reading_frame in ('TAA', 'TAG', 'TGA')


def startstop_codon(zero_index_pos, seq):
    for i in range(zero_index_pos, len(seq), 3):
        codon1 = seq[i:i + 3]
        if is_START_codon(codon1):
            position1 = i
            for j in range(position1, len(seq), 3):
                codon2 = seq[j:j + 3]
                if is_STOP_codon(codon2):
                    position2 = j
                    yield (zero_index_pos + 1, position2 - position1 + 3, seq[position1:position2 + 3])
                    break


def longest_orfs_in_seq_per_pos(seq):
    """ The input is the sequence
        Returns a dict of reading_frame_pos : Tuple with (ORF_length, ORF) for a sequence """
    orfs = {1: [], 2: [], 3: []}
    for p in range(0, 3):
        for pos, orflen, orf in startstop_codon(p, seq):
            orfs[pos].append((orflen, orf))
    longest_orfs = {}
    count_total_orfs_found = 0
    for pos, orf_list in orfs.items():
        count_total_orfs_found += len(orf_list)
        longest_orfs[pos] = max(orf_list, default=(0, ''), key=lambda item: item[0])
    # print('Total orfs found to check against NCBI ORF Finder', count_total_orfs_found)
    return longest_orfs


def longest_orf_in_seq(seq):
    longest_orfs_per_position = longest_orfs_in_seq_per_pos(seq)
    longest_length = 0
    orf_found = (1, 0, '')
    for pos, (length, orf) in longest_orfs_per_position.items():
        if length > longest_length:
            orf_found = pos, length, orf
    return orf_found


def records_with_longest_orfs(f):
    """ What is the identifier of the sequence containing the longest ORF throughout the sequences in the whole fasta file?
        What is the starting position of the longest ORF in the sequence that contains it ? Starting position should be 1/2/3
        Returns dict(identifier, seq, reading frame starting position, longest ORF length and the ORF)"""
    records = read_fasta.read(f)
    records_of_longest_orfs = {}
    for id, seq in records.items():
        records_of_longest_orfs[id] = {
            "seq": seq,
            "orfs": longest_orf_in_seq(seq)
        }
    return records_of_longest_orfs

def records_with_longest_orfs_at_pos (f,pos):
    """ What is the identifier of the sequence containing the longest ORF throughout the sequences in the whole fasta file?
        What is the starting position of the longest ORF in the sequence that contains it ? Starting position should be 1/2/3
        Returns dict(identifier, seq, reading frame starting position, longest ORF length and the ORF)"""
    records = read_fasta.read(f)
    records_of_longest_orfs = {}
    for id, seq in records.items():
        records_of_longest_orfs[id] = {
            "seq": seq,
            "orfs": longest_orfs_in_seq_per_pos(seq)[pos]
        }
    return records_of_longest_orfs


def longest_ORF_in_file(f):
    records = records_with_longest_orfs(f)
    longest_length = 0
    longest = None
    for id, info in records.items():
        if info['orfs'][1] > longest_length:
            longest = info['orfs']
    return id, info['seq'], longest


def longest_ORF_in_file_at_pos(f, pos):
    records = records_with_longest_orfs_at_pos(f, pos)
    longest_length = 0
    longest = None
    for id, info in records.items():
        print(info['orfs'])
        if info['orfs'][0] > longest_length:
            longest = info['orfs']
    return id, info['seq'], pos, longest




def longest_ORF_for_given_id(f, identifier):
    """For a given sequence identifier, what is the longest ORF contained in the sequence represented by that identifier?"""
    records = read_fasta.read(f)
    desired_seq = records[identifier]
    return longest_orf_in_seq(desired_seq)


if __name__ == "__main__":
    # Longest ORF in whole file

    # l = longest_ORF_in_file(sys.argv[1])
    # id = l[0]
    # seq = l[1]
    # start_pos = l[2][0]
    # l_length = l[2][1]
    # orf = l[2][2]
    # index = seq.index(orf)
    # print(
    #     "\nLongest ORF in the whole file has id", id, 
    #     "has sequence", seq,  "has reading frame start position",
    #     start_pos,
    #     "and length of",
    #     l_length,
    #     "and is",
    #     orf, "and the orf starts inside the sequence at pos: index + 1", index + 1)

    # Longest ORF in whole file at the given start position
    l = longest_ORF_in_file_at_pos(sys.argv[1], 2)
    id = l[0]
    seq = l[1]
    start_pos = l[2]
    l_length = l[3][0]
    orf = l[3][1]
    print(
        "\nLongest ORF in the whole file at start position", start_pos, " has id",  id, 
        "has sequence", seq,  
        "and length of",
        l_length,
        "and is",
        orf)




   
    # Longest ORF in a given identifier i.e. gi|142022655|gb|EQ086233.1|16

    # given_id = 'gi|142022655|gb|EQ086233.1|16'
    # print('Longest ORF for ', given_id, 'is', longest_ORF_for_given_id(sys.argv[1], given_id))
