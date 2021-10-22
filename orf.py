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
        Returns a list of Tuples with (start, ORF) for a sequence """
    reading_frames = [
        seq_as_reading_frames(seq),
        seq_as_reading_frames(seq, 1),
        seq_as_reading_frames(seq, 2)
    ]
    orfs_found = []
    for i in range(0, 3):
        orf = []
        for frame in reading_frames[i]:
            if is_START_codon(frame):
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
                    orf.append(frame)
                    # reset orf
                    orfs_found.append((i, "".join(orf)))
                    orf = []
            else:  # if neither START nor STOP codons
                orf.append(frame)
    # do not forget to filter out those who dont start by the START codon and end with any of the STOP codons
    filtered = [o for o in orfs_found if o[1].startswith('ATG') and (len(o[1]) >= 9) and is_STOP_codon(o[1][len(o[1]) - 3: len(o[1])])]
    return filtered


def longest_orf_for_seq(orfs):
    return max(orfs, key=lambda item: len(item[1]))


def is_START_codon(reading_frame):
    """Returns True is ATG"""
    return reading_frame == "ATG"


def is_STOP_codon(reading_frame):
    """Returns True if is TAA, TAG or TGA"""
    return reading_frame in ('TAA', 'TAG', 'TGA')


def longest_ORF_in_file(f):
    """ What is the identifier of the sequence containing the longest ORF throughout the sequences in the whole fasta file?
        What is the starting position of the longest ORF in the sequence that contains it ? Starting position should be 1/2/3
        Returns dict(identifier, position, longest orf found and its length)"""
    records = read_fasta.read(f)
    records_of_longest_orfs = build_id_pos_longest_orfs(records)
    records_of_orfs_lengths = {identifier: (pos, len(orf)) for (identifier, (pos, orf)) in records_of_longest_orfs.items()}
    longest = max(records_of_orfs_lengths.values(), key=lambda item: item[1])
    identifier = [(id, longest) for id in records_of_orfs_lengths if records_of_orfs_lengths[id] == longest][0][0]
    return {"id": identifier,
            "start position": records_of_longest_orfs[identifier][0] + 1,  # !!!! atention to add 1 here cause i used the zero index
            "longest_ORF": records_of_longest_orfs[identifier][1],
            "longest_ORF_length": len(records_of_longest_orfs[identifier][1])}


def longest_ORF_in_file_at_pos(position, f):
    """ What is the identifier of the sequence containing the longest ORF throughout the sequences in the whole fasta file?
        What is the starting position of the longest ORF in the sequence that contains it ? Starting position should be 1/2/3
        Returns dict(identifier, position, longest orf found and its length)"""
    records = read_fasta.read(f)
    records_of_longest_orfs = build_id_pos_longest_orfs(records)
    records_of_orfs_lengths = {identifier: (pos, len(orf)) for (identifier, (pos, orf)) in records_of_longest_orfs.items() if pos == position - 1}
    longest = max(records_of_orfs_lengths.values(), key=lambda item: item[1])
    identifier = [(id, longest) for id in records_of_orfs_lengths if records_of_orfs_lengths[id] == longest][0][0]
    return {"id": identifier,
            "start position": records_of_longest_orfs[identifier][0] + 1,  # !!!! atention to add 1 here cause i used the zero index
            "longest_ORF": records_of_longest_orfs[identifier][1],
            "longest_ORF_length": len(records_of_longest_orfs[identifier][1])}


def build_id_pos_longest_orfs(records):
    """ Returns dict of id: (pos, longest_orf_found) """
    result = {}
    for identifier, seq, in records.items():
        orfs = find_orfs_in_seq(seq)
        # print("Sequence", identifier, "has ", len(orfs), "orfs")

        if orfs:
            result[identifier] = longest_orf_for_seq(orfs)
    return result


def longest_ORF_length_in_file(f):
    """What is the length of the longest ORF in the file?
        Returns integer """
    return longest_ORF_in_file(f)['longest_ORF_length']


def longest_ORF_in_given_seq(f, identifier):
    """For a given sequence identifier, what is the longest ORF contained in the sequence represented by that identifier?"""
    records = read_fasta.read(f)
    records_of_longest_orfs = build_id_pos_longest_orfs(records)
    return records_of_longest_orfs[identifier]  # !!!! atention to add 1 to get the position cause i used the zero index


if __name__ == "__main__":
    # Longest ORF in whole file

    # longest_orf_at_pos = longest_ORF_in_file(sys.argv[1])
    # identifier = longest_orf_at_pos['id']
    # start_pos = longest_orf_at_pos['start position']
    # l_length = longest_orf_at_pos['longest_ORF_length']
    # print("Longest ORF in file is at identifier", identifier, "with reading frame start position", start_pos, "and length of", l_length)

    # Longest ORF in whole file at the given start position

    # longest_orf_at_pos = longest_ORF_in_file_at_pos(2, sys.argv[1])
    # identifier = longest_orf_at_pos['id']
    # start_pos = longest_orf_at_pos['start position']
    # l_length = longest_orf_at_pos['longest_ORF_length']
    # print("Longest ORF in file with reading frame start position", start_pos, "is at identifier", identifier, "and has the length of", l_length)

    # Longest ORF in a given identifier i.e. gi|142022655|gb|EQ086233.1|16
    given_id = 'gi|142022655|gb|EQ086233.1|16'
    pos, orf = longest_ORF_in_given_seq(sys.argv[1], given_id)
    l_length = len(orf)
    print('Longest ORF for ', given_id, ' has length', l_length)
