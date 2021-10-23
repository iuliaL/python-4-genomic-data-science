import read_fasta
import sys
import seq_lengths


"""A repeat is a substring of a DNA sequence that occurs in multiple copies (more than one) somewhere in the sequence.
   Although repeats can occur on both the forward and reverse strands of the DNA sequence, we will only consider repeats
   on the forward strand here. Also we will allow repeats to overlap themselves.
   For example, the sequence ACACA contains two copies of the sequence ACA - once at position 1 (index 0 in Python), and once at position 3.
   Given a length N, the program should be able to identify all repeats of length N in all sequences in the FASTA file.
   The program should also determine how many times each repeat occurs in the file, and which is the most frequent repeat of a given length.
   """


def all_repeats_in_all_seqs(N, f):
    """Identify all repeats of length N in all sequences in the FASTA file
    Returns dict of repeats: times"""
    # read fasta file
    seqs_dict = read_fasta.read(f)
    all_seqs = seqs_dict.values()
    # initialize empty dict for storage
    candidate_repeats = {}
    # go through all the positions
    for seq in all_seqs:
        for i in range(len(seq) - N + 1):
            current_Nmer = seq[i:i + N]
            # store all N-mers as keys of dict and store the number of occurence in values of dict
            candidate_repeats[current_Nmer] = candidate_repeats.get(current_Nmer, 0) + 1
    # filter out the ones with value < 2, cause they occur once only so they aren't repeats
    filtered = {repeat: times for repeat, times in candidate_repeats.items() if times >= 2}
    return filtered


def most_freq_repeat(N, f):
    repeats = all_repeats_in_all_seqs(N, f)
    max_freq = max(repeats.values())
    return [k for k in repeats if repeats[k] == max_freq], max_freq


N = 12
# print('Repeats dict', all_repeats_in_all_seqs(N, sys.argv[1]))
most_freq_repeats, times = most_freq_repeat(N, sys.argv[1])
print('Most frequent repeats are:', most_freq_repeats, " and they occur", times, "times")
