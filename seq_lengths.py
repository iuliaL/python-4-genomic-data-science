import read_fasta
import sys


def compute_lengths(recs):
    """ What are the lengths of the sequences in the fasta file?
        Returns a <id, sequence_length> dict"""
    records_length = {identifier: len(seq) for (identifier, seq) in recs.items()}
    return records_length


def shortest_seq(records_lengths):
    """ What are the identifiers of the shortest sequences?
        Returns a list"""
    min_value = min(records_lengths.values())
    return [k for k in records_lengths if records_lengths[k] == min_value]


def longest_seq(records_lengths):
    """ What are the identifiers of the longest sequences?
        Returns a list"""
    max_value = max(records_lengths.values())
    return [k for k in records_lengths if records_lengths[k] == max_value]


if __name__ == "__main__":
    f = sys.argv[1]
    parsed_records = read_fasta.read(f)
    # Longest sequences
    sys.stdout.write("Longest sequences are: " + " ".join(longest_seq(parsed_records)) + '\n')

    # Shortest sequences
    sys.stdout.write("Shortest sequences are: " + " ".join(shortest_seq(parsed_records)) + '\n')
