import read_fasta
import sys


def compute_lengths(recs):
    """ What are the lengths of the sequences in the fasta file?
        Returns a <id, sequence_length> dict"""
    records_length = {identifier: len(seq) for (identifier, seq) in recs.items()}
    return records_length


def shortest_seq(records):
    """ What are the identifiers of the shortest sequences?
        Returns a tuple with a list of all shortest sequences and the length"""
    lengths_dict = compute_lengths(records)
    min_value = min(lengths_dict.values())
    return [k for k in lengths_dict if lengths_dict[k] == min_value], min_value


def longest_seq(records):
    """ What are the identifiers of the longest sequences?
        Returns a tuple with a list of all longest sequences and the length"""
    lengths_dict = compute_lengths(records)
    max_value = max(lengths_dict.values())
    return [k for k in lengths_dict if lengths_dict[k] == max_value], max_value


if __name__ == "__main__":
    f = sys.argv[1]
    parsed_records = read_fasta.read(f)
    # Longest sequences
    long_sequences, l_length = longest_seq(parsed_records)
    sys.stdout.write(
        "Longest sequences and their length : " +
        " ".join(
            long_sequences) +
        '\n -> ' +
        str(l_length) + '\n')

    # Shortest sequences
    short_sequences, s_length = shortest_seq(parsed_records)
    sys.stdout.write("Shortest sequences and their length: " + " ".join(
        short_sequences) +
        '\n -> ' +
        str(s_length) + '\n')
