
import sys
import getopt


def usage():
    print(
        """Returns a dictionary of all sequences read from a FASTA file where the keys are the sequences names.
Suppose I want an extra condition and I want to optionally store the sequences that are more than 250 bases long.

i.e. read_fasta.py -l 250 fasta_example.py

Usage:
    read_fasta.py [-h] [-l <length>] <filename>

    -h                  means print this message
    -l <length>         filter out sequences shorter than <length> (default = 0)
    <filename>          FASTA file name
""")


def read(f):
    try:
        file = open(f)
    except IOError:
        sys.exit("Error opening file")
    seq_dictionary = {}
    for line in file:
        line = line.rstrip()
        if line[0].startswith(">"):
            # means is a info line not nucleotide line
            seq_name = line.split()[0][1:]
            seq_dictionary[seq_name] = ''
        else:  # means i have nucleotides
            seq_dictionary[seq_name] += line
    file.close()
    return seq_dictionary


if __name__ == "__main__":
    optional_args, required_args = getopt.getopt(sys.argv[1:], 'hl:')
    # sys.argv starts at 1 because 0 is the very script that is intended to run in this case read_fasta.py
    # ':' after l means -l needs a parameter after typing the flag
    options = {}
    min_seq_len = 0

    for k, v in optional_args:
        options[k] = v  # this will look like { '-h': '', '-l': 250 }

    if '-h' in options.keys():
        usage()
        sys.exit()

    if len(required_args) < 1:
        usage()
        sys.stderr.write("No input FASTA file provided\n")
        sys.exit()

    if "-l" in options.keys():
        if int(options['-l']) < 0:
            sys.exit("Positive integer length required")
        else:
            min_seq_len = int(options['-l'])

    seq_dict = read(required_args[0])
    filtered_dict = {name: seq for (name, seq) in seq_dict.items() if len(seq) >= min_seq_len}  # dict comprehension

    print("The sequences with at least", min_seq_len, "bases are: \n", filtered_dict)
