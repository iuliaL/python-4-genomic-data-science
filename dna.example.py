
import read_fasta
import sys

"""(1) How many records are in the file? A record in a FASTA file is defined as a single-line header,
followed by lines of sequence data. The header line is distinguished from the sequence data by a greater-than (">") symbol in the first column.
The word following the ">" symbol is the identifier of the sequence, and the rest of the line is an optional description of the entry.
There should be no space between the ">" and the first letter of the identifier. """


def count_fasta_records(f):
    records = read_fasta.read(f)
    return len(records)
    
if __name__ == "__main__":
    sys.stdout.write("Records count: %s\n"  %count_fasta_records(sys.argv[1]))
