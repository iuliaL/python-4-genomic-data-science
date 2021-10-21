
#!/usr/bin/python

import os
from time import process_time
import random

dna = "atgcaaagtaccggt"
c = dna.count("c")
g = dna.count("g")
sumcg = c + g
length = len(dna)
result = sumcg * 100 / length


"""Nicer outputting: what i want to print is 45.455%. So 3f means round with 3 decimal points.
5 means the number of digits but if there are more needed then they will print either way.
If less are needed the result is padded with spaces in front. Double %% means I want to escape the literal % I need in my output.
The % result part means feeding the argument"""
to_print = "The GC percentage is %5.3f%%" % result + ".\n"
# print(to_print)


lst = ["a", "b", "c", "d", "e"]
lst_copy = lst[:]
# print("lst == lst_copy", lst == lst_copy, "\nlst 'is' lst_copy", lst is lst_copy)  # !! Identity comparison -> is / is not
# so == returns True for copies of objects, and is returns True only if the objects point to the same object in memory
# print(lst[2:4]) # -> ["c", "d"]
# lst[2:4] = "x"  # mutates the list (replaces c, d with x)
# lst[:] = []  # clears a list also lst = [] "clears" a list
# del lst[2:4]
# i s the same as
# lst[2:4] = []


# print("List printed %s" % lst)

thisdict = {
    "year": 1964,
    "brand": "Ford",
    "model": "Mustang"
}
# print(sorted(thisdict.keys()))
# print(list(thisdict.values()))


# print("Dna length", len(dna))
""" Donor splice means gt seq. """
# pos = dna.find("gt", 0) # or simply dna.find("gt")
# while pos > -1: # meaning if i find at least one occurence
#     print("Donor splice site candidate at position %d" % pos)
#     pos = dna.find("gt", pos + 1)
""" this loop will eval to false and it will stop when dna.find will output -1"""
# the very same problem done with for loop
# for pos,v in enumerate(dna[:-1]):
#     if (dna[pos] + dna[pos + 1]) == "gt":
#         print("Donor splice site candidate at position %d" % pos)
# protein = "ADSDVHIRKUUUGW"
# corrected_protein = ""
# for i in range(len(protein)):
#     if protein[i] not in "ABCD":
#         # print("Invalid aminoacid at pos %d: '%s'" % (i, protein[i]))
#         continue # do nothing
#     corrected_protein += protein[i]
# print("Corrected protein", corrected_protein)

""" Find all prime numbers smaller than a given integer"""
""" a prime number means that it must divide to itself only (and to 1)"""
N = 14
primes = []
for y in range(2, N):
    for x in range(2, y):
        if y % x == 0:
            break
    else:  # !!! because i used a break
        primes.append(y)

# uglier, more verbose and less efficient
# for y in range(2, N):
#     isPrime = True
#     for x in range(2, y):
#         if y%x == 0:
#             isPrime = False
#     if isPrime:
#         primes.append(y)
# print(primes)

# seq = dna
# print("Sequence is", seq)
# for i in range(len(seq) + 1):     # line 1
#     for j in range(i):        # line 2
#         # print("j", j, "i", i)
#         print(seq[j:i])     # line 3

# i=0
# while i<len(seq) :
#       j=0
#       while(j<i+1) :
#                 print(seq[j:i])
#                 j=j+1
#       i=i+1


def gc(dna):
    "this function should the percentage of the occurences of 'G' and 'C in dna string"
    # remove possible undefined bases
    number_of_undefined_bases = dna.count("N") + dna.count("n")
    bases_to_check = len(dna) - number_of_undefined_bases
    g = dna.count('g')
    c = dna.count('c')
    return (g + c) * 100 / bases_to_check
# print(gc(dna))
# print(help(gc))


def has_STOP_codon(dna, frame=0):
    "Checks if a dna seq has an in-frame STOP codon. Frame means position where to start"
    stop_codons = ["TGA", "TAG", "TAA"]
    for i in range(frame, len(dna), 3):  # step 3 because I don't want to check each codon only stepping
        curr_codon = dna[i:i + 3].upper()
        if curr_codon in stop_codons:
            return True
    return False


# print(has_STOP_codon("atgagcggccggct", 1)) => True
# print(has_STOP_codon("atgagcggccggct")) => False
def reverse_string(dna):
    return dna[::-1]


dna = "AIGTGTGGGGCG"


def reverse_complement(dna):
    "Return the complementary string sequence"
    complements = dict(a="t", t="a", g="c", c="g", n="n")
    list_of_complementary_bases = [complements[base.lower()] for base in reverse_string(dna) if base.lower() in complements]
    return "".join(list_of_complementary_bases)


# print(reverse_complement(dna))


def compute(n, x, y):
    if n == 0:
        return x
    return compute(n - 1, x + y, y)


# Fastest counts

def create_long_dna(n, alphabet="ACGT"):
    return "".join([random.choice(alphabet) for i in range(n)])


def count1(dna, base):
    t0 = process_time()
    i = 0
    for c in dna:
        if c == base:
            i += 1
    elapsed = (process_time() - t0)
    return i, elapsed, "count1"


def count2(dna, base):
    t0 = process_time()
    i = 0
    for j in range(len(dna)):
        if dna[j] == base:
            i += 1
    elapsed = (process_time() - t0)
    return i, elapsed, "count2"


def count3(dna, base):
    t0 = process_time()
    match = [c == base for c in dna]
    elapsed = (process_time() - t0)
    return sum(match), elapsed, "count3"


def count4(dna, base):
    t0 = process_time()
    count = dna.count(base)
    elapsed = (process_time() - t0)
    return count, elapsed, "count4"


def count5(dna, base):
    t0 = process_time()
    count = len([i for i in range(len(dna)) if dna[i] == base])
    elapsed = (process_time() - t0)
    return count, elapsed, "count5"


def count6(dna, base):
    t0 = process_time()
    count = sum(c == base for c in dna)
    elapsed = (process_time() - t0)
    return count, elapsed, "count6"


Dna = create_long_dna(1000000)

data = [count1(Dna, "T"),
        count3(Dna, "T"),
        count2(Dna, "T"),
        count4(Dna, "T"),
        count5(Dna, "T"),
        count6(Dna, "T")]

fastest_function = min(data, key=lambda t: t[1])
whichFn, seconds = fastest_function[2], fastest_function[1]
# print("Min is {} with {} seconds".format(whichFn, seconds))
# count 4 wins. Basically the python implementation of count wins


def read_FASTA(file_name):
    """Returns a dictionary of the sequences read where the keys are the seq. names"""
    try:
        file = open(file_name)
    except IOError:
        print('No file found', file_name)
        return None
    seq_dictionary = {}
    for line in file:
        line = line.rstrip()
        if line.startswith(">"):
            # means is a info line not nucleotide line
            seq_name = line[1:]
            seq_dictionary[seq_name] = ''
        else:  # means i have nucleotides
            seq_dictionary[seq_name] += line
    file.close()
    return seq_dictionary

# print(read_FASTA(os.path.expanduser("~/Desktop/fasta_example.fa")))


# print to file
file = open(os.path.expanduser("~/Desktop/my-test.txt"), "w+")
file.write("%s" % to_print)
file.close()
