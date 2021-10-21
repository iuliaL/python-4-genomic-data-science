from Bio.Blast import NCBIWWW, NCBIXML

fasta_str = open("lowest_e_val.fa").read()
result_handle = NCBIWWW.qblast("blastn", "nt", fasta_str)  # query blast
# program: 'blastn' searches nucleotides against nucleotides
# database: 'nt'
blast_record = NCBIXML.read(result_handle)

lowest_e_value = 9999
lowest_e_val_alignment = None
for alignment in blast_record.alignments:
    for hsp in alignment.hsps:  # hsp -> High scoring pairs
        if hsp.expect < lowest_e_value:
            lowest_e_value = hsp.expect
            lowest_e_val_alignment = alignment

for hsp in lowest_e_val_alignment.hsps:  # hsp -> High scoring pairs
    print('*** Lowest e alignment is: ***')
    print('sequence:', lowest_e_val_alignment.title)
    print('length:', lowest_e_val_alignment.length)
    print('e value:', hsp.expect)
    print(hsp.query)
    print(hsp.match)
    print(hsp.sbjct)

"""
Output:
sequence: gi|1783584753|gb|MN651324.1| Nicotiana tabacum strain zhongyan90 cytoplasmic male sterility(CMS) line cultivar MSzhongyan90 mitochondrion, complete genome
length: 530869
e value: 9.96777e-96
"""
