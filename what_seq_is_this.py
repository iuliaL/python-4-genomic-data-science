from Bio.Blast import NCBIWWW, NCBIXML

fasta_str = open("unknown_seq.fa").read()
result_handle = NCBIWWW.qblast("blastn", "nt", fasta_str)  # query blast
# program: 'blastn' searches nucleotides against nucleotides
# database: 'nt'
blast_record = NCBIXML.read(result_handle)

E_value_threshold = 0.01
for alignment in blast_record.alignments:
    for hsp in alignment.hsps:  # hsp -> High scoring pairs
        if hsp.expect < E_value_threshold:
            print('*** Alignment ***')
            print('sequence:', alignment.title)
            print('length:', alignment.length)
            print('e value:', hsp.expect)
            print(hsp.query)
            print(hsp.match)
            print(hsp.sbjct)
