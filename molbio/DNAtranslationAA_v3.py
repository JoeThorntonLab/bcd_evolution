#! usr/bin/env python
from sys import argv
script, filename = argv

#This is used to translate any give DNA sequences to protein sequences according
#to the standard genetic codons. It generates translated protein sequences from all
#three frames and two directions.
#Qinwen Liu (qinwen@uchicago.edu)

code = {'ttt': 'F', 'tct': 'S', 'tat': 'Y', 'tgt': 'C',
       'ttc': 'F', 'tcc': 'S', 'tac': 'Y', 'tgc': 'C',
       'tta': 'L', 'tca': 'S', 'taa': '*', 'tga': '*',
       'ttg': 'L', 'tcg': 'S', 'tag': '*', 'tgg': 'W',
       'ctt': 'L', 'cct': 'P', 'cat': 'H', 'cgt': 'R',
       'ctc': 'L', 'ccc': 'P', 'cac': 'H', 'cgc': 'R',
       'cta': 'L', 'cca': 'P', 'caa': 'Q', 'cga': 'R',
       'ctg': 'L', 'ccg': 'P', 'cag': 'Q', 'cgg': 'R',
       'att': 'I', 'act': 'T', 'aat': 'N', 'agt': 'S',
       'atc': 'I', 'acc': 'T', 'aac': 'N', 'agc': 'S',
       'ata': 'I', 'aca': 'T', 'aaa': 'K', 'aga': 'R',
       'atg': 'M', 'acg': 'T', 'aag': 'K', 'agg': 'R',
       'gtt': 'V', 'gct': 'A', 'gat': 'D', 'ggt': 'G',
       'gtc': 'V', 'gcc': 'A', 'gac': 'D', 'ggc': 'G',
       'gta': 'V', 'gca': 'A', 'gaa': 'E', 'gga': 'G',
       'gtg': 'V', 'gcg': 'A', 'gag': 'E', 'ggg': 'G'}

DNA_match = {'a': 't', 't': 'a', 'c': 'g', 'g': 'c'}

def FormatFile(filename):
    seqfile = open(filename, 'r')
    seq = ''
    for line in seqfile.readlines():
        if line[0] != '>':
            line = line.lower()
            seq = seq + line.replace('\n', '')
    seqfile.close()
    return seq

seq = FormatFile(filename)

def ReverseComplement(seq, DNA_match):
    revSeq = ''
    n = len(seq) - 1
    while n >= 0:
        locus = seq[n]
        if locus in 'atgc':
            revSeq = revSeq + DNA_match[locus]
        else:
            revSeq = revSeq + '-'
        n = n - 1
    return revSeq

def translate(seqDNA, code):
    prot = ''
    for i in range(0, len(seqDNA), 3):
        if (i+3) > len(seqDNA):
    	    break
        else:
            codon = seqDNA[i:i+3]
            if codon in code.keys():
                prot = prot + code[codon]
            else:
                prot = prot + '-'
    return prot

output = open('translatedAA_Pco-zen.fasta', 'w')

for i in range(3):
    name = '>frame_' + str(i+1)
    output.write(name + '\n')
    DNA = seq[i:]
    trans = translate(DNA, code)
    output.write(trans + '\n')    

for i in range(3):
    name = '>reverseframe_' + str(i+1)
    output.write(name + '\n')
    rev_seq = ReverseComplement(seq, DNA_match)
    DNA = rev_seq[i:]
    trans = translate(DNA, code)
    output.write(trans + '\n')    
   
output.close()

    
