#! usr/local/bin/python
import sys

# This is written for generating site-directed mutagenic primers (forward strand)
# and mutant proteins for gene synthesis purpose
# from a list of diagnostic residues (qinwen@uchicago.edu)

AncZenG = "accaaacgctctcgcaccgcttttacctcgctgcaactgatcgaactggaacgtgaatttcatattaacaagtatctgtgtcgtccgcgtcgcattgaaatcagccagcgtctgtctctgtccgaacgccaagtgaaaatctggttccaaaatcgccgtatgaaatcaaagaaggactcgtga"
AncBcdG = "ccgcgtcgcacccgcaccacctttacgtccgcccaaatcgctgaactggaacaacacttcctgcaaggccgctacctgaccgcatcacgtctggcagaactgagcgcaaaactggcactgggcaccgcacaggtgaaaatttggtttaaaaaccgtcgccgtcgccataaaatccaatcttga"

AncZenP = "TKRSRTAFTSLQLIELEREFHINKYLCRPRRIEISQRLSLSERQVKIWFQNRRMKSKKDS"
AncBcdP = "PRRTRTTFTSAQIAELEQHFLQGRYLTASRLAELSAKLALGTAQVKIWFKNRRRRHKIQS"

# a list containing locations of the substitutions.
# in this particular case, I want to mutate residues at following locations in AncBcd
# to the corresponding ones in AncZenBcd

## Salt-bridge and backbone contacts restoration mutations (AncBcdHD_re7) or
## Ancestral salt-bridge and backbone contacts disruption mutations (AncZenHD_f7)
##Mutation_loc = [19, 23, 28, 31, 42, 43, 58]

## 6 diagnostic residues on the recognition helix (AncBcdHD_bH3 and AncZenHD_fH3)
##Mutation_loc = [42, 43, 50, 54, 55, 58]

## All 13 diagnostic residues
##Mutation_loc = [2, 7, 13, 19, 23, 28, 31, 42, 43, 50, 54, 55, 58]

## Mutation combinations
Mutation_loc = [54]

def CreatePrimer(FromSeq, ToSeq, Mutations):
    # sequences at 5' and 3' outside the gene of interest (in LIC-MBP vector)
    Prime5 = "TACTTCCAATCCAATGCG"
    Prime3 = "CGCATTGGAAGTGGATAA"
    output = open('mutant_primers.txt', 'a')
    for i in Mutations:
        primerSeq = ''
        if i == 1:
            primerSeq = Prime5 + ToSeq[0:3] + FromSeq[3:24]
        if i > 1 and i < 8:
            primerSeq = Prime5[(3*i-6):18] + FromSeq[0:(3*(i-1))] + ToSeq[(3*i-3):(3*i)] + FromSeq[(3*i):(3*i+21)]
        if i > 7 and i < 54:
            primerSeq = FromSeq[3*(i-8):(3*(i-1))] + ToSeq[(3*i-3):(3*i)] + FromSeq[(3*i):(3*i+21)]
        if i > 53 and i < 61:
            primerSeq = FromSeq[3*(i-8):(3*(i-1))] + ToSeq[(3*i-3):(3*i)] + FromSeq[(3*i):183] + Prime3[0:3*(i-54)]
        output.write("backMut_" + str(i) + ";")
        output.write(primerSeq)
        output.write("\n")

CreatePrimer(AncZenG, AncBcdG, Mutation_loc)

# From a list of mutant positions, mutate the residues at the background "FromSeq"
# to the residues of "ToSeq"
def CreateProtein(FromSeq, ToSeq, backbone, Mutations):
    output = open('mutant_proteins.txt', 'a')
    for i in Mutations:
        FromSeq = FromSeq[0:(i-1)] + ToSeq[i-1] + FromSeq[i:60]
    output.write(">" + backbone + "\n")
    output.write(FromSeq)
    output.write("\n")

##CreateProtein(AncBcdP, AncZenP, "AncBcdHD_bH3", Mutation_loc)
##CreateProtein(AncZenP, AncBcdP, "AncZenHD_fH3", Mutation_loc)

## Create tandem point mutation primer
def CreatePrimer(FromSeq, ToSeq, Mutations):
    # sequences at 5' and 3' outside the gene of interest (in LIC-MBP vector)
    Prime5 = "TACTTCCAATCCAATGCG"
    Prime3 = "CGCATTGGAAGTGGATAA"
    output = open('mutant_primers.txt', 'a')
    for i in Mutations:
        primerSeq = ''
        if i == 1:
            primerSeq = Prime5 + ToSeq[0:3] + FromSeq[3:24]
        if i > 1 and i < 8:
            primerSeq = Prime5[(3*i-6):18] + FromSeq[0:(3*(i-1))] + ToSeq[(3*i-3):(3*i)] + FromSeq[(3*i):(3*i+21)]
        if i > 7 and i < 54:
            primerSeq = FromSeq[3*(i-8):(3*(i-1))] + ToSeq[(3*i-3):(3*i)] + FromSeq[(3*i):(3*i+21)]
        if i > 53 and i < 61:
            primerSeq = FromSeq[3*(i-8):(3*(i-1))] + ToSeq[(3*i-3):(3*i)] + FromSeq[(3*i):183] + Prime3[0:3*(i-54)]
        FromSeq = FromSeq[0:3*(i-1)] + ToSeq[3*(i-1):(3*i)] + FromSeq[(3*i):183]
        output.write("ForwardMut_" + str(i) + ";")
        output.write(primerSeq)
        output.write("\n")

##CreatePrimer(AncZenG, AncBcdG, Mutation_loc)










