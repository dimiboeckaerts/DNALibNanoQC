import os
import sys
import pandas as pd
import numpy as np

from Bio import pairwise2
from Bio import SeqIO

from Bio.pairwise2 import format_alignment

PROJDIR = os.path.dirname(os.path.realpath(__file__))
FASTA   = os.path.join(PROJDIR, "sequences")
FASTQ   = os.path.join(PROJDIR, "reads")
READS   = sys.argv[1]

OUTPUT  = os.path.join(sys.argv[2])
f = open(OUTPUT, "x")

def fasta_upper(pos):
    for fasta in pos:
        fasta.seq = fasta.seq.upper()
    return pos

def print_best_match(counter, alignment, tile, position):
    f.write(str(counter) + ",")                   # read number
    f.write(str(tile.name) + ",")                 # tile
    f.write(str(position) + ",")                  # tile position
    f.write(str(len(tile.seq)) + ",")             # lenght of tile
    f.write(str(alignment[4]-alignment[3]) + ",") # length alignment
    f.write(str(alignment[3]) + ",")              # start
    f.write(str(alignment[4]) + ",")              # end
    f.write(str(alignment[2]) + "\n")             # score
    f.flush()

    
# open each file of sequences and make sure they are upper case (pairwise2 case sensitive)
pos1 = fasta_upper(list(SeqIO.parse(os.path.join(FASTA, "tiles_position1.fasta"), "fasta")))
pos2 = fasta_upper(list(SeqIO.parse(os.path.join(FASTA, "tiles_position1.fasta"), "fasta")))
pos3 = fasta_upper(list(SeqIO.parse(os.path.join(FASTA, "tiles_position3.fasta"), "fasta")))
pos4 = fasta_upper(list(SeqIO.parse(os.path.join(FASTA, "tiles_position4.fasta"), "fasta")))

# the tiles have position specific, linker sequences 
for fasta in pos1:
    fasta.seq = "ACCATG" + fasta.seq + "GGTGCA"

for fasta in pos2:
    fasta.seq = "GGTGCA" + fasta.seq + "GCAGGT"

for fasta in pos3:
    fasta.seq = "GCAGGT" + fasta.seq + "GGAAGC"
                    
for fasta in pos4:
    fasta.seq = "GGAAGC" + fasta.seq + "AAGTAA"

results = []
counter = 1
read_count = sum(1 for line in open(os.path.join(FASTQ, READS), "r"))/4

# configuration of localms parameters
match = 5
mismatch = -4
gapopen = -8
gapext = -4

# open fastq file and iterate through each read to examine tile content
with open(os.path.join(FASTQ, READS), "r") as handle:
    f.write("ReadId,TileName,TilePosition,TileLength,AlignmentLength,Start,Stop,AlignmentScore,Position\n")
    
    for read in SeqIO.parse(handle, "fastq"):
        print("Read %s of %s" % (str(counter), str(read_count)))

        best_score_pos1 = 0
        best_align_pos1 = []
        best_tile1 = ""

        for fasta1 in pos1:
            aln = pairwise2.align.localms(fasta1.seq, read.seq, match, mismatch, gapopen, gapext)
            for alignment in aln:
                if (alignment and alignment[2] > best_score_pos1):
                    best_score_pos1 = alignment[2]
                    best_align_pos1 = alignment
                    best_tile1 = fasta1
                
        best_score_pos2 = 0.0
        best_align_pos2 = []
        best_tile2 = ""

        for fasta2 in pos2:
            aln = pairwise2.align.localms(fasta2.seq, read.seq, match, mismatch, gapopen, gapext)
            for alignment in aln:
                if (alignment and alignment[2] > best_score_pos2):
                    best_score_pos2 = alignment[2]
                    best_align_pos2 = alignment
                    best_tile2 = fasta2

        best_score_pos3 = 0.0
        best_align_pos3 = []
        best_tile3 = ""
        
        for fasta3 in pos3:
            aln = pairwise2.align.localms(fasta3.seq, read.seq, match, mismatch, gapopen, gapext)
            for alignment in aln:
                if (alignment and alignment[2] > best_score_pos3):
                    best_score_pos3 = alignment[2]
                    best_align_pos3 = alignment
                    best_tile3 = fasta3
         
        best_score_pos4 = 0.0
        best_align_pos4 = []
        best_tile4 = ""
        for fasta4 in pos4:
            aln = pairwise2.align.localms(fasta4.seq, read.seq, match, mismatch, gapopen, gapext)
            for alignment in aln:
                if (alignment and alignment[2] > best_score_pos4):
                    best_score_pos4 = alignment[2]
                    best_align_pos4 = alignment
                    best_tile4 = fasta4
                    
 
        # processing results
        print_best_match(counter, best_align_pos1, best_tile1, 1)
        print_best_match(counter, best_align_pos2, best_tile2, 2)
        print_best_match(counter, best_align_pos3, best_tile3, 3)
        print_best_match(counter, best_align_pos4, best_tile4, 4)       

        counter += 1

f.close()

# Adding a few calculated value and filtering the result file (csv)
blocks = pd.read_csv(OUTPUT, header = 0, sep=",")
blocks["Distance"] = 0
blocks["Suspicious"] = 0

print(blocks.head())

pd.set_option("display.max_rows", None)

# for each group of 4 blocks (belonging to a read), we calculate
# the distance between them, and flag suspicious reads.
for i in range(len(blocks)):
    if(blocks.iloc[i, 3] == 4):
        continue
    # distance between two consecutive tiles in a read
    blocks.iloc[i, 8] = blocks.iloc[i+1, 5] - blocks.iloc[i, 6]
    # if distance is greater than 50 nt, flagged as suspicious
    if(abs(blocks.iloc[i, 8]) > 50):
        blocks.iloc[i, 9] = 1


# we filter out to a file the reads that have a suspicious
# alignment
reads_pass = blocks.iloc[0:0,:].copy()
reads_fail = blocks.iloc[0:0,:].copy()
for i, g in blocks.groupby(np.arange(len(blocks)) // 4):
    flag = 0
    for j in range(0,4):
        if (g.iloc[j, 9] == 1):
            flag = 1
        break
    if(flag == 1):
        reads_fail = reads_fail.append(g)
    else:
        reads_pass = reads_pass.append(g)

reads_fail.to_csv('reads_fail.csv')
reads_pass.to_csv('reads_pass.csv')
