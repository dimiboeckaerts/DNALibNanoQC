"""
NANOPORE ANALYSIS

Created on 19/10/2020

@author: dimiboeckaerts

TO DO & QUESTIONS:
- how much faster is Julia compared to Python?
- how to go from Cedric's version to something we can use?
- check the positions: count number of '-' in alignment? Test a single
    read with tiles in and see if they are detected and how (check alignment).
"""

# 0 - LIBRARIES
# --------------------------------------------------
#Pkg.add("DataFrames")
using BioAlignments
using BioSequences
using DataFrames
using FASTX
using ProgressMeter


# 1 - FUNCTIONS
# --------------------------------------------------
# convert file to a list
function fasta_to_array(file)
    """
    Function that reads a FASTA file and puts its sequences in an array.
    """
    sequences = []
    ids = []
    reader = FASTA.Reader(open(file, "r"))
    for record in reader
        seq = FASTA.sequence(record)
        identifier = FASTA.identifier(record)
        push!(sequences, seq)
        push!(ids, identifier)
    end
    return sequences, ids
end

function fastq_to_array(file)
    """
    Function that reads a FASTQ file and puts its sequences in an array.
    """
    sequences = []
    reader = FASTQ.Reader(open(file, "r"))
    for record in reader
        seq = FASTQ.sequence(record)
        push!(sequences, seq)
    end
    return sequences
end

# calculate alingment and its score
function local_alignment_score(sequence1, sequence2)
    """
    This function calculates the score of matches between two aligned DNA/RNA
    sequences.
    """
    scoremodel = AffineGapScoreModel(EDNAFULL, gap_open=-10, gap_extend=-1)
    res = pairalign(LocalAlignment(), sequence1, sequence2, scoremodel);
    aln_score = score(res)

    return res, aln_score
end


# 2 - ANALYSIS
# --------------------------------------------------
# position markers
PM1 = "ACCATG"
PM2 = "GGTGCA"
PM3 = "GCAGGT"
PM4 = "GGAAGC"
PM5 = "AAGTAA"

# Load all tiles and FASTQ into arrays
tiles_p1, ids_p1 = fasta_to_array("/Users/dimi/Documents/GitHub_Local/DNALibNanoQC/sequences/tiles_position1.fasta")
tiles_p2, ids_p2 = fasta_to_array("/Users/dimi/Documents/GitHub_Local/DNALibNanoQC/sequences/tiles_position2.fasta")
tiles_p3, ids_p3 = fasta_to_array("/Users/dimi/Documents/GitHub_Local/DNALibNanoQC/sequences/tiles_position3.fasta")
tiles_p4, ids_p4 = fasta_to_array("/Users/dimi/Documents/GitHub_Local/DNALibNanoQC/sequences/tiles_position4.fasta")
#tiles = fasta_to_array("tiles.fasta") # more general

nanopore_reads = fastq_to_array("/Users/dimi/Documents/GitHub_Local/DNALibNanoQC/test/nanopore_reads.fastq")
nanopore_mini = nanopore_reads[12000:14000]
nanopore_toy, toy_ids = fasta_to_array("/Users/dimi/Documents/GitHub_Local/DNALibNanoQC/test/toy_reads.fasta")

# add position markers
tiles_p1 = [PM1*string(tile)*PM2 for tile in tiles_p1]
tiles_p2 = [PM2*string(tile)*PM3 for tile in tiles_p2]
tiles_p3 = [PM3*string(tile)*PM4 for tile in tiles_p3]
tiles_p4 = [PM4*string(tile)*PM5 for tile in tiles_p4]

# do local alignments per position
data = nanopore_reads
id_list_p1 = []
alignments_list_p1 = []; alignments_start_p1 = []; alignments_stop_p1 = []
id_list_p2 = []
alignments_list_p2 = []; alignments_start_p2 = []; alignments_stop_p2 = []
id_list_p3 = []
alignments_list_p3 = []; alignments_start_p3 = []; alignments_stop_p3 = []
id_list_p4 = []
alignments_list_p4 = []; alignments_start_p4 = []; alignments_stop_p4 = []
bar = Progress(length(data))

for nano_read in data
    # position 1
    best_score_p1 = 0
    best_tile_p1 = ""; best_id_p1 = ""
    best_alignment_p1 = "" # the alignment of the tile
    start_p1 = -1; stop_p1 = -1
    for (i, tile) in enumerate(tiles_p1)
        tile_res, tile_score = local_alignment_score(LongDNASeq(tile), nano_read)
        if tile_score > best_score_p1
            best_score_p1 = tile_score
            best_tile_p1 = tile
            best_id_p1 = ids_p1[i]
            best_alignment_p1 = tile_res
            start_p1 = tile_res.aln.a.aln.anchors[1].refpos
            stop_p1 = tile_res.aln.a.aln.anchors[2].refpos
        end
    end
    push!(alignments_list_p1, best_alignment_p1)
    push!(id_list_p1, best_id_p1)
    push!(alignments_start_p1, start_p1)
    push!(alignments_stop_p1, stop_p1)

    # position 2
    best_score_p2 = 0
    best_tile_p2 = ""; best_id_p2 = ""
    best_alignment_p2 = "" # the alignment of the tile
    start_p2 = -1; stop_p2 = -1
    for (i, tile) in enumerate(tiles_p2)
        tile_res, tile_score = local_alignment_score(LongDNASeq(tile), nano_read)
        if tile_score > best_score_p2
            best_score_p2 = tile_score
            best_tile_p2 = tile
            best_id_p2 = ids_p2[i]
            best_alignment_p2 = tile_res
            start_p2 = tile_res.aln.a.aln.anchors[1].refpos
            stop_p2 = tile_res.aln.a.aln.anchors[2].refpos
        end
    end
    push!(alignments_list_p2, best_alignment_p2)
    push!(id_list_p2, best_id_p2)
    push!(alignments_start_p2, start_p2)
    push!(alignments_stop_p2, stop_p2)

    # position 3
    best_score_p3 = 0
    best_tile_p3 = ""; best_id_p3 = ""
    best_alignment_p3 = "" # the alignment of the tile
    start_p3 = -1; stop_p3 = -1
    for (i, tile) in enumerate(tiles_p3)
        tile_res, tile_score = local_alignment_score(LongDNASeq(tile), nano_read)
        if tile_score > best_score_p3
            best_score_p3 = tile_score
            best_tile_p3 = tile
            best_id_p3 = ids_p3[i]
            best_alignment_p3 = tile_res
            start_p3 = tile_res.aln.a.aln.anchors[1].refpos
            stop_p3 = tile_res.aln.a.aln.anchors[2].refpos
        end
    end
    push!(alignments_list_p3, best_alignment_p3)
    push!(id_list_p3, best_id_p3)
    push!(alignments_start_p3, start_p3)
    push!(alignments_stop_p3, stop_p3)

    # position 4
    best_score_p4 = 0
    best_tile_p4 = ""; best_id_p4 = ""
    best_alignment_p4 = "" # the alignment of the tile
    start_p4 = -1; stop_p4 = -1
    for (i, tile) in enumerate(tiles_p4)
        tile_res, tile_score = local_alignment_score(LongDNASeq(tile), nano_read)
        if tile_score > best_score_p4
            best_score_p4 = tile_score
            best_tile_p4 = tile
            best_id_p4 = ids_p4[i]
            best_alignment_p4 = tile_res
            start_p4 = tile_res.aln.a.aln.anchors[1].refpos
            stop_p4 = tile_res.aln.a.aln.anchors[2].refpos
        end
    end
    push!(alignments_list_p4, best_alignment_p4)
    push!(id_list_p4, best_id_p4)
    push!(alignments_start_p4, start_p4)
    push!(alignments_stop_p4, stop_p4)

    # update bar
    next!(bar)
end

# Build DataFrame
results = DataFrame(
    id1=id_list_p1, id2=id_list_p2, id3=id_list_p3, id4=id_list_p4,
    aln1=alignments_list_p1, start1=alignments_start_p1, stop1=alignments_stop_p1,
    aln2=alignments_list_p2, start2=alignments_start_p2, stop2=alignments_stop_p2,
    aln3=alignments_list_p3, start3=alignments_start_p3, stop3=alignments_stop_p3,
    aln4=alignments_list_p4, start4=alignments_start_p4, stop4=alignments_stop_p4)


# 3 - TEST
# --------------------------------------------------
s1 = dna"CCTAGGAGGG";
s2 = dna"ACCTGGTATGATAGCG";

#get shortest sequence of an alignment:
first_seq = res.aln.a
first_seq = split(string(first_seq), "\n")[2]



nano_r = LongDNASeq("ATCAATTGCTTCGTTCGATTTAAATATTGCTAAGGTTAAGGATTCATTCCCACGGTAACACCAGCACCTGTTAGCAGAATGAATCACCGATGCGAGCGAACGTGAAGCGACTGCTGCTCTAAAACGTCTGCGACCTGAGCAACAACATAGATGGTCTTCGATTTCCGTGTTTCGTAAGAGAAGTCTGGAAGCATGAAGTCCCCACGTGCTGGTTCTTATTCTGAGTTACAACAGTCACACCGCTGCCGGTAGCATTTCCCTTCCGGTGGGCGCGGGGGCATGACCTCGTCGCCGCACTTATGACTGTCTTCTTTATCATGCAACTCGTAGGACAGGAGATGTTCCGGATGCCATCGCGGGATGCTGCTGGCTACCTGTGGAGAACACCTACATCTCTGTGATAACGAAGCGCTAACCACCGTTTTATCAGGCTCTGGGGAGGCAGAATAAATGATCATATCGTCAATTATTACCTCCACGAGGAGGCCTGAGCAAACTGGCCTCAGGCATTTGAGAAGCACACGGTCACACTGCTTCGGTAGTCAATAAACCGGTAAACCAGCAATGAACATAAGCGGCTGTTAACGACCCTGCTGAACCGACGACCGGGTCGGGCCCTTTCAATTTCTGCCATTCATCCGCTTATTATCACTTATTCAGGCGTAGCAACCAGCGGGCGTTTAAGGGCACCAATAACTGCCTTAAAATTACGCCCGCCCTGCCACTCATCACTGGTACTGTTGTAATTCATTAAGCATTCTGCTT")
tile = LongDNASeq("ACCATGATTCTGGGCAAAATTTGGAAAGGCATTAAAAGCCTGTTTGGTGCA")

scoremodel = AffineGapScoreModel(EDNAFULL, gap_open=-10, gap_extend=-1)
res = pairalign(LocalAlignment(), tile, nano_r, scoremodel);

local_alignment_score(tile, nano_r)
