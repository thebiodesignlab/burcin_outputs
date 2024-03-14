#!/bin/bash

ml MMseqs2/13-45111-gompi-2020b

#mmseqs createdb 1CEZ_A_formatted.fa 1CEZDB

#mmseqs cluster 1CEZDB clusters tmp --min-seq-id 0.7

#mmseqs createtsv 1CEZDB clusters clusters.tsv

mmseqs createseqfiledb 1CEZDB clusters clu_seq
mmseqs result2flat 1CEZDB 1CEZDB clu_seq clu_seq.fasta

#mmseqs createsubdb 1CEZDB DB DB_clu_rep
#mmseqs convert2fasta DB_clu_rep DB_clu_rep.fasta
