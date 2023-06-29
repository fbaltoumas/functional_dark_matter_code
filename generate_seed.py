#!/usr/bin/env python3
# generate_seed.py
# run as:
#  python3 generate_seed.py list_of_fastas
# where ‘list_of_fastas’ contains the names of the input MSAs in this format:
# F105000.fasta
# F045821.fasta
# F000388.fasta
from sys import argv,stderr
import multiprocessing as mp
from Bio import AlignIO
from Bio.Align import AlignInfo
from prody import *
#global variables to be changed by the user
SEQ_ID=0.9 #sequence identity cutoff, 90%
COV=0.75 # alignment coverage, 75%
CPUS=16 # number of CPUs to use for multiprocessing

def create_seed_and_consensus(aln_file):
  #parse the input MSA
  msa=parseMSA(aln_file, format="FASTA")
  name=aln_file.rstrip(".fasta")
  #produce seed alignment based on SEQ_ID and COV
  seed_msa=refineMSA(msa, seqid=SEQ_ID, rowocc=COV)
  #if for some reason, the applied COV cut-off removes ALL sequences,
  #recalculate the alignment using only the SEQID and warn the user
  if seed_msa.numSequences()==0:
    print("Applying alignment coverage %s resulted in empty alignment, using only the seq-id cutoff instead..." %COV, file=stderr)
    seed_msa=refineMSA(msa, seqid=SEQ_ID)
  writeMSA("%s_seed.fasta" %name, seed_msa)
  #parse the seed MSA as an AlignIO object
  alignment=AlignIO.read("%s_seed.fasta" %name,"fasta")
  #calculate the alignment summary
  summary=AlignInfo.SummaryInfo(alignment)
  #create a consensus sequence
  cons=summary.dumb_consensus()
  #write the consensus in a fasta file
  cons_out=open("%s_consensus.fasta" %name,"w")
  cons_out.write(">%s_consensus\n%s\n"%(name,cons))
  cons_out.close()
  print(">%s_consensus\n%s"%(name,cons), file=stderr)
#code to run the create_seed_and_consensus method in parallel
pool = mp.Pool(CPUS)
jobs = []
f=open(argv[1],"r")
files=f.readlines()
f.close()
dataset=[i.rstrip() for i in files]
run=pool.map(create_seed_and_consensus, dataset)