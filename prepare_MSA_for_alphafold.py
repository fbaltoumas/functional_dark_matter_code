import os
import subprocess as sp
import argparse
import string, sys, getopt

#numpy and sklearn
import numpy as np
from sklearn.metrics import pairwise_distances

# Sergei's custom script utils.py
#from utils import *
alpha_1 = list("ARNDCQEGHILKMFPSTWYV-")
states = len(alpha_1)
alpha_3 = ['ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE',
           'LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL','GAP']
aa_1_N = {a:n for n,a in enumerate(alpha_1)}
aa_3_N = {a:n for n,a in enumerate(alpha_3)}
aa_N_1 = {n:a for n,a in enumerate(alpha_1)}
aa_1_3 = {a:b for a,b in zip(alpha_1,alpha_3)}
aa_3_1 = {b:a for a,b in zip(alpha_1,alpha_3)}

def AA_to_N(x):
  # ["ARND"] -> [[0,1,2,3]]
  x = np.array(x);
  if x.ndim == 0: x = x[None]
  return [[aa_1_N.get(a, states-1) for a in y] for y in x]

def N_to_AA(x):
  # [[0,1,2,3]] -> ["ARND"]
  x = np.array(x);
  if x.ndim == 1: x = x[None]
  return ["".join([aa_N_1.get(a,"-") for a in y]) for y in x]

def parse_fasta(filename, a3m=False):
  '''function to parse fasta file'''
  if a3m:
    # for a3m files the lowercase letters are removed
    # as these do not align to the query sequence
    rm_lc = str.maketrans(dict.fromkeys(string.ascii_lowercase))
  header, sequence = [],[]
  lines = open(filename, "r")
  for line in lines:
    line = line.rstrip()
    if len(line) > 0:
      if line[0] == ">":
        header.append(line[1:])
        sequence.append([])
      else:
        if a3m: line = line.translate(rm_lc)
        else: line = line.upper()
        sequence[-1].append(line)
  lines.close()
  sequence = [''.join(seq) for seq in sequence]
  return header, sequence

def mk_msa(seqs):
  '''one hot encode msa'''
  alphabet = list("ARNDCQEGHILKMFPSTWYV-")
  states = len(alphabet)

  alpha = np.array(alphabet, dtype='|S1').view(np.uint8)
  msa = np.array([list(s) for s in seqs], dtype='|S1').view(np.uint8)
  for n in range(states):
    msa[msa == alpha[n]] = n
  msa[msa > states] = states-1

  return np.eye(states)[msa]





#define the location of the hhfilter exe from hhsuite
hhfilter="/usr/bin/hhfilter"







##### Command Line Interface (CLI) #####
parser = argparse.ArgumentParser(prog="prepare_MSA_for_alphafold.py", description="script to refine a multiple sequence alignment")
parser.version=0.1
help=parser._action_groups.pop()
required = parser.add_argument_group('required arguments')
required.add_argument("-i", type=str, action='store', metavar="input.fasta", help="input multiple sequence alignment in FASTA or A2M format", required=True)
parser.add_argument("--prehhfilter", action='store_true', help="Run  hhfilter BEFORE the main refinement, to make things faster")
parser.add_argument("-id", action='store', metavar=90, default= 90, help="Sequence Identity cut-off.  Default: 90%%")
parser.add_argument("-cov", action='store', metavar=75, default= 75, help="Alignment Coverage cut-off.  Default: 75%%")
parser.add_argument("-v", action='version', help="display version and exit")
parser._action_groups.append(help)
args = parser.parse_args()







def detect_gaps(seq):
    positions=[]
    for i in range(len(seq)):
        if seq[i]=="-" or seq[i]==".":
            positions.append(i)
    return positions

def remove_positions(seq, positions):
    new_seq=""
    for i in range(len(seq)):
        if i not in positions:
            new_seq+=seq[i]
    return new_seq


if args:
    inp=args.i
    cov=args.cov
    id=args.id
    run_hhfilter=args.prehhfilter
else:
    exit()







name=inp.split(".")
ext=name[-1]
prefix=inp.rstrip(".%s" %ext)


if run_hhfilter==False:
    if ext.lower() in ["a2m", "a3m"]:
        a3m_check==True
    else:
        a3m_check=False
    input_msa=inp
else:
    a3m_check=True
    tmp_filename="%s_tmp.a3m" %prefix
    print("Running preliminary hhfilter filtering...", file=sys.stderr)
    print("%s -i %s -o %s -id %s -cov %s -M 50" %(hhfilter, inp, tmp_filename, id, cov), file=sys.stderr)
    sp.call("%s -i %s -o %s -id %s -cov %s -M 50" %(hhfilter, inp, tmp_filename, id, cov), shell=True)
    input_msa=tmp_filename

print("Parsing alignment...", file=sys.stderr)
headers, sequences = parse_fasta(input_msa, a3m=a3m_check)

print("Finding pairwise distances...", file=sys.stderr)

idx_of_center_seq = pairwise_distances(np.array(AA_to_N(sequences)),metric="hamming").mean(0).argmin()


center_seq_header=headers[idx_of_center_seq]
center_seq_sequence=sequences[idx_of_center_seq]
pos_to_remove=detect_gaps(center_seq_sequence)
center_seq_sequence=remove_positions(center_seq_sequence, pos_to_remove)

reorder_seqs=""

print("Reorder sequences...", file=sys.stderr)

for i in range(len(sequences)):
    if headers[i]!=center_seq_header:
        seq=remove_positions(sequences[i], pos_to_remove)
        reorder_seqs+=">%s\n%s\n" %(headers[i], seq)

reorder_seqs=">%s\n%s\n%s" %(center_seq_header, center_seq_sequence, reorder_seqs)


f=open("tmp_%s" %inp, "w")
f.write(reorder_seqs)
f.close()

tmp_filename2="tmp_%s" %inp
output_filtered="%s_filtered.fasta" %prefix

if run_hhfilter==False:
    print("Run final hhfilter for refinement...", file=sys.stderr)
    print("%s -i %s -o %s -id %s -cov %s -M first" %(hhfilter, tmp_filename2, output_filtered, id, cov), file=sys.stderr)
    sp.call("%s -i %s -o %s -id %s -cov %s -M first" %(hhfilter, tmp_filename2, output_filtered, id, cov), shell=True)
else:
    sp.call("cp %s %s" %(tmp_filename2, output_filtered), shell=True)

#remove tmp file
if a3m_check == True:
    os.remove(tmp_filename)
os.remove(tmp_filename2)
