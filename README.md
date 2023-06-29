# functional_dark_matter_code


This repository contains script used for the manipulation of sequence alignments, preparation for 3D modeling and mass-scale 3D structure search of PDB models against reference datasets (e.g. PDB, SCOP) with TMalign/MMalign.

The aforementioned code was used in Pavlopoulos et al (2023), Unraveling the functional dark matter through global metagenomics. https://doi.org/10.1038/s41586-019-0000-0

### Contents:
- `generate_seed.py` : parses a multiple sequence alignment in FASTA format and produces its consensus sequence and a non-redundant dataset ("seed" alignment)

- `prepare_MSA_for_alphafold.py`: parses an alignment in FASTA or A3M format and produces a refined MSA output that can be used as input for AlphaFold2 or ColabFold / LocalColabFold.

- `structure_search.py`: performs structural alignment for a set of query PDB/mmCIF files against a set of reference PDB/mmCIF files with TMalign and/or MMalign.

### System Requirements:

- OS: Linux, MacOS or Windows 10/11 with WSL (Windows Subsystem for Linux) installed
- Python v. 3.8 or newer
- Python modules: numpy, scipy, scikit-learn, biopython, prody
- TMalign & MMalign
- HH-suite 3


### Installation instructions:
**Note:** Examples are shown for Ubuntu/Debian and similar Linux distributions. For rpm-based distributions, please use the equivalent package managers
1. Download and install software dependencies:
```
sudo apt install python3 python3-pip
sudo apt install hhsuite
sudo apt install tmalign
sudo apt install mmalign

pip3 install numpy scipy scikit-learn biopython prody

```

2.


### Example runs:

- Generation of seed alignments

```
cd example_data/MSAs/

python ../../generate_seed.py MSA_list.txt

```

- Prepare an alignment for 3D modeling with AlphaFold2:
```
cd example_data/MSAs/seed/
python ../../../prepare_MSA_for_alphafold.py --input F040820.fasta --hhfilter_bin /opt/hh-suite/bin/hhfilter --id 90 --cov 75

```


- Pairwise structural alignment of NMPF 3D models against SCOP domains

```
cd example_data/structures

python ../../structure_search.py --query NMPF_list.txt --target SCOP_list.txt --method tmalign --cpus 4 1> results_table.tsv 2>stderr.log

```
