# meta_decoder
# This is a temp version of meta_decoder
## Introduction
* meta_decoder: a reference-based metagenomics workflow
* input: metagenomes (-i)
* optional input: reference genomes (-r)
* requirement: bowtie samtools bedtools emboss  bwa bedops
* requirement: python 3 (branch py3); or python < 3 (branch master)
* requirement: python-pip python-pandas biopython python-numpy python-FuncDesigner python-DerApproximator
* to install all requirements:\
`apt-get install bowtie samtools bedtools emboss python-pip python-pandas`\
`wget https://github.com/bedops/bedops/releases/download/v2.4.35/bedops_linux_x86_64-v2.4.35.tar.bz2 && tar jxf bedops_linux_x86_64-v2.4.35.tar.bz2 && sudo cp bin/* /usr/local/bin/`\
`pip install openopt`\
`pip install biopython`\
`pip install pandas`\
`pip install numpy`\
`pip install FuncDesigner`\
`pip install DerApproximator`\
`git clone https://github.com/lh3/bwa.git && cd bwa && make`

## Install meta_decoder
`git clone https://github.com/caozhichongchong/meta_decoder.git`\
`git clone git@github.com:caozhichongchong/meta_decoder.git`\
in preparation: `pip install meta_decoder`\
in preparation: `anaconda download caozhichongchong/meta_decoder`

## Test (any of these two commands)
`python meta_decoder`\
compare the results in 'result_decoder' to the results in 'example'

## How to use it

`python meta_decoder -i your_dir_metagenomes -inf your_format_metagenomes`\
`python meta_decoder -i your_dir_metagenomes -inf your_format_metagenomes --r your_dir_reference_genomes --rf your_format_genomes`

## Results
The result dir of 'result_decoder':
### Phase variation
* `out.ratio.txt`: containing the reads that supporting either R or F orientation of invertible DNA
* `test.ID.fasta`: containing inverted (R) and non-inverted (F) putative invertible DNA regions flanked by sequences of specified length (bowtie indexed)
* `test.ID.fasta.info.tab`:  describing the location of inverted repeats in the above fasta file
* `test.ID.fasta.info.tab`:  describing the location of inverted repeats in the above fasta file
* `test.einverted.tab`:  containing the postion information of invereted repeats in the genome
### Mobile genetic elements and split read insertion detection (SRID)
* `out.tab`:  containing the Split Read Insertion Detection (SRID) in the genomes and metagenomes
### Strain finder
* `EM file`:  a binary cPickled file. This object holds: (1) the input alignment, (2) simulated data (a Data object), (3) the strain genotypes, and the strain frequencies (Estimate objects).
* `OTU table`:  This writes the strain genotypes and strain frequencies as an OTU table. The strain genotypes are included in the OTU names.

## Copyright
Copyright: An-Ni Zhang, Christopher Smillie, Xiaofang Jiang, Eric Alm, Alm Lab in MIT\
Contact: anniz44@mit.edu or caozhichongchong@gmail.com\
Citation:
1. Jiang X, Hall AB, et al. Invertible promoters mediate bacterial phase variation, antibiotic resistance, and host adaptation in the gut, Science (2019) DOI: 10.1126/science.aau5238\
2. Smillie, C. S., Sauk, J., Gevers, D., Friedman, J., Sung, J., Youngster, I., ... & Allegretti, J. R. (2018). Strain tracking reveals the determinants of bacterial engraftment in the human gut following fecal microbiota transplantation. Cell host & microbe, 23(2), 229-240.
