# meta_decoder
# This is a temp version of meta_decoder
## Introduction
* meta_decoder: a reference-based metagenomics workflow
* input: metagenomes (-i)
* optional input: reference genomes (-r)
* requirement: bowtie samtools bedtools vcftools emboss bwa bedops
* requirement: python 3 (branch py3); or python < 3 (branch master)
* requirement: python-pip python-pandas biopython python-numpy python-FuncDesigner python-DerApproximator
* to install all requirements:\
`apt-get install bowtie samtools bedtools emboss python-pip python-pandas`\
`wget https://github.com/bedops/bedops/releases/download/v2.4.35/bedops_linux_x86_64-v2.4.35.tar.bz2 && tar jxf bedops_linux_x86_64-v2.4.35.tar.bz2 && sudo cp bin/* /usr/local/bin/`\
`pip install vcfstats openopt biopython pandas numpy FuncDesigner DerApproximator`\
`git clone https://github.com/lh3/bwa.git && cd bwa && make`\
Require R and package ggplot2\
https://ggplot2.tidyverse.org/\
Unfortunately you may have to change the DerApproximator/__init__.py manually for python >= 3\
https://github.com/PythonCharmers/OOSuite/issues/1\
`from .DerApproximator import DerApproximatorException, get_d1, check_d1, get_d2`

## Install meta_decoder
`git clone https://github.com/caozhichongchong/meta_decoder.git`\
`git clone git@github.com:caozhichongchong/meta_decoder.git`\
in preparation: `pip install meta_decoder`\
in preparation: `anaconda download caozhichongchong/meta_decoder`

## Test (any of these two commands)
`python meta_decoder`\
compare the results in 'result_decoder' to the results in 'example'

## How to use it
#### Step1: alignment and statistics calculation
`python meta_decoder -i your_dir_metagenomes -inf your_format_metagenomes`\
`python meta_decoder -i your_dir_metagenomes -inf your_format_metagenomes --r your_dir_reference_genomes --rf your_format_genomes`
#### Step2 (optional): convert result into html formats
`python meta_decoder --html T`

## Statistics
* Measures nucleotide diversity 
1. on a per-site basis: output with the suffix ".sites.pi"
2. in windows (default of 1000bp): output with the suffix ".windowed.pi"
* Calculate allele frequency for each site: \
output with the suffix ".frq"
* Count raw allele counts for each site: \
output with the suffix ".frq.count"
* Calculate Transition / Transversion ratio in bins of size 1000bp:\
output with the suffix ".TsTv"
* Calculate a simple summary of all Transitions and Transversions:\
output with the suffix ".TsTv.summary"
* Calculates the Transition / Transversion ratio as a function of alternative allele count:\
output with the suffix ".TsTv.count"
* Calculates a measure of heterozygosity on a per-individual basis.\
Specfically, the inbreeding coefficient, F, is estimated for each individual using a method of moments:\
output with the suffix ".het"
* Reports a p-value for each site from a Hardy-Weinberg Equilibrium test (as defined by Wigginton, Cutler and Abecasis (2005)).\
The resulting file also contains the Observed numbers of Homozygotes and Heterozygotes and the corresponding Expected numbers under HWE:\
output with the suffix ".hwe"
* Calculate Tajima D statistic in bins with size of 1000bp:\
output with the suffix ".Tajima.D"
* Calculate a relatedness statistic of unadjusted Ajk statistic based on the method of Yang et al, Nature Genetics 2010 (doi:10.1038/ng.608).\
Specifically, calculate the unadjusted Ajk statistic.\
Expectation of Ajk is zero for individuals within a populations, and one for an individual with themselves:\
output with the suffix ".relatedness"
* Calculate number and density of SNPs in bins of size of 1000bp:\
output with the suffix ".snpden"
* Generate a histogram file of the length of all indels (including SNPs).\
It shows both the count and the percentage of all indels for indel lengths that occur at least once in the input file.\
SNPs are considered indels with length zero:\
output with the suffix ".indel.hist"

## Visualization in folder "plot"
* Alternative allele frequency on each contig: "Allele_frequency_on_each_contig.violin.png"
* Alternative allele frequency on each contig: "Allele_frequency_on_each_contig_(boxplot).boxplot.png"
* Counting type of mutant genotypes (HET, HOM_ALT) on each contig: "Mutant_genotypes_on_each_contig.col.png"
* Number of substitutions of SNPs pass all filters: "Number_of_substitutions_of_SNPs_(passed).col.png"
* Number of variants on each chromosome: "Number_of_variants_on_genome.col.png"
* Overall distribution of allele frequency: "Overall_allele_frequency_distribution.histogram.png"
* Counting types of variants on each contig: "Types_of_variants_on_each_contig.col.png"

## Results
The result dir of 'result_decoder':
### Phase variation
* `out.ratio.txt`: containing the metagenomic reads that supporting either R or F orientation of invertible DNA, in the format of \
`ID(contig:pos_A:pos_B:pos_C:pos_D) 	Pe_F(reads supporting F orientation)	Pe_R(reads supporting R orientation)	Pe_ratio(Pe_R/(Pe_F + Pe_R))	Span_F(reads supporting F orientation spanning the inverted repeat)	Span_R(reads supporting R orientation spanning the inverted repeat)	Span_ratio(Span_R/(Span_F + Span_R))`
* `test.ID.fasta`: containing inverted (R) and non-inverted (F) putative invertible DNA regions flanked by sequences of specified length (bowtie indexed)
* `test.ID.fasta.info.tab`:  describing the location of inverted repeats in the above fasta file, in the format of `contig:pos_A:pos_B:pos_C:pos_D`
* `test.einverted.tab`:  containing the position information of inverted repeats in the genome, in the format of \
`contig  Split_Site_1(The start of putative MGEs)   Split_Site_2(The end of putative MGEs)   Supporting_Split_Reads  Supporting_Read_Pairs`
### Mobile genetic elements and split read insertion detection (SRID)
* `out.tab`:  containing the Split Read Insertion Detection (SRID) in the genomes and metagenomes, in the format of `Scaffolds  pos_A   pos_B   pos_C   pos_D   IR  Sequence_within    IR`
### Strain finder
* `EM file`:  a binary cPickled file. This object holds: (1) the input alignment, (2) simulated data (a Data object), (3) the strain genotypes, and the strain frequencies (Estimate objects).
* `OTU table`:  This writes the strain genotypes and strain frequencies as an OTU table. The strain genotypes are included in the OTU names.

## Copyright
Copyright: An-Ni Zhang, Christopher Smillie, Xiaofang Jiang, Eric Alm, Alm Lab in MIT\
https://github.com/cssmillie/StrainFinder \
https://github.com/XiaofangJ/PhaseFinder \
https://github.com/XiaofangJ/SRID \
Contact: anniz44@mit.edu or caozhichongchong@gmail.com\
Citation:
1. Jiang X, Hall AB, et al. Invertible promoters mediate bacterial phase variation, antibiotic resistance, and host adaptation in the gut, Science (2019) DOI: 10.1126/science.aau5238
2. Smillie, C. S., Sauk, J., Gevers, D., Friedman, J., Sung, J., Youngster, I., ... & Allegretti, J. R. (2018). Strain tracking reveals the determinants of bacterial engraftment in the human gut following fecal microbiota transplantation. Cell host & microbe, 23(2), 229-240.
3. Xiaofang Jiang, Andrew Brantley Hall, Ramnik J Xavier, Eric J Alm (2016) Comprehensive analysis of mobile genetic elements in the gut microbiome reveals a phylum-level niche-adaptive gene pool bioRxiv 214213; https://www.biorxiv.org/content/10.1101/214213v2.
