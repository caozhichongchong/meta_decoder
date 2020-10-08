#!/bin/bash
python PhaseFinder.py ratio -i /scratch/users/anniz44/ENIGMA/genome/SRX4123255_contigs.fasta.ID.fasta -1 /scratch/users/anniz44/ENIGMA/metagenomes/2old_S37_R1_001.fastq -2 /scratch/users/anniz44/ENIGMA/metagenomes/2old_S37_R2_001.fastq -p 16 -o /scratch/users/anniz44/ENIGMA/strain_finder/2old_S37_R1_001.fastq_SRX4123255_contigs.fasta.out
bwa index /scratch/users/anniz44/ENIGMA/genome/SRX4123255_contigs.fasta
bwa mem /scratch/users/anniz44/ENIGMA/genome/SRX4123255_contigs.fasta /scratch/users/anniz44/ENIGMA/metagenomes/2old_S37_R1_001.fastq /scratch/users/anniz44/ENIGMA/metagenomes/2old_S37_R2_001.fastq |samtools view -F 4 -Sb >/scratch/users/anniz44/ENIGMA/strain_finder/2old_S37_R1_001.fastq_SRX4123255_contigs.fasta.SRID.bam
python SRID.py -b /scratch/users/anniz44/ENIGMA/strain_finder/2old_S37_R1_001.fastq_SRX4123255_contigs.fasta.SRID.bam -p 12 -r 100 -m 200 -s 71 -n 4 -o /scratch/users/anniz44/ENIGMA/strain_finder/2old_S37_R1_001.fastq_SRX4123255_contigs.fasta.SRID.out.tab -t tmp
ls -l /scratch/users/anniz44/ENIGMA/strain_finder/2old_S37_R1_001.fastq_SRX4123255_contigs.fasta.SRID.bam
rm -rf /scratch/users/anniz44/ENIGMA/strain_finder/2old_S37_R1_001.fastq_SRX4123255_contigs.fasta.SRID.bam
