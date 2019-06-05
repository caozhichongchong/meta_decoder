#!/bin/bash
cat /scratch/users/anniz44/ENIGMA/genome/*.fasta> /scratch/users/anniz44/ENIGMA/strain_finder/all.ref.genomes.fasta
python 0.run.py --fastqs /scratch/users/anniz44/ENIGMA/strain_finder/metagenome.list  --ref /scratch/users/anniz44/ENIGMA/strain_finder/all.ref.genomes.fasta  --map /scratch/users/anniz44/ENIGMA/strain_finder/ref.map.txt
python StrainFinder.py --aln.cpickle  -N 5 --max_reps 10 --dtol 1 --ntol 2 --max_time 3600 --converge --em em.cpickle --em_out em.cpickle --otu_out otu_table.txt --log log.txt --n_keep 3 --force_update --merge_out --msg
mv em.cpickle otu_table.txt log.txt > /scratch/users/anniz44/ENIGMA/strain_finder
mv /scratch/users/anniz44/ENIGMA/genome/*.fasta.* > /scratch/users/anniz44/ENIGMA/strain_finder
