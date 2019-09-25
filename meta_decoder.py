from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import glob
import os
import argparse

############################################ Arguments and declarations ##############################################
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("-i",
                    help="input dir of metagenomes", type=str,
                    default='example',metavar='current dir (.)')
parser.add_argument("-inf",
                    help="input file format",
                    type=str, default='.fq', metavar='.fastq')
parser.add_argument("--r",
                    help="input dir of reference genomes", type=str,
                    default='example',metavar='current dir (.)')
parser.add_argument("--rf",
                    help="reference genome file format",
                    type=str, default='.fa', metavar='.fasta')
parser.add_argument("--s",
                    help="1: single orienation or 2: both orienations",
                    type=int, default=2, metavar='1 or 2', choices = [1,2])
parser.add_argument("--o",
                    help="output directory",
                    type=str, default='result_decoder',
                    metavar='result_decoder')
parser.add_argument("--t",
                    help="number of threads you would like to use",
                    type=int, default=1,
                    metavar='1 or more')

################################################## Definition ########################################################
args = parser.parse_args()
try:
    os.mkdir(args.o)
except OSError:
    pass
################################################## Function ########################################################
def reverse_metagenome(filename):
    cmd = 'mv %s %s\n' %(filename,filename.split(args.inf)[0]+'1'+args.inf)
    f1=open(filename.split(args.inf)[0]+'2'+args.inf,'w')
    f1.close()
    return cmd

################################################## Programme ########################################################
f0=open('meta.decoder.sh','w')
f0.write('#!/bin/bash\n')
# number of bash scripts
task = int(args.t)
# generate map file
genome_files=glob.glob(os.path.join(args.r,'*'+args.rf))
f1=open(os.path.join(args.o,'ref.map.txt'),'w')
for filename in genome_files:
    #if '.ID.fasta' not in filename:
        for record in SeqIO.parse(filename, "fasta"):
            f1.write(os.path.split(filename)[-1]+'\t'+str(record.id)+'\n')
f1.close()
cmd = 'cat '+os.path.join(args.r,'*'+args.rf)+'> %s\n' %(os.path.join(args.o,'all.ref.genomes.fasta'))
f0.write(cmd)

# generate metagenome list file
metagenome_files=glob.glob(os.path.join(args.i,'*'+args.inf))
f1=open(os.path.join(args.o,'metagenome.list'),'w')
cmd = ''
for filename in metagenome_files:
    if args.s == 1:
        # create reverse orienation metagenomes
        cmd += reverse_metagenome(filename)
        f1.write(filename.split(args.inf)[0]+'2'+args.inf+'\n')
        f1.write(filename.split(args.inf)[0] + '1' + args.inf + '\n')
    else:
        f1.write(str(filename) + '\n')
f1.close()
f0.write(cmd)

# run strain finder preprocess

cmd = 'python bin/0.run.py --fastqs %s  --ref %s  --map %s\n' % (
            os.path.join(args.o, 'metagenome.list'), os.path.join(args.o,'all.ref.genomes.fasta'),
            os.path.join(args.o, 'ref.map.txt'))
f0.write(cmd)
# run PhaseFinder

i=1
for genomes in genome_files:
    if '.ID.fasta' not in genomes:
        try:
            f1=open(genomes+'.ID.fasta','r')
        except IOError:
            fsub = open(str(int(i % task)) + '.sh', 'a')
            fsub.write('#!/bin/bash\n')
            cmd = ''
            try:
                ftest = open(
                    genomes+'.einverted.tab', 'r')
            except IOError:
                cmd += 'python bin/PhaseFinder.py locate -f %s -t %s -g 15 85 -p\n'\
                       %(genomes,genomes+'.einverted.tab')
            try:
                ftest = open(
                    genomes + '.ID.fasta', 'r')
            except IOError:
                cmd += 'python bin/PhaseFinder.py create -f %s -t %s -s 1000 -i %s\n' \
                       %(genomes,genomes+'.einverted.tab',genomes+'.ID.fasta')
            i+=1
            fsub.write(cmd)
            fsub.close()

# run PhaseFinder and SRID
i=1
for genomes in genome_files:
    if '.ID.fasta' not in genomes:
        for metagenomes in metagenome_files:
            if '1'+args.inf in metagenomes or args.s == 1:
                if args.s == 1:
                    metagenomes = metagenomes.split(args.inf)[0]+'1'+args.inf
                fsub = open(str(int(i % task)) + '.sh', 'a')
                fsub.write('#!/bin/bash\n')
                cmd = ''
                try:
                    ftest=open(os.path.join(args.o,os.path.split(metagenomes)[-1]+'_'+os.path.split(genomes)[-1]+'.out'),'r')
                except IOError:
                    cmd += 'python bin/PhaseFinder.py ratio -i %s -1 %s -2 %s -p 16 -o %s\n' % (
                        genomes+'.ID.fasta', metagenomes,
                        metagenomes.replace('1'+args.inf,'2'+args.inf),
                        os.path.join(args.o,os.path.split(metagenomes)[-1]+'_'+os.path.split(genomes)[-1]+'.out'))
                cmd += 'bwa index %s \nbwa mem %s %s %s |samtools view -F 4 -Sb >%s\n'% (
                    genomes,
                    genomes, metagenomes,
                    metagenomes.replace('1'+args.inf,'2'+args.inf),
                    os.path.join(args.o,os.path.split(metagenomes)[-1]+'_'+os.path.split(genomes)[-1]+'.SRID.bam'))
                cmd += 'python bin/SRID.py -b %s -p 12 -r 100 -m 200 -s 71 -n 4 -o %s -t tmp\n' % (
                    os.path.join(args.o, os.path.split(metagenomes)[-1] + '_' + os.path.split(genomes)[-1] + '.SRID.bam'),
                    os.path.join(args.o, os.path.split(metagenomes)[-1] + '_' + os.path.split(genomes)[-1] + '.SRID.out.tab'))
                cmd += '#ls -l %s\n#rm -rf %s\n' %(os.path.join(args.o, os.path.split(metagenomes)[-1] + '_' +
                       os.path.split(genomes)[-1] + '.SRID.bam'),
                                                 os.path.join(args.o,
                                                              os.path.split(metagenomes)[-1] + '_' + os.path.split(genomes)[
                                                                  -1] + '.SRID.bam'))
                i+=1
                fsub.write(cmd)
                fsub.close()

# run strain finder
i=1
for genomes in genome_files:
    #if '.ID.fasta' not in genomes:
        fsub = open(str(int(i % task)) + '.sh', 'a')
        fsub.write('#!/bin/bash\n')
        genomes = os.path.split(genomes)[-1]
        genome_alignments = genomes +'.np.cPickle'
        cmd = 'python bin/StrainFinder.py --aln %s  -N 5 --max_reps 10 --dtol 1 --ntol 2 --max_time 3600 --converge --em %s.cpickle'%(genome_alignments, genomes)+\
               ' --em_out %s.cpickle --otu_out %s.otu_table.txt --log %s.log.txt --n_keep 3 --force_update --merge_out --msg\n' %(genomes, genomes,genomes)
        cmd += 'mv *.cpickle *.otu_table.txt *.log.txt > '+str(args.o)+'\n'
        cmd += 'mv %s > %s \n' %(args.r+'/*'+args.rf+'.*',args.o)
        i += 1
        fsub.write(cmd)
        fsub.close()

shfiles = glob.glob('*.sh')
for files in shfiles:
    if 'meta.decoder.sh' not in files:
        #f0.write(('nohup sh %s > %s.nohup.out&\n')%(files,files))
        f0.write('sbatch -p sched_mem1TB -c 40 -t 5-00:00:00 --mem=500000 -J %straits -o %s.out -e %s.err %s\n' %(files,files,files,files))
f0.close()
