from Bio import SeqIO
import glob
import os
import argparse
import pandas as pd


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
parser.add_argument('--bwa',
                    help="Optional: complete path to bwa if not in PATH,",
                    metavar="/usr/local/bin/bwa",
                    action='store', default='bwa', type=str)
parser.add_argument('--bcf',
                    help="Optional: complete path to bcftools if not in PATH,",
                    metavar="/usr/local/bin/bcftools",
                    action='store', default='bcftools', type=str)
parser.add_argument('--sam',
                    help="Optional: complete path to bwa if not in PATH,",
                    metavar="/usr/local/bin/samtools",
                    action='store', default='samtools', type=str)
parser.add_argument('--vcf',
                    help="Optional: complete path to bwa if not in PATH,",
                    metavar="/usr/local/bin/vcftools",
                    action='store', default='vcftools', type=str)
parser.add_argument("--html",
                    help="convert output into html (--html T)",
                    type=str, default='F',
                    metavar='F or T')

################################################## Definition ########################################################
args = parser.parse_args()
try:
    os.mkdir(args.o)
except OSError:
    pass
################################################## Function ########################################################
def reverse_metagenome(filename):
    f1=open(filename.split('1'+args.inf)[0].split(args.inf)[0]+'2'+args.inf,'w')
    f1.close()
    return cmd

def bowtie(genomes, metagenomes):
    # bowtie alignment
    cmds = ''
    tempbamoutput = os.path.join(args.o, os.path.split(metagenomes)[-1] + '_' + os.path.split(genomes)[-1])
    try:
        ftest = open(genomes + '.bwt', 'r')
    except IOError:
        # build bwa library
        cmds += args.bwa + ' index %s \n' % (genomes)
    # run alignment
    cmds += args.bwa + ' mem %s %s |%s view -S -b >%s.bam\n%s sort %s.bam -o %s.sorted.bam\n%s index %s.sorted.bam\n' % (
        genomes, metagenomes,args.sam,
        tempbamoutput, args.sam, tempbamoutput, tempbamoutput, args.sam, tempbamoutput)
    cmds += '%s mpileup -Ou -f %s %s.sorted.bam  | %s call -mv > %s.raw.vcf' % (
        args.bcf,genomes, tempbamoutput, args.bcf,tempbamoutput)
    cmds += '\n%s filter -s LowQual -e \'%s || DP>100\' %s.raw.vcf > %s.flt.vcf \n' % (
        args.bcf,'QUAL<20', tempbamoutput, tempbamoutput)
    # statistics
    # output nucleotide_diversity per site with the suffix ".sites.pi"
    cmds += '%s --vcf %s.flt.vcf --site-pi --out %s.flt.vcf \n' % (
        args.vcf, tempbamoutput, tempbamoutput)
    # output nucleotide_diversity per 1000bp with the suffix ".windowed.pi"
    cmds += '%s --vcf %s.flt.vcf --window-pi 1000 --out %s.flt.vcf \n' % (
        args.vcf, tempbamoutput, tempbamoutput)
    # output allele frequency for each site with the suffix ".frq"
    cmds += '%s --vcf %s.flt.vcf --freq --out %s.flt.vcf \n' % (
        args.vcf, tempbamoutput, tempbamoutput)
    # output raw allele counts for each site with the suffix ".frq.count"
    cmds += '%s --vcf %s.flt.vcf --counts --out %s.flt.vcf \n' % (
        args.vcf, tempbamoutput, tempbamoutput)
    # output Transition / Transversion ratio  in bins of size 1000bp with the suffix ".TsTv"
    cmds += '%s --vcf %s.flt.vcf --TsTv 1000 --out %s.flt.vcf \n' % (
        args.vcf, tempbamoutput, tempbamoutput)
    # output a simple summary of all Transitions and Transversions with the suffix ".TsTv.summary"
    cmds += '%s --vcf %s.flt.vcf --TsTv-summary --out %s.flt.vcf \n' % (
        args.vcf, tempbamoutput, tempbamoutput)
    # output Transition / Transversion ratio as a function of alternative allele count with the suffix ".TsTv.count"
    cmds += '%s --vcf %s.flt.vcf --TsTv-by-count --out %s.flt.vcf \n' % (
        args.vcf, tempbamoutput, tempbamoutput)
    # output a measure of heterozygosity on a per-individual basis with the suffix ".het"
    cmds += '%s --vcf %s.flt.vcf --het --out %s.flt.vcf \n' % (
        args.vcf, tempbamoutput, tempbamoutput)
    # output p-value for each site from a Hardy-Weinberg Equilibrium test with the suffix ".hwe"
    cmds += '%s --vcf %s.flt.vcf --hardy --out %s.flt.vcf \n' % (
        args.vcf, tempbamoutput, tempbamoutput)
    # output Tajima D statistic in bins with size of 1000bp with the suffix ".Tajima.D"
    cmds += '%s --vcf %s.flt.vcf --TajimaD 1000 --out %s.flt.vcf \n' % (
        args.vcf, tempbamoutput, tempbamoutput)
    # output a relatedness statistic of unadjusted Ajk statistic with the suffix ".relatedness"
    # Expectation of Ajk is zero for individuals within a populations, and one for an individual with themselves
    cmds += '%s --vcf %s.flt.vcf --relatedness --out %s.flt.vcf \n' % (
        args.vcf, tempbamoutput, tempbamoutput)
    # output number and density of SNPs in bins of size of 1000bp with the suffix ".snpden"
    cmds += '%s --vcf %s.flt.vcf --SNPdensity 1000 --out %s.flt.vcf \n' % (
        args.vcf, tempbamoutput, tempbamoutput)
    # output a histogram file of the length of all indels (including SNPs) with the suffix ".indel.hist"
    cmds += '%s --vcf %s.flt.vcf --hist-indel-len --out %s.flt.vcf \n' % (
        args.vcf, tempbamoutput, tempbamoutput)
    return cmds

def tsv_to_html(file):
    # Read the csv file in
    df = pd.read_csv(file, sep='\t',header=0,error_bad_lines=False)
    # Save to file
    df.to_html(file+'.html')

################################################## Programme ########################################################
if args.html == 'F':
    # process alignment and statistic
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
    #f0.write(cmd)

    # generate metagenome list file
    metagenome_files=glob.glob(os.path.join(args.i,'*'+args.inf))
    f1=open(os.path.join(args.o,'metagenome.list'),'w')
    cmd = ''
    for filename in metagenome_files:
            # create reverse orienation metagenomes
            try:
                ftest = open(filename.split('1'+args.inf)[0].split('2'+args.inf)[0]+'2'+args.inf,'r')
                f1.write(str(filename) + '\n')
            except IOError:
                cmd += reverse_metagenome(filename)
                f1.write(filename.split(args.inf)[0]+'2'+args.inf+'\n')
                f1.write(filename.split(args.inf)[0] + '1' + args.inf + '\n')
    f1.close()
    f0.write(cmd)

    # run strain finder preprocess
    cmd = 'python bin/0.run.py --fastqs %s  --ref %s  --map %s\n' % (
                os.path.join(args.o, 'metagenome.list'),
                os.path.join(args.o,'all.ref.genomes.fasta'),
                os.path.join(args.o, 'ref.map.txt'))
    f0.write(cmd)

    # run PhaseFinder
    i=1
    for genomes in genome_files:
        if '.ID.fasta' not in genomes:
            try:
                f1=open(genomes+'.ID.fasta','r')
            except IOError:
                #fsub = open(str(int(i % task)) + '.sh', 'a')
                #fsub.write('#!/bin/bash\n')
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
                #fsub.write(cmd)
                #fsub.close()

    # run bowtie, PhaseFinder and SRID
    i=1
    for genomes in genome_files:
        if '.ID.fasta' not in genomes:
            for metagenomes in metagenome_files:
                if '1'+args.inf in metagenomes or args.s == 1:
                    fsub = open(str(int(i % task)) + '.sh', 'a')
                    fsub.write('#!/bin/bash\n')
                    cmd = ''
                    try:
                        ftest=open(os.path.join(args.o,os.path.split(metagenomes)[-1]+'_'+os.path.split(genomes)[-1]+'.out'),'r')
                    except IOError:
                        cmd += ''
                        #cmd += 'python bin/PhaseFinder.py ratio -i %s -1 %s -2 %s -p 16 -o %s\n' % (
                        #    genomes+'.ID.fasta', metagenomes,
                        #    metagenomes.replace('1'+args.inf,'2'+args.inf),
                        #    os.path.join(args.o,os.path.split(metagenomes)[-1]+'_'+os.path.split(genomes)[-1]+'.out'))
                    cmd += bowtie(genomes, metagenomes)
                    metagenomes = metagenomes.replace('1' + args.inf, '2' + args.inf)
                    cmd += bowtie(genomes, metagenomes)
                    #cmd += 'python bin/SRID.py -b %s -p 12 -r 100 -m 200 -s 71 -n 4 -o %s -t tmp\n' % (
                    #    os.path.join(args.o, os.path.split(metagenomes)[-1] + '_' + os.path.split(genomes)[-1] + '.bam'),
                    #    os.path.join(args.o, os.path.split(metagenomes)[-1] + '_' + os.path.split(genomes)[-1] + '.SRID.out.tab'))
                    #cmd += '#ls -l %s\n#rm -rf %s\n' %(os.path.join(args.o, os.path.split(metagenomes)[-1] + '_' +
                    #       os.path.split(genomes)[-1] + '.SRID.bam'),
                    #                                 os.path.join(args.o,
                    #                                              os.path.split(metagenomes)[-1] + '_' + os.path.split(genomes)[
                    #                                                  -1] + '.SRID.bam'))
                    i+=1
                    fsub.write(cmd)
                    fsub.close()

    # run strain finder
    i=1
    for genomes in genome_files:
        if '.ID.fasta' not in genomes:
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
            f0.write(('nohup sh %s > %s.nohup.out&\n')%(files,files))
            #f0.write('sbatch -p sched_mem1TB -c 40 -t 5-00:00:00 --mem=500000 -J %straits -o %s.out -e %s.err %s\n' %(files,files,files,files))
    f0.close()
else:
    # convert results to html files
    output_files=glob.glob(os.path.join(args.o,'*.flt.vcf.*'))
    for output_filename in output_files:
        if '.html' not in output_filename:
            tsv_to_html(output_filename)