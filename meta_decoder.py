from Bio import SeqIO
import glob
import os
import argparse
import pandas as pd
import subprocess
import numpy as np
import random


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
parser.add_argument("--os",
                    help="output directory for scripts",
                    type=str, default='sub_metadecoder',
                    metavar='single_cell')
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
parser.add_argument('--vcfstats',
                    help="Optional: complete path to vcfstats if not in PATH,",
                    metavar="/usr/local/bin/vcfstats",
                    action='store', default='vcfstats', type=str)
parser.add_argument("--html",
                    help="convert output into html (--html T)",
                    type=str, default='F',
                    metavar='F or T')
parser.add_argument("--cal",
                    help="calculate coverage (--cal T)",
                    type=str, default='F',
                    metavar='F or T')
parser.add_argument('--strainfinder',
                    help="Optional: complete path to strainfinder",
                    metavar="/scratch/users/anniz44/bin/miniconda3/bin/strainfinder",
                    action='store',
                    default='None',
                    type=str)


################################################## Definition ########################################################
args = parser.parse_args()
try:
    os.mkdir(args.o)
except OSError:
    pass
try:
    os.mkdir(args.os)
except OSError:
    pass
try:
    os.mkdir(os.path.join(args.o,'plots'))
except OSError:
    pass
SNP_outputdir = os.path.join(args.o,'snp_dynamics')
try:
    os.mkdir(SNP_outputdir)
except OSError:
    pass

# Set up A T G C
Allels = dict()
Allels['A']=0
Allels['T']=1
Allels['G']=2
Allels['C']=3
Allels_set = 'ATGC'
# import strainfinder
import os, sys
if args.strainfinder!='None':
    sys.path.append(args.strainfinder)
    from strainFinder import fitStrains, genomes_given_fracs


################################################## Function ########################################################
def reverse_metagenome(filename):
    f1=open(filename.split('1'+args.inf)[0].split(args.inf)[0]+'2'+args.inf,'w')
    f1.close()
    return cmd

def visual(vcffile,outputdir):
    suboutputdir = os.path.join(args.o, 'plots/'+os.path.split(vcffile)[1].split('.flt.vcf')[0])
    try:
        os.mkdir(suboutputdir)
    except OSError:
        pass
    cmds = ''
    cmds += "%s --vcf %s --outdir %s/ --formula \'COUNT(1) ~ CONTIG\' --title \'Number of variants on genome\' --ggs \'ylab(\"# Variants\")\'" \
            %(args.vcfstats,vcffile,suboutputdir)
    cmds += "\n%s --vcf %s --outdir %s/ --formula \'COUNT(1, VARTYPE[snp]) ~ SUBST[A>T,A>G,A>C,T>A,T>G,T>C,G>A,G>T,G>C,C>A,C>T,C>G]\' --title \'Number of substitutions of SNPs (passed)\' --passed --ggs \'ylab(\"# Substitutions of SNPs\")\'"\
            % (args.vcfstats, vcffile, suboutputdir)
    cmds += "\n%s --vcf %s --outdir %s/ --formula \'AAF ~ CONTIG\' --title \'Allele frequency on each contig\' --ggs \'ylab(\"# Allele frequency\")\'" \
            % (args.vcfstats, vcffile, suboutputdir)
    cmds += "\n%s --vcf %s --outdir %s/ --formula \'AAF ~ CONTIG\' --title \'Allele frequency on each contig (boxplot)\' --figtype boxplot --ggs \'ylab(\"# Allele frequency\")\'" \
            % (args.vcfstats, vcffile, suboutputdir)
    cmds += "\n%s --vcf %s --outdir %s/ --formula \'AAF ~ 1\' --title \'Overall allele frequency distribution\' --ggs \'ylab(\"Density\")\'" \
            % (args.vcfstats, vcffile, suboutputdir)
    cmds += "\n%s --vcf %s --outdir %s/ --formula \'COUNT(1, group=VARTYPE) ~ CHROM\' --title \'Types of variants on each contig\' --ggs \'ylab(\"Percentage\")\'" \
            % (args.vcfstats, vcffile, suboutputdir)
    cmds += "\n%s --vcf %s --outdir %s/ --formula \'COUNT(1, group=GTTYPEs[HET,HOM_ALT]{0}) ~ CHROM\' --title \'Mutant genotypes on each contig\' --ggs \'ylab(\"# Variants\")\'\n" \
            % (args.vcfstats, vcffile, suboutputdir)
    return cmds

def statistics_vcf(vcffile):
    cmds = ''
    # output nucleotide_diversity per site with the suffix ".sites.pi"
    cmds += '%s --vcf %s --site-pi --out %s \n' % (
        args.vcf,  vcffile,  vcffile)
    # output nucleotide_diversity per 1000bp with the suffix ".windowed.pi"
    cmds += '%s --vcf %s --window-pi 1000 --out %s \n' % (
        args.vcf,  vcffile,  vcffile)
    # output allele frequency for each site with the suffix ".frq"
    cmds += '%s --vcf %s --freq --out %s \n' % (
        args.vcf,  vcffile,  vcffile)
    # output raw allele counts for each site with the suffix ".frq.count"
    cmds += '%s --vcf %s --counts --out %s \n' % (
        args.vcf,  vcffile,  vcffile)
    # output Transition / Transversion ratio  in bins of size 1000bp with the suffix ".TsTv"
    cmds += '#%s --vcf %s --TsTv 1000 --out %s \n' % (
        args.vcf,  vcffile,  vcffile)
    # output a simple summary of all Transitions and Transversions with the suffix ".TsTv.summary"
    cmds += '#%s --vcf %s --TsTv-summary --out %s \n' % (
        args.vcf,  vcffile,  vcffile)
    # output Transition / Transversion ratio as a function of alternative allele count with the suffix ".TsTv.count"
    cmds += '#%s --vcf %s --TsTv-by-count --out %s \n' % (
        args.vcf,  vcffile,  vcffile)
    # output a measure of heterozygosity on a per-individual basis with the suffix ".het"
    cmds += '#%s --vcf %s --het --out %s \n' % (
        args.vcf,  vcffile,  vcffile)
    # output p-value for each site from a Hardy-Weinberg Equilibrium test with the suffix ".hwe"
    cmds += '#%s --vcf %s --hardy --out %s \n' % (
        args.vcf,  vcffile,  vcffile)
    # output Tajima D statistic in bins with size of 1000bp with the suffix ".Tajima.D"
    cmds += '#%s --vcf %s --TajimaD 1000 --out %s \n' % (
        args.vcf,  vcffile,  vcffile)
    # output a relatedness statistic of unadjusted Ajk statistic with the suffix ".relatedness"
    # Expectation of Ajk is zero for individuals within a populations, and one for an individual with themselves
    cmds += '#%s --vcf %s --relatedness --out %s \n' % (
        args.vcf,  vcffile,  vcffile)
    # output number and density of SNPs in bins of size of 1000bp with the suffix ".snpden"
    cmds += '%s --vcf %s --SNPdensity 1000 --out %s \n' % (
        args.vcf,  vcffile,  vcffile)
    # output a histogram file of the length of all indels (including SNPs) with the suffix ".indel.hist"
    cmds += '#%s --vcf %s --hist-indel-len --out %s \n' % (
        args.vcf,  vcffile,  vcffile)
    return cmds

def bowtie(genomes, metagenomes):
    print('processing bowtie mapping %s for %s' % (metagenomes, genomes))
    # bowtie alignment
    cmds = ''
    tempbamoutput = os.path.join(args.o, os.path.split(metagenomes)[-1] + '_' + os.path.split(genomes)[-1])
    try:
        ftest = open(genomes + '.bwt', 'r')
    except IOError:
        # build bwa library
        os.system(args.bwa + ' index %s \n' % (genomes))
    # run alignment
    try:
        ftest = open('%s.sorted.bam' %(tempbamoutput),'r')
    except IOError:
        cmds += args.bwa + ' mem -t %s %s %s |%s view -@ %s -S -b -F 4 >%s.bam\n%s sort -@ %s %s.bam -o %s.sorted.bam\n%s index -@ %s %s.sorted.bam\n' % (
            args.t, genomes, metagenomes, args.sam, args.t,
            tempbamoutput, args.sam, args.t, tempbamoutput, tempbamoutput, args.sam, args.t, tempbamoutput)
    cmds += '%s mpileup --threads %s -q30 -B -Ou -d3000 -f %s %s.sorted.bam  | %s call --ploidy 1 --threads %s -mv > %s.raw.vcf' % (
        args.bcf, args.t, genomes, tempbamoutput, args.bcf,args.t, tempbamoutput)
    cmds += '\n%s filter --threads %s -s LowQual -e \'DP>100\' %s.raw.vcf > %s.flt.vcf \n' % (
        args.bcf,args.t, tempbamoutput, tempbamoutput)
    # calculate coverage
    cmds += '%s depth -Q 10 %s.sorted.bam > %s.sorted.bam.cov\n' % (
        args.sam, tempbamoutput, tempbamoutput)
    cmds += 'echo -e "Ref_ID\\tCov_length\\tAverage\\tStdev" > %s.sorted.bam.avgcov\n' % (tempbamoutput)
    cmds += '%s depth %s.sorted.bam |  awk \'{sum[$1]+=$3; sumsq[$1]+=$3*$3; count[$1]++} END { for (id in sum) { print id,"\t",count[id],"\t",sum[id]/count[id],"\t",sqrt(sumsq[id]/count[id] - (sum[id]/count[id])**2)}}\' >> %s.sorted.bam.avgcov\n' % (
        args.sam, tempbamoutput, tempbamoutput)
    # statistics
    #cmds += statistics_vcf('%s.flt.vcf' % tempbamoutput)
    # optional cleanup
    cmds += 'rm -rf %s.bam\n' % (tempbamoutput)
    return cmds


def tsv_to_html(file):
    # Read the csv file in
    df = pd.read_csv(file, sep='\t',header=0,error_bad_lines=False)
    # Save to file
    df.to_html(file+'.html')


def cal_cov(output_filename):
    os.system('cat %s | cut -f 3 | sort -n > %s.sort' %(output_filename,output_filename))
    #print("sed \"s/.*\\t//g\" %s | sort - n | head - n \"$(($(wc -l %s | sed \'s/ .*//g\')*9/10))\" | > awk \'{total = total + $1}END{print total}\' > %s.cut"
    #      %(output_filename,output_filename,output_filename))
    Coverage = 0
    Total_bp = len(open('%s.sort' %(output_filename)).readlines())
    Lines_num = 0
    Genome_length = 3000000
    genomes = output_filename.split(args.inf+'_')[1].split(args.rf)[0]
    for lines in open('%s.sort' %(output_filename),'r'):
        Lines_num += 1
        new_cov = int(lines.replace('\r','').replace('\n',''))
        if Lines_num <= float(Total_bp)*0.9 and new_cov >=2:
            Coverage += new_cov
    for genomes in Length:
        if genomes in output_filename:
            Genome_length = Length[genomes]*0.9
            break
    os.system('rm -rf %s.sort' %output_filename)
    all_output.write('%s\t%s\t%s\t%s\n'%(output_filename.split(args.inf)[0],genomes,
                     '%.3f'%(float(Lines_num/Genome_length)),
                                         '%.3f'%(float(Coverage/Genome_length))))


def snpsum(output_filename):
    SNP = 0
    genomes = output_filename.split(args.inf+'_')[1].split(args.rf)[0]
    Lines_num = 0
    for lines in open((output_filename),'r'):
        Lines_num += 1
        try:
            SNP += int(lines.split('\t')[2])
        except ValueError:
            pass
    if Lines_num > 1:
        all_output.write('%s\t%s\t%s\t%s\n' % (output_filename.split(args.inf)[0], genomes,
                                               '%.3f' % (float(SNP / Lines_num)),
                                               int(SNP)))


def strain_finder(SNP_file):
    try:
        foutput = open(SNP_file + '.abu', 'r')
    except IOError:
        try:
            counts = np.genfromtxt(SNP_file, int, skip_header=1)
            foutput = open(SNP_file + '.abu','w')
            ## fit strains and relative abundances
            fracs, e, ll, dof = fitStrains(counts)
            best_fracs = list(fracs.values())[-2]
            best_fracs_new = []
            for abu in best_fracs:
                best_fracs_new.append(str('%.3f'%(abu)))
            foutput.write('%s\n' % ('\t'.join(best_fracs_new)))
            ## get ML genomes
            n = len(best_fracs)  # fitted number of strains
            k = counts.shape[1]  # number of alleles in alphabet
            perm_file = os.path.join(args.strainfinder, 'presence_matrices/strains_%d.alleles_%d.npy' % (n, k))
            allele_perm = np.load(perm_file)
            genomes = genomes_given_fracs(counts, best_fracs, e[n], alleles=['A', 'C', 'G', 'T'])
            ## output genomes
            genome_file = SNP_file + '.genome'
            try:
                np.savetxt(genome_file, genomes, '%s', '\t')
            except IOError:
                pass
        except IndexError:
            pass



def freq_call_sub(vcf_file_list):
    SNP = dict()
    Position_diff = []
    Chr = ''
    Position = 0
    for vcf_file in vcf_file_list:
        for lines in open(vcf_file, 'r'):
            if not lines.startswith("#"):
                lines_set = lines.split('\t')
                Chr_new = lines_set[0]
                Position_new = int(lines_set[1])
                Chr_position = '%s_%s' % (Chr_new, Position_new)
                Allels_frq = [0, 0, 0, 0]
                SNP.setdefault(Chr_position, Allels_frq)
                Allels_frq = SNP[Chr_position]
                if "INDEL" not in lines_set[7] and lines_set[6] != 'LowQual':
                    # calculate position apart
                    if Chr_new == Chr:
                        if Position_new != Position:
                            Position_diff.append(abs(Position_new - Position) % 3)
                    Chr = Chr_new
                    Position = Position_new
                    Poly = lines_set[7].split('DP4=')[1].split(';')[0].split(',')
                    # Depth = int(lines_set[7].split('DP=')[1].split(';')[0])
                    Sum_ref = int(Poly[0]) + int(Poly[1])
                    Sum_alt = int(Poly[2]) + int(Poly[3])
                    Allels_frq[Allels[lines_set[3]]] += Sum_ref
                    try:
                        Allels_frq[Allels[lines_set[4]]] += int(Sum_alt)
                    except KeyError:
                        allels_set = lines_set[4].split(',')
                        for allels in allels_set:
                            if allels in Allels:
                                Allels_frq[Allels[allels]] += int(Sum_alt * 0.5)
                            else:
                                pass
                    SNP[Chr_position] = Allels_frq
    return [SNP,Position_diff]

def freq_call(vcf_file,cov_file):
    Output = []
    Output2 = []
    Output3 = []
    try:
        foutput = open(cov_file + '.frq.clean.snp', 'r')
    except IOError:
        try:
            SNP, Position_diff = freq_call_sub([vcf_file])
            for Chr_position in SNP:
                Depth = sum(SNP[Chr_position])
                if Depth > 2 and Depth <= 3000 and not any(Depth == allels for allels in SNP[Chr_position]):
                    Output3.append('%s\t%s\t%s\t%s\n' % (SNP[Chr_position][0],
                                                         SNP[Chr_position][1],
                                                         SNP[Chr_position][2],
                                                         SNP[Chr_position][3]))
            Position_diff_new =[]
            for diff in Position_diff:
                Position_diff_new.append(str(diff))
            foutput = open(cov_file + '.frq.clean.snp', 'w')
            foutput.write('#A\tT\tG\tC\n')
            foutput.write(''.join(Output3))
            foutput.close()
            foutput = open(cov_file + '.frq.clean.snp.position.diff', 'w')
            foutput.write('\n'.join(Position_diff_new))
            foutput.close()
        except IOError:
            pass
    if args.strainfinder!= 'None':
        strain_finder(cov_file + '.frq.clean.snp')

def SNP_dynamics(output_files, SNP_outputfile, genomes):
    SNP_dynamics_pair = dict()
    SNP_dynamics_pair_output = []
    Total = len(output_files)
    for repeat_time in range(0,10):
        for i in range(1,Total):
            temp_list = set()
            while len(temp_list) < i:
                temp_list.add(random.choice(output_files))
            SNP, Position_diff = freq_call_sub(temp_list)
            SNP_dynamics_pair.setdefault(i,[])
            SNP_dynamics_pair[i].append(len(SNP))
    for i in SNP_dynamics_pair:
        for numbers in SNP_dynamics_pair[i]:
            SNP_dynamics_pair_output.append('%s\t%s\t%s\n' %(genomes,i,numbers))
    SNP_outputfile.write(''.join(SNP_dynamics_pair_output))


################################################## Programme ########################################################
if args.html == 'F':
    if args.cal == 'F':
        # process alignment and statistic
        f0=open('meta.decoder.sh','w')
        f0.write('#!/bin/bash\nsource ~/.bashrc\n')
        # number of bash scripts
        task = int(args.t)
        # generate map file
        genome_files=glob.glob(os.path.join(args.r,'*'+args.rf))
        f1=open(os.path.join(args.o,'ref.map.txt'),'w')
        for filename in genome_files:
            if '.ID.fasta' not in filename:
                for record in SeqIO.parse(filename, "fasta"):
                    f1.write(os.path.split(filename)[-1]+'\t'+str(record.id)+'\n')
        f1.close()
        cmd = 'cat '+os.path.join(args.r,'*'+args.rf)+'> %s\n' %(os.path.join(args.o,'all.ref.genomes.fasta'))
        #f0.write(cmd)

        # generate metagenome list file
        metagenome_files=glob.glob(os.path.join(args.i,'*'+args.inf))
        f1=open(os.path.join(args.o,'metagenome.list'),'w')
        cmd = ''
        if args.s != 1:
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
            #f0.write(cmd)
        # run strain finder preprocess
        cmd = 'python bin/0.run.py --fastqs %s  --ref %s  --map %s\n' % (
                    os.path.join(args.o, 'metagenome.list'),
                    os.path.join(args.o,'all.ref.genomes.fasta'),
                    os.path.join(args.o, 'ref.map.txt'))
        #f0.write(cmd)

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

        # run bowtie, PhaseFinder and SRID
        i=1
        print('processing all genomes %s' %(' '.join(genome_files)))
        print('processing all metagenomes %s' % (' '.join(metagenome_files)))
        for genomes in genome_files:
            if '.ID.fasta' not in genomes:
                for metagenomes in metagenome_files:
                    if '1'+args.inf in metagenomes or (args.s == 1 and '2'+args.inf not in metagenomes):
                        print('generating bowtie codes for mapping %s to %s' %(metagenomes, genomes))
                        fsub = open(str(int(i % task)) + '.sh', 'a')
                        fsub = open(args.os + '/' +str(int(i % task)) + '.sh', 'w')
                        fsub.write('#!/bin/bash\n')
                        cmd = ''
                        try:
                            ftest=open(os.path.join(args.o,os.path.split(metagenomes)[-1]+'_'+os.path.split(genomes)[-1]+'.out'),'r')
                        except IOError:
                            cmd += ''
                            cmd += 'python bin/PhaseFinder.py ratio -i %s -1 %s -2 %s -p 16 -o %s\n' % (
                                genomes+'.ID.fasta', metagenomes,
                                metagenomes.replace('1'+args.inf,'2'+args.inf),
                                os.path.join(args.o,os.path.split(metagenomes)[-1]+'_'+os.path.split(genomes)[-1]+'.out'))
                        try:
                            ftest=open(os.path.join(args.o,os.path.split(metagenomes)[-1]+
                                                    '_'+os.path.split(genomes)[-1]+'.sorted.bam.cov'),'r')
                        except IOError:
                            cmd += bowtie(genomes, metagenomes)
                            if args.s != 1:
                                metagenomes = metagenomes.replace('1' + args.inf, '2' + args.inf)
                                cmd += bowtie(genomes, metagenomes)
                        cmd += 'python bin/SRID.py -b %s -p 12 -r 100 -m 200 -s 71 -n 4 -o %s -t tmp\n' % (
                            os.path.join(args.o, os.path.split(metagenomes)[-1] + '_' + os.path.split(genomes)[-1] + '.sorted.bam'),
                            os.path.join(args.o, os.path.split(metagenomes)[-1] + '_' + os.path.split(genomes)[-1] + '.SRID.out.tab'))
                        cmd += '#ls -l %s\n#rm -rf %s\n' %(os.path.join(args.o, os.path.split(metagenomes)[-1] + '_' +
                               os.path.split(genomes)[-1] + '.SRID.bam'),
                                                         os.path.join(args.o,
                                                                      os.path.split(metagenomes)[-1] + '_' + os.path.split(genomes)[
                                                                          -1] + '.SRID.bam'))
                        i += 1
                        fsub.write(cmd)
                        print('output bowtie codes for mapping %s to %s' % (metagenomes, genomes))
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
                cmd += 'mv %s > %s \n' %(os.path.join(args.r,'/*'+args.rf+'.*'),args.o)
                i += 1
                #fsub.write(cmd)
                fsub.close()

        shfiles = glob.glob('%s/*.sh'%(args.os))
        for files in shfiles:
            if 'meta.decoder.sh' not in files:
                f0.write(('sh %s\n')%(files))
                #f0.write('jobmit %s %s.single\n' %(files,os.path.split(files)[1]))
        f0.close()
    else:
        genome_files = glob.glob(os.path.join(args.r, '*' + args.rf))
        Length = dict()
        try:
            for lines in open(os.path.join(args.r, 'Genome_size.txt')):
                Length.setdefault(lines.split('\t')[0], float(lines.split('\t')[1].replace('\r', '').replace('\n', '')))
        except IOError:
            pass
        # calculate average coverage
        all_output = open((os.path.join(args.o, 'all.bam.cov.sum')), 'w')
        all_output.write('metagenome\tgenome\tcoverage\tdepth\n')
        output_files = glob.glob(os.path.join(args.o, '*.bam.cov'))
        for output_filename in output_files:
            cal_cov(output_filename)
        all_output.close()
        # calculate average SNPs
        all_output = open((os.path.join(args.o, 'all.flt.vcf.snpden.sum')), 'w')
        all_output.write('metagenome\tgenome\tavgSNP_per_1kb\ttotal_SNPs\n')
        output_files = glob.glob(os.path.join(args.o, '*.flt.vcf.snpden'))
        for output_filename in output_files:
            snpsum(output_filename)
        all_output.close()
        # calculate allel frequency
        output_files = glob.glob(os.path.join(args.o, '*.flt.vcf'))
        for vcf_file in output_files:
            cov_file = vcf_file.replace('.flt.vcf', '.sorted.bam.cov')
            freq_call(vcf_file, cov_file)
        # calculate SNPs dynamics among samples
        foutput = open(os.path.join(SNP_outputdir, 'all.snp.dynamics.txt'), 'w')
        foutput.write('genome\tsample_number\ttotal_snps\n')
        for genomes in genome_files:
            genomes = os.path.split(genomes)[1]
            output_files = glob.glob(os.path.join(args.o, '*%s*.flt.vcf') % (genomes))
            SNP_dynamics(output_files, foutput, genomes)
        foutput.close()
        # calculate SNPs diversity
        foutput = open(os.path.join(SNP_outputdir, 'all.strain.diversity.txt'), 'w')
        foutput.write('genome\tsample\ttotal_strains\tabu_list\n')
        output_files = glob.glob(os.path.join(args.o, '*.abu'))
        foutput_list = []
        for abu_file in output_files:
            sample = os.path.split(abu_file)[1].split(args.inf + '_')[0]
            genomes = abu_file.split(args.inf + '_')[1].split(args.rf)[0]
            abu_list = []
            strain_num = 0
            for lines in open(abu_file, 'r'):
                lines = lines.replace('\r', '').replace('\n', '')
                abu_list.append(lines)
                strain_num+=len(lines.split('\t'))
            foutput_list.append('%s\t%s\t%s\t%s\n' % (genomes, sample,strain_num, ''.join(abu_list)))
        foutput.write(''.join(foutput_list))
        foutput.close()
else:
    # convert results to html files
    output_files=glob.glob(os.path.join(args.o,'*.flt.vcf.*'))
    for output_filename in output_files:
        if '.html' not in output_filename:
            tsv_to_html(output_filename)
    output_files = glob.glob(os.path.join(args.o, '*.flt.vcf'))
    cmds = ''
    for output_filename in output_files:
        #visualization
        cmds += visual(output_filename, os.path.join(args.o,'plots'))
    f0 = open('meta.decoder.visual.sh', 'w')
    f0.write('#!/bin/bash\n'+cmds)
    f0.close()