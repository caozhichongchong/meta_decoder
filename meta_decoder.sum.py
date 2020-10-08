import glob
import os
import argparse
import subprocess

# Arguments and declarations
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("-i",
                    help="input dir of vcfs output by meta_decoder", type=str,
                    default='result_decoder',metavar='result_decoder')
parser.add_argument("-fq",
                    help="file extension of fastq #1 files", type=str,
                    default='_1.fastq',metavar='_1.fastq')
parser.add_argument('-R',
                          help="Optional: complete path to R if not in PATH",
                          metavar="/usr/local/bin/R",
                          action='store', default='R', type=str)
args = parser.parse_args()

workingdir=os.path.abspath(os.path.dirname(__file__))
# function
def freq_call_sub(vcf_file,REF):
    SNP = dict()
    SNP_count = dict()
    Total = 0
    SNP_major = dict()
    for lines in open(vcf_file, 'r'):
        if not lines.startswith("#"):
            lines_set = lines.split('\t')
            Depth = int(lines_set[7].split('DP=')[1].split(';')[0])
            Chr = lines_set[0]
            SNP_count.setdefault(Chr, 0)
            Position = int(lines_set[1])
            Chr_position = '%s--%s' % (Chr, Position)
            Allels_frq = [0, 0, 0, 0]
            temp = 0
            if Total == 0:
                Total = len(lines_set) - 10
            if "INDEL" not in lines_set[7] and lines_set[6] != 'LowQual':
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
                for numbers in Allels_frq:
                    if numbers > 0:
                        temp += 1
                if temp >= 2:
                    SNP.setdefault(Chr_position, Allels_frq)
                    REF.setdefault(Chr_position, lines_set[3])
                    SNP_count[Chr] += 1
                    if Sum_ref < Sum_alt:
                        SNP_major.setdefault(Chr_position, [lines_set[4],
                                                            lines_set[3]])
                    else:
                        SNP_major.setdefault(Chr_position, [lines_set[3],
                                                            lines_set[4]])
    return [SNP, SNP_count, SNP_major, Total,REF]

# set up
vcf = glob.glob(os.path.join(args.i,'*flt.vcf'))
snp_profile = dict()
snp_number = dict()
REF = dict()
#snp_profile_sample = dict()
Set=['major','minor']
Allels = dict()
Allels['A']=0
Allels['T']=1
Allels['G']=2
Allels['C']=3
Allels_order = ['A','T','G','C']
print('run example: python meta_decoder.sum.py -i example/ -fq 1.fastq -R R')
for files in vcf:
    print('processing vcf %s'%(files))
    SNP, SNP_count, SNP_major, Total, REF = freq_call_sub(files, REF)
    samplename = os.path.split(files)[-1].split(args.fq)[0]
    i = 0
    for Chr_i in SNP_major:
        temp = 0
        for major in SNP_major[Chr_i]:
            snp_profile.setdefault('%s_%s' % (Chr_i, major), [])
            snp_profile['%s_%s' % (Chr_i, major)].append(samplename + ':' + Set[temp])
            # snp_profile_sample.setdefault(samplename, [])
            # snp_profile_sample[samplename].append(['%s_%s' % (Chr_i, major), SNP[Chr_i][Allels[major]]])
            temp += 1
        i += 1
    for Chr in SNP_count:
        snp_number.setdefault(samplename, SNP_count[Chr])

print('output all snp numbers %s'%(os.path.join(args.i,'all.snp.number')))
f1=open(os.path.join(args.i,'all.snp.number'),'w')
for samplename in snp_number:
  f1.write('%s\t%s\n'%(samplename,snp_number[samplename]))

f1.close()

print('output all snp profiles %s'%(os.path.join(args.i,'all.snp.profile')))
f1=open(os.path.join(args.i,'all.snp.profile'),'w')
for Chr_i in snp_profile:
  f1.write('%s\t%s\t%s\n'%(Chr_i,REF['_'.join(Chr_i.split('_')[:-1])],'\t'.join(snp_profile[Chr_i])))

f1.close()

# run heatmap.R
print('generating heatmap.R')
newRcode_temp = os.path.join(args.i,'heatmap_temp.R')
newRcode = os.path.join(args.i,'heatmap.R')
f1 = open(newRcode_temp,'w')
f1.write('inputpath = \'%s/\'\n'%(args.i))
f1.close()
os.system('cat %s %s > %s'%(newRcode_temp,os.path.join(workingdir,'heatmap.R'),newRcode))
print('plot all snp profiles %s/*.snp.profile.SNP.pdf'%(args.i))
os.system('%s --vanilla < %s'%(args.R, newRcode))

#f1=open(os.path.join(args.i,'all.snp.profile.sample'),'w')
#for position in snp_profile_sample:
#    for Chr_i_major in snp_profile_sample[position]:
#            for samplename_allele in snp_profile[Chr_i_major[0]]:
#                f1.write('%s\t%s\n'%(Chr_i_major[0],(samplename_allele)))

#f1.close()