import glob
import os
import argparse

# Arguments and declarations
parser.add_argument("-i",
                    help="input dir of vcfs output by meta_decoder", type=str,
                    default='result_decoder',metavar='result_decoder')
args = parser.parse_args()

# function
def freq_call_sub(vcf_file):
    SNP = dict()
    SNP_count = dict()
    Total = 0
    SNP_major = dict()
    for lines in open(vcf_file, 'r'):
        if not lines.startswith("#"):
            lines_set = lines.split('\t')
            Depth = int(lines_set[7].split('DP=')[1].split(';')[0])
            print(vcf_file,Depth)
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
                            print('multiple')
                            Allels_frq[Allels[allels]] += int(Sum_alt * 0.5)
                        else:
                            pass
                for numbers in Allels_frq:
                    if numbers > 0:
                        temp += 1
                if temp >= 2:
                    SNP.setdefault(Chr_position, Allels_frq)
                    SNP_count[Chr] += 1
                    if Sum_ref < Sum_alt:
                        SNP_major.setdefault(Chr_position, [lines_set[4],
                                                            lines_set[3]])
                    else:
                        SNP_major.setdefault(Chr_position, [lines_set[3],
                                                            lines_set[4]])
    return [SNP, SNP_count, SNP_major, Total]

# set up
vcf = glob.glob(os.path.join(args.i,'*flt.vcf'))
snp_profile = dict()
snp_number = dict()
snp_profile_sample = dict()
Set=['major','minor']
Allels = dict()
Allels['A']=0
Allels['T']=1
Allels['G']=2
Allels['C']=3
Allels_order = ['A','T','G','C']

for files in vcf:
  SNP,SNP_count,SNP_major,Total = freq_call_sub(files)
  samplename = files.split('.fa.blast.txt.filter.aa.flt.vcf')[0]
  for Chr_i in SNP_major:
    temp = 0
    for major in SNP_major[Chr_i]:
        snp_profile.setdefault('%s_%s' % (Chr_i, major), [])
        snp_profile['%s_%s' % (Chr_i, major)].append(samplename + ':' + Set[temp])
        snp_profile_sample.setdefault(samplename, [])
        snp_profile_sample[samplename].append(['%s_%s' % (Chr_i, major), SNP[Chr_i][Allels[major]]])
        temp += 1
  for Chr in SNP_count:
    snp_number.setdefault(samplename,SNP_count[Chr])

f1=open(os.path.join(args.i,'all.snp.number'),'w')
for samplename in snp_number:
  f1.write('%s\t%s\n'%(samplename,snp_number[samplename]))

f1.close()

f1=open(os.path.join(args.i,'all.snp.profile'),'w')
for Chr_i in snp_profile:
  f1.write('%s\t%s\n'%(Chr_i,'\t'.join(snp_profile[Chr_i])))

f1.close()

f1=open(os.path.join(args.i,'all.snp.profile.sample'),'w')
for position in snp_profile_sample:
    for Chr_i_major in snp_profile_sample[position]:
            for samplename_allele in snp_profile[Chr_i_major[0]]:
                f1.write('%s\t%s\n'%(Chr_i_major[0],(samplename_allele)))

f1.close()