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
                if 'GW715' in vcf_file:
                    if Sum_ref < Sum_alt:
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
                else:
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


import glob
vcf = glob.glob('*.vcf')
snp_profile = dict()
snp_number = dict()
snp_profile_sample = dict()
Set=['major','minor']
for files in vcf:
  SNP,SNP_count,SNP_major,Total = freq_call_sub(files)
  samplename = files.replace('.inter1.fastq_it3_600_bin-contigs.fa.blast.txt.filter.aa.flt.vcf','')
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

f1=open('all.snp.number','w')
for samplename in snp_number:
  f1.write('%s\t%s\n'%(samplename,snp_number[samplename]))

f1.close()

f1=open('all.snp.profile','w')
for Chr_i in snp_profile:
  f1.write('%s\t%s\n'%(Chr_i,'\t'.join(snp_profile[Chr_i])))

f1.close()

Genotype = set()
f1=open('GW715.all.snp.profile.sample','w')
for Chr_i_major in snp_profile_sample['GW715-65-2-12-13_45461_9_1']:
        for samplename_allele in snp_profile[Chr_i_major[0]]:
                Genotype.add(samplename_allele)
                f1.write('%s\t%s\n'%(Chr_i_major[0],(samplename_allele)))

f1.close()
Genotype.add('GW715-65-2-12-13_45461_9_1:major')

f1=open('Other.all.snp.profile.sample','w')
for position in snp_profile_sample:
    for Chr_i_major in snp_profile_sample[position]:
            for samplename_allele in snp_profile[Chr_i_major[0]]:
                if samplename_allele in Genotype:
                    f1.write('%s\t%s\n'%(Chr_i_major[0],(samplename_allele)))

f1.close()