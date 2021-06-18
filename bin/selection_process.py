################################################### END ########################################################
################################################### SET PATH ########################################################
import glob
import os
from Bio import SeqIO
from Bio.Seq import Seq
import statistics
import itertools
import random
# set up path
import argparse
############################################ Arguments and declarations ##############################################
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('optional arguments')
required.add_argument("-vcf",
                      help="path to flt vcf file",
                      type=str, default='sample.variants.flt.vcf',
                      metavar='sample.variants.flt.vcf')
required.add_argument("-ref",
                      help="path to reference genome",
                      type=str, default='ref.fasta',
                      metavar='ref.fasta')
# requirement for software calling
optional.add_argument('-pro',
                          help="Optional: complete path to prodigal if not in PATH",
                          metavar="/usr/local/bin/prodigal",
                          action='store', default='prodigal', type=str)
optional.add_argument('-fasttree',
                          help="Optional: complete path to fasttree if not in PATH",
                          metavar="/usr/local/bin/FastTreeMP",
                          action='store', default='FastTreeMP', type=str)

################################################## Definition ########################################################
args = parser.parse_args()
# set up parameters
vcf_name = '.flt.snp.vcf'
Qual = 29.5
min_minor_ALT = 3 # at least 3 reads mapped
outputname_set = ['filtered']
# Set up A T G C
Allels = dict()
Allels['A']=0
Allels['T']=1
Allels['G']=2
Allels['C']=3
Allels_order = ['A','T','G','C']
################################################### Function ########################################################
# set up functions
def ALT_freq(Allels_count,all_ALT):
    Major_ALT = []
    Minor_ALT = []
    ALT_set = dict()
    ALT_frq_set = set()
    for alleles in range(0, 4):
        ALT_frq = int(Allels_count[alleles])
        ALT_set.setdefault(ALT_frq, set())
        ALT_set[ALT_frq].add(alleles)
        ALT_frq_set.add(ALT_frq)
    ALT_frq_set = sorted(ALT_frq_set,reverse=True)
    for ALT_frq in ALT_frq_set:
        for alleles in ALT_set[ALT_frq]:
            if Major_ALT == []:
                Major_ALT = [Allels_order[alleles],ALT_frq]
            else:
                Minor_ALT.append([Allels_order[alleles],ALT_frq])
            if ALT_frq > 0:
                all_ALT.add(Allels_order[alleles])
    return [Major_ALT,Minor_ALT,all_ALT]

def vcf_to_txt(lines,output_list):
    lines_set = lines.split('\n')[0].split('\t')
    if len(lines_set) >9:
        CHR = lines_set[0]
        POS = int(lines_set[1])
        temp_line = []
        temp_line.append(CHR)
        temp_line.append(str(POS))
        i = 9
        for Subdepth_all in lines_set[9:]:
            Subdepth = Subdepth_all.split(':')[-1].replace('\n', '').split(',')
            total_sub_depth = sum(int(Subdepth_sub) for Subdepth_sub in Subdepth)
            temp_line.append(str(total_sub_depth))
            i += 1
        output_list.append('\t'.join(temp_line)+'\n')
    else:
        print(lines)

def curate_REF(allels_set,Depth4):
    Subdepth = Depth4.split(':')[-1].replace('\n', '').split(',')
    Subdepth_REF = int(Subdepth[0]) + int(Subdepth[1])
    Subdepth_ALT = int(Subdepth[2]) + int(Subdepth[3])
    if Subdepth_REF <= Subdepth_ALT:
        return [allels_set[1],1]
    else:
        return [allels_set[0],0]

def outputvcf(output_name):
    vcf_file_filtered = open(vcf_file + '.%s.snpfreq.txt' % (output_name), 'w')
    vcf_file_filtered.write('CHR\tPOS\tMajor_ALT\tMinor_ALT\tQuality\tGene\tGene_POS\tN_or_S\tAA_change\tallele_freq\n'+''.join(vcf_file_list_freq))
    vcf_file_filtered.close()

def SNP_seq(seq1, seq2, POS_info,POS_info_CHR,POS_info_CHR_LEN,POS_info_output,G1,G2):
    SNP_total = 0
    j = 0
    POS_DIS = []
    total_length = len(seq1)
    for i in range(0, total_length):
        if seq1[i] != seq2[i]:
            # a SNP
            SNP_total += 1
            CHR = POS_info_CHR[i]
            POS = POS_info[i]
            LEN = POS_info_CHR_LEN[CHR]
            if CHR == POS_info_CHR[j]:  # same CHR
                DIS = abs(POS - POS_info[j])
                POS_DIS.append(DIS)  # POS diff
                POS_info_output.append('%s\t%s\t%s\t%s\t%s\t%s\t\n' % (G1, G2, CHR, POS, DIS, LEN))
            else:  # new CHR
                POS_info_output.append('%s\t%s\t%s\t%s\t%s\t%s\t\n' % (G1, G2, CHR, POS, 0, LEN))
            j = i
    return SNP_total

def translate(seq):
    seq = Seq(seq)
    try:
        return seq.translate()
    except ValueError:
        try:
            return seq.translate(seq.complement())
        except ValueError:
            return ['None']

def dnORds(amino1, amino2):
    if amino1 == amino2:
        return 'S'
    else:
        return 'N'

def causeSNP(seq,position,ALT,Reverse_chr):
    if Reverse_chr == 1:
        ALT=str(Seq(ALT).reverse_complement())
    seq = list(seq)
    seq[position - 1]=ALT
    return ''.join(seq)

def loaddatabase(database):
    # load database seq
    Mapping = dict()
    Mapping_loci = dict()
    reference_database = os.path.split(database)[-1]
    print('reference database set as %s' % (reference_database))
    Ref_seq = dict()
    Reverse = []
    for record in SeqIO.parse(database, 'fasta'):
        record_id = str(record.id)
        record_seq = str(record.seq)
        Ref_seq.setdefault(record_id, record_seq)
        Mapping.setdefault(record_id, len(record_seq))
        description = str(record.description).replace(' ', '').split('#')
        contig = '_'.join(record_id.split('_')[0:-1])
        Mapping_loci.setdefault(contig, [])
        if float(description[3]) == -1.0: # reverse str
            Reverse.append(record_id)
        Mapping_loci[contig].append([float(description[1]),
                                     float(description[2]),
                                     record_id])
    return [Ref_seq,Mapping,Mapping_loci,Reverse]

def contig_to_gene(CHR, POS):
    all_genes = Mapping_loci.get(CHR,[])
    Reverse_chr = 0
    for a_gene in all_genes:
        POS1, POS2, GENE = a_gene
        if POS >= POS1 and POS <= POS2:
            Ref_seq_chr = Ref_seq.get(GENE, 'None')
            Gene_length = len(Ref_seq_chr)
            if GENE in Reverse:  # reversed
                POS_gene = Gene_length-(int(POS-POS1))
                Reverse_chr = 1
            else:
                POS_gene = int(POS-POS1)+1
            codon_start = POS_gene - 1 - int((POS_gene - 1) % 3)
            return [GENE,POS_gene,codon_start,Ref_seq_chr,Reverse_chr]
    return []

def SNP_check_all_fq(lines_set):
    temp_snp_line = []
    temp_snp_line_frq2 = []
    temp_snp_line_NS = ['None', 'None', 'None']
    temp_snp_line_AA = ''
    REF = lines_set[3]
    allels_set = [REF]
    lines_set_sub = lines_set[9:]
    CHR = lines_set[0]
    POS = int(lines_set[1])
    SNP_quality = lines_set[5]
    if '.' not in lines_set[4]:
        allels_set += lines_set[4].split(',')
    Total_alleles = len(allels_set)
    genome_order = 0
    Depth4 = lines_set[7].split('DP4=')[1].split(';')[0]
    REF, REF_where = curate_REF(allels_set, Depth4)  # as the major alt in the population
    all_ALT = set()
    for Subdepth_all in lines_set_sub:
        genome_order += 1
        Allels_frq = [0, 0, 0, 0]
        Subdepth = Subdepth_all.split(':')[-1].replace('\n', '').split(',')
        #total_sub_depth = sum(int(Subdepth_sub) for Subdepth_sub in Subdepth)
        for num_allels in range(0, Total_alleles):
            allels = allels_set[num_allels]
            Subdepth_alleles = int(Subdepth[num_allels])
            if Subdepth_alleles < min_minor_ALT:
                # remove small reads
                Subdepth_alleles = 0
            if allels in Allels:
                Allels_frq[Allels[allels]] += Subdepth_alleles
            else:
                pass
        # find major alt and calculate frequency
        Major_ALT, Minor_ALT,all_ALT = ALT_freq(Allels_frq,all_ALT)
        temp_snp_line_frq2.append(';'.join(str(frq_sub) for frq_sub in Allels_frq))
    if len(all_ALT) > 1: # with SNPs
        allels_set = list(all_ALT)
        # calculate NS
        gene_info = contig_to_gene(CHR, POS)
        if gene_info != []:
            Chr_gene, POS_gene, codon_start, Ref_seq_chr, Reverse_chr = gene_info
            if Ref_seq_chr != 'None':
                #  observed NS ratio calculated
                HS_gene.setdefault(Chr_gene,[])
                temp_snp_line_NS = [Chr_gene, str(POS_gene), '']
                if codon_start <= POS_gene - 1:
                    Ref_seq_chr = causeSNP(Ref_seq_chr, POS_gene, REF, Reverse_chr)
                    Ref_seq_codon = Ref_seq_chr[codon_start:(codon_start + 3)]
                    SNP_seq_chr = Ref_seq_chr
                    if len(Ref_seq_codon) == 3:
                        Ref_seq_aa = translate(Ref_seq_codon)[0]
                        temp_snp_line_AA += Ref_seq_aa
                        ALT_set = allels_set
                        ALT_set.remove(REF)
                        for ALT in ALT_set:
                            SNP_seq_chr = causeSNP(SNP_seq_chr, POS_gene, ALT, Reverse_chr)
                            SNP_seq_codon = SNP_seq_chr[codon_start:(codon_start + 3)]
                            SNP_seq_aa = translate(SNP_seq_codon)[0]
                            temp_snp_line_AA += SNP_seq_aa
                            temp_NorS = dnORds(Ref_seq_aa, SNP_seq_aa)
                            temp_snp_line_NS[-1] += temp_NorS
                            if 'N' in temp_NorS:
                                HS_gene[Chr_gene].append('N')
                            else:
                                HS_gene[Chr_gene].append('S')
        temp_snp_line.append(CHR)
        temp_snp_line.append(str(POS))
        temp_snp_line.append(REF)
        temp_snp_line.append(','.join([ALT for ALT in allels_set if ALT != REF]))
        vcf_file_list_freq.append(
            '\t'.join(temp_snp_line) + '\t%s\t%s\t%s\t%s\n' % (
                SNP_quality, '\t'.join(temp_snp_line_NS),
                temp_snp_line_AA, '\t'.join(temp_snp_line_frq2)))

def HS_filtering():
    HS_list = []
    for gene in HS_gene:
        num_N = HS_gene[gene].count('N')
        num_S = HS_gene[gene].count('S')
        if num_N + num_S > 1 and (num_S == 0 or float(num_N)/float(num_S) > 1):
            HS_list.append(gene)
    return HS_list

def outputvcf_HS(output_name,HS_list):
    vcf_file_list_freq_HS = []
    for lines in open(vcf_file + '.%s.snpfreq.txt' % (output_name), 'r'):
        if lines.split('\t')[5] in HS_list:
            vcf_file_list_freq_HS.append('True\t%s' %lines)
        else:
            vcf_file_list_freq_HS.append('False\t%s' % lines)
    vcf_file_filtered = open(vcf_file + '.%s.snpfreq.HS.txt' % (output_name), 'w')
    vcf_file_filtered.write(
        'High_select\tCHR\tPOS\tMajor_ALT\tMinor_ALT\tQuality\tGene\tGene_POS\tN_or_S\tAA_change\tallele_freq\n' + ''.join(
            vcf_file_list_freq_HS))
    vcf_file_filtered.close()

################################################### Main ########################################################
# run vcf filtering
vcf_file = args.vcf
output_name = outputname_set[0]
filesize = 0
try:
    filesize = int(os.path.getsize(vcf_file + '.%s.snpfreq.txt' % (output_name)))
except FileNotFoundError:
    pass
if filesize == 0:
    Total = 0
    donor_species = os.path.split(vcf_file)[-1].split(vcf_name)[0]
    vcf_file_list_freq = []
    HS_gene = dict()
    Ref_seq = dict()
    Mapping = dict()
    Mapping_loci = dict()
    database_file = args.ref
    print('running %s' % (donor_species))
    if database_file.split('.')[-1] != '.fna':
        # not gene file
        try:
            f1 = open(database_file + '.fna', 'r')
        except FileNotFoundError:
            os.system('%s -q -i %s -d %s.fna' % (args.pro, database_file, database_file))
        database_file = database_file + '.fna'
    Ref_seq, Mapping, Mapping_loci, Reverse = loaddatabase(database_file)
    for lines in open(vcf_file, 'r'):
        if not lines.startswith("#"):
            lines_set = lines.split('\n')[0].split('\t')
            CHR = lines_set[0]
            POS = int(lines_set[1])
            # a SNP confirmed in WGS mapping
            CHR_POS = '%s__%s' % (CHR, POS)
            Depth = int(lines_set[7].split('DP=')[1].split(';')[0])
            if Total == 0:
                Total = len(lines_set) - 9
            if "INDEL" not in lines_set[7] \
                    and (lines_set[6] != 'LowQual' or float(lines_set[5]) >= Qual):
                SNP_check_all_fq(lines_set)
    outputvcf(output_name)
    HS_list = HS_filtering()
    outputvcf_HS(output_name,HS_list)
    print('finished processing %s' % (donor_species))

