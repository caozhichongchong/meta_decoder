Allele1=read.delim('all.snp.profile',header=F)
Ref='GACCAATTGACAACAACAAAGACTCATCCGAGCTTTAAAGCCGGTGATAATATTACCGTTAACTACAAGATCGTAGAGGGTGGTAAAGAGCGTATCCAGAGTTTCCGTGGAGATGTGATCAAAACACAGGGAAATGGTGCTACTGCATCTTTCACTGTACGTAAGATCTCTGATGGAGTAGGTGTAGAACGTTTGTTCCCTATCCTTTCTCCAAATATCGACTCTATCGTATTTAACAAAGCCGGTAAATTTCGTCGCGCTAAATTATACTACCTGCGTGACCGCAGCGGTAAGAGCGCACGTATCAAGGAAAAA'
Allele1=Allele1[!duplicated(paste(Allele1$V1,Allele1$V2)),]
library(stringr)
Allele1$allele=str_split_fixed(Allele1$V1, "--", 2)[,2]
Allele1$position=str_split_fixed(Allele1$allele, "_", 2)[,1]
Allele1$allele=str_split_fixed(Allele1$allele, "_", 2)[,2]
Allele1$genotypes = paste(str_split_fixed(Allele1$V2, "_", 2)[,1],
                          str_split_fixed(Allele1$V2, ":", 2)[,2],
                          sep = ':')
genotypes = unique(Allele1$genotypes)
position = unique(Allele1$position)
Allele = matrix('ref',ncol=length(genotypes),
                nrow=length(position))
colnames(Allele)=genotypes
row.names(Allele)=position
for(i in 1:nrow(Allele1))
{
  n=which(colnames(Allele)==Allele1$genotypes[i])
  m=which(row.names(Allele)==Allele1$position[i])
  Allele[m,n]=Allele1$allele[i]
}
nchar(Ref)
temp = 234
substr(Ref, temp, temp)
Allele_correct=Allele
for(i in 1:nrow(Allele))
{for (j in 1:ncol(Allele))
{
  if(Allele_correct[i,j]=='ref')
    Allele_correct[i,j]=substr(Ref, as.numeric(row.names(Allele_correct)[i]), 
                               as.numeric(row.names(Allele_correct)[i]))
}}

library(pheatmap)
library(RColorBrewer)
mat= Allele
mat[which(mat=='ref')]=0
mat[which(mat=='A')]=1
mat[which(mat=='T')]=2
mat[which(mat=='G')]=3
mat[which(mat=='C')]=4
mat=matrix(as.numeric(mat),nrow=nrow(mat))
row.names(mat)=row.names(Allele)
colnames(mat)=colnames(Allele)
mat=mat[order(as.numeric(row.names(mat))),]

pheatmap(mat= mat,
         color = colorRampPalette((brewer.pal(n = 8, name = "Purples")))(5),
         show_rownames= TRUE,
         show_colnames = TRUE,
         cluster_cols = TRUE,
         cluster_rows = FALSE)