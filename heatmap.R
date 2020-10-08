#inputpath = 'example'
Allele1=read.delim(paste(inputpath,'all.snp.profile',sep = '/'),header=F)
Allele1=Allele1[!duplicated(paste(Allele1$V1,Allele1$V3)),]
Allele1 = Allele1[which(Allele1$V2!=''),]
library(stringr)
Allele1$allele=str_split_fixed(Allele1$V1, "--", 2)[,2]
Allele1$position=str_split_fixed(Allele1$allele, "_", 2)[,1]
Allele1$allele=str_split_fixed(Allele1$allele, "_", 2)[,2]
Allele1$genotypes = paste(str_split_fixed(Allele1$V3, "_", 2)[,1],
                          str_split_fixed(Allele1$V3, ":", 2)[,2],
                          sep = ':')
Allele1$refgenome = str_split_fixed(Allele1$V1, "_", 2)[,1]
Alleleall = Allele1
allgenome = unique(Alleleall$refgenome)
for (genome in allgenome)
{Allele1=Alleleall[which(Alleleall$refgenome == genome),]
  genotypes = unique(Allele1$genotypes)
position = unique(Allele1$position)
Allele = matrix('ref',ncol=length(genotypes),
                nrow=length(position))
colnames(Allele)=genotypes
row.names(Allele)=position
for(i in 1:nrow(Allele1))
{
  if(Allele1$V2[i]!=toString(Allele1$allele[i]))
  {n=which(colnames(Allele)==Allele1$genotypes[i])
  m=which(row.names(Allele)==Allele1$position[i])
  Allele[m,n]=Allele1$allele[i]}
}

library(pheatmap)
library(RColorBrewer)
mat= Allele
Color_set = c('#e5e5e5','#ca0020',
              '#f4a582','#92c5de','#0571b0')

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
         color = Color_set,
         show_rownames= TRUE,
         show_colnames = TRUE,
         cluster_cols = TRUE,
         legend_labels = c('Major allele','A','T','G','C'),
         filename = paste(inputpath,paste(genome,'snp.profile.SNP.pdf',sep = '.'),sep = '/'),
         cluster_rows = FALSE)}

