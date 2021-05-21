inputpath = '/Users/caozhichongchong/Desktop/AnniZhang/project/github/meta_decoder/example/'
#inputpath = 'example'
Allele1=read.delim(paste(inputpath,'all.snp.profile',sep = '/'),header=F)
Allele1=Allele1[!duplicated(paste(Allele1$V1,Allele1$V3)),]
Allele1 = Allele1[which(Allele1$V2!=''),]
library(stringr)
Allele1$allele=str_split_fixed(Allele1$V1, "--", 2)[,2]
Allele1$position=str_split_fixed(Allele1$allele, "_", 2)[,1]
Allele1$allele=str_split_fixed(Allele1$allele, "_", 2)[,2]
Allele1$genotypes = paste(str_split_fixed(str_split_fixed(Allele1$V3, "_", 2)[,1],":",2)[1],
                          str_split_fixed(Allele1$V3, ":", 2)[,2],
                          sep = ':')
Allele1$refgenome = str_split_fixed(Allele1$V1, "__", 2)[,1]
Allele1$contig = str_split_fixed(Allele1$V1, "__", 2)[,2]
Allele1$contig = str_split_fixed(Allele1$contig, "--", 2)[,1]

Alleleall = Allele1
allgenome = unique(Alleleall$refgenome)
for (genome in allgenome)
{
  Allele1=Alleleall[which(Alleleall$refgenome == genome),]
  genotypes = unique(Allele1$genotypes)
  contig = unique(Allele1$contig)
  position = unique(paste(Allele1$contig,Allele1$position))
  Allele = matrix('ref',ncol=length(genotypes)+2,
                  nrow=length(position))
  colnames(Allele)=c('POS','REF',genotypes)
  row.names(Allele)=position
  Allele[,1]=str_split_fixed(position, " ", 2)[,2]
  for(i in 1:nrow(Allele1))
  {
    if(Allele1$V2[i]!=toString(Allele1$allele[i]))
    {n=which(colnames(Allele)==Allele1$genotypes[i])
    m=which(row.names(Allele)==paste(Allele1$contig[i],Allele1$position[i]))
    Allele[m,n]=Allele1$allele[i]}
    else
    {
      m=which(row.names(Allele)==paste(Allele1$contig[i],Allele1$position[i]))
      Allele[m,2]=Allele1$allele[i]
    }
  }
  row.names(Allele)=str_split_fixed(position, " ", 2)[,1]
  allcontigs = unique(Allele1$contig)
  empty_Allele=matrix('empty',nrow = length(allcontigs),ncol=ncol(Allele))
  row.names(empty_Allele)=allcontigs
  empty_Allele[,1]=0
  Allele=rbind(Allele,
               empty_Allele)
  Allele=Allele[order(as.numeric(Allele[,1])),]
  Allele=Allele[order(row.names(Allele)),]
  row.names(Allele)[which(Allele[,1]!='0')]=Allele[which(Allele[,1]!='0'),1]
  library(pheatmap)
  library(RColorBrewer)
  mat= Allele[,-1]
  Color_set = c('#ffffff','#e5e5e5','#ca0020',
                '#f4a582','#92c5de','#0571b0')
  
  mat[which(mat=='empty')]=-1
  mat[which(mat=='ref')]=0
  mat[which(mat=='A')]=1
  mat[which(mat=='G')]=2
  mat[which(mat=='T')]=3
  mat[which(mat=='C')]=4
  mat=matrix(as.numeric(mat),nrow=nrow(mat))
  row.names(mat)=row.names(Allele)
  colnames(mat)=colnames(Allele)[-1]
  pheatmap(mat= mat,
           color = Color_set,
           show_rownames= TRUE,
           show_colnames = TRUE,
           cluster_cols = FALSE,
           legend_breaks =   c(-1,0,1,2,3,4),
           fontsize_col = 10,
           fontsize_row = 8,
           legend_labels = c('','Major allele','A','G','T','C'),
           cluster_rows = FALSE,
           width = 5+2*ncol(mat),
           height = min(5+0.1*nrow(mat),50),
           filename = paste(inputpath,paste(genome,'snp.profile.SNP.pdf',sep = '.'),sep = '/')
  )
}
