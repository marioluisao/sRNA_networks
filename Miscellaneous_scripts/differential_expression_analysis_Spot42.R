#Perform differential expression analysis of spf transcriptional data (Beisel and Storz 2011)
#Work directory should be the Miscellaneous_scripts folder - change accordingly
#Read table mapping E. coli genes and microarray probes 
genes.probes.table<-read.table("../Miscellaneous_data/Spot 42_transcriptional_data/ecoli_genechip2_genes_to_probes.txt",header=2)
#Verify that there is a single  probe per gene (limited to genes present in the original expression compendium used for inference)
#Load expression matrix
load("../Inferelator_output_files/Ecoli_8sRNAs/params_and_input.RData")
ecoli.genes<-rownames(IN$exp.mat)
probe.counts<-c()
for(i in ecoli.genes)
{
  current.gene<-strsplit(i,"_")[[1]][2]
  gene.positions<-grep(current.gene,genes.probes.table$gene)
  probe.counts<-c(probe.counts,length(gene.positions))
}  
#There are 3818 genes (out of 4297) with one probe in the array
#Read the expression data (downloaded from NCBI GEO - accesion number GSE24875)
spf.matrix<-read.csv("../Miscellaneous_data/Spot 42_transcriptional_data/Spot 42_microarray_normalized_raw_data_Beisel_and _Storz_2011_GSE24875.csv",row.names=1)
#Replace probes IDs with gene names
new.spf.expression.matrix<-c()
genes.present.in.matrix<-c()
for(i in ecoli.genes)
{
  current.gene<-strsplit(i,"_")[[1]][2]
  gene.position<-grep(current.gene,genes.probes.table[,1])
  #If the current gene was found
  if(length(gene.position)>0)
  {
    current.probe<-as.character(genes.probes.table[gene.position,2])
    current.probe.position<-which(rownames(spf.matrix)==current.probe)
    new.spf.expression.matrix<-rbind(new.spf.expression.matrix,spf.matrix[current.probe.position,])
    genes.present.in.matrix<-c(genes.present.in.matrix,current.gene)
  }
}
#Name rows in the new matrix
rownames(new.spf.expression.matrix)<-genes.present.in.matrix
#Remove genes that were absent (indicated with an "A") in any of the six experiments 
genes.presence<-sapply(1:nrow(new.spf.expression.matrix),function(x){length(which(new.spf.expression.matrix[x,seq(2,12,by=2)]=="P"))})
new.spf.expression.matrix2<-new.spf.expression.matrix[which(genes.presence == 6),]
#Delete the columns with the presence/absence information
new.spf.expression.matrix2<-new.spf.expression.matrix2[,-1*seq(from=2,to=12,by=2)]
#Perform differential expression analysis using bayesT (Baldi and Long, 2001)
#This step requires the cyberT R source code available at http://cybert.ics.uci.edu
source("../Rcode_Bayesian_Ttest/cyberTtest.R")
library("multtest")
#Run differential expression analysis 
diff.exp.analysis<-bayesT(new.spf.expression.matrix2,numC=3,numE=3,bayes=T,conf=7,ppde=T, doMulttest = T)
#Define diff. exp. genes (DEGs)
DEGs_spf<-rownames(diff.exp.analysis)[which(diff.exp.analysis$pVal <= 0.01)]
DEGs_spf<-DEGs_spf[order(DEGs_spf)]
write.csv(file="../Miscellaneous_data/Spot 42_transcriptional_data/spf_diff_exp_analysys_output.csv",diff.exp.analysis[DEGs_spf,],quote = F)