#This script was used to create the shuffled and noisy E. coli sRNA priors
#Working directory should be the Miscellaneous_scripts folder - change to your own path
#Create shuffled sRNA priors
#Load Inferelator output for standard run with manually selected sRNA priors
load("../Inferelator_output_files/Ecoli_8sRNAs/betas_frac_tp_100_perm_1--frac_fp_0_perm_1_1.1.RData")
load("../Inferelator_output_files/Ecoli_8sRNAs/combinedconf_frac_tp_100_perm_1--frac_fp_0_perm_1_1.1.RData")
load("../Inferelator_output_files/Ecoli_8sRNAs/params_and_input.RData")
#Save sRNAs names
sRNAs<-colnames(IN$priors.mat)[140:147]
#Loop to create ten random versions of shuffled sRNA priors
for(j in 1:10)
{
  #Compile the names of all sRNA priors in the original run
  sRNA.target.genes<-c()
  for(k in sRNAs)
  {
    sRNA.target.genes<-c(sRNA.target.genes,names(which(IN$priors.mat[,k]!=0)))
  }
  #Create matrix with TFs priors. Leave all sRNA columns as zeroes
  shuffled.sRNA.priors.matrix<-cbind(IN$priors.mat[,1:139],matrix(0,ncol=length(sRNAs),nrow=nrow(IN$priors.mat)))
  colnames(shuffled.sRNA.priors.matrix)[140:147]<-sRNAs
  #Randomly select the priors for each sRNA
  for(i in sRNAs)
  {
    #Randomly select the sRNA priors (the same number of priors originally selected)
    current.sRNA.random.sample<-sample(1:length(sRNA.target.genes),length(which(IN$priors.mat[,i]!=0)))
    random.targets<-sRNA.target.genes[current.sRNA.random.sample]
    #Remove the selected genes from the set of candidate sRNA targets
    sRNA.target.genes<-sRNA.target.genes[-1*current.sRNA.random.sample]
    #Replace entries for the selected genes in the corresponding column in the priors matrix
    shuffled.sRNA.priors.matrix[random.targets,i]<- -1
  }
  #Due to the existence of eight genes controlled by more than one sRNA, the resulting sRNA prior network may contain a slightly smaller number of total priors
  #Among the ten shuffled sRNA networks used for Fig 3A, five sets of sRNA priors contained the same number of original sRNA priors (75 interactions)
  #The other five contained between one and three less priors
  #Create the .csv file required to run the Inferelator
  write.table(shuffled.sRNA.priors.matrix,sep="\t",col.names=NA,quote=F,
              file=paste("../Input_files/inf_shuffle_sRNA_priors_",j,".csv",sep=""))
}

###Create noisy priors
source("sRNA_candidate_targets_compilation.R")
#Specify the ratio of true:false sRNA priors
noise<-c(1,2,5)
#Loop for creating noisy sRNA priors networks (merged with curated transcriptional prior network) for each noise level
for(k in noise)
{
  #Loop to create ten noisy prior networks
  for(u in 1:10)
  {
    #Create matrix with TFs priors. Leave all sRNA columns as zeroes
    noisy.sRNA.priors.matrix<-cbind(IN$priors.mat[,1:139],matrix(0,ncol=length(sRNAs),nrow=nrow(IN$priors.mat)))
    colnames(noisy.sRNA.priors.matrix)[140:147]<-sRNAs
    #Randomly select false positives for each sRNA
    for(i in sRNAs)
    {
      true.sRNA.targets<-names(which(IN$priors.mat[,i]!=0))
      #Set of candidate sRNA targets
      supported.targets<-expanded.sRNA.gs[which(expanded.sRNA.gs[,1]==i),2]
      supported.targets<-unlist(sapply(1:length(supported.targets),function(x){rownames(noisy.sRNA.priors.matrix)[grep(supported.targets[x],rownames(noisy.sRNA.priors.matrix))]}))
      #Define set of not supported targets
      potential.false.targets<-setdiff(rownames(noisy.sRNA.priors.matrix),supported.targets)
      #Randomly select false sRNA priors
      false.positives<-potential.false.targets[sample(1:length(potential.false.targets),length(true.sRNA.targets)*k)]
      #Combine true targets and false targets
      noisy.priors<-union(true.sRNA.targets,false.positives)
      noisy.sRNA.priors.matrix[noisy.priors,i]<- -1
    }
    write.table(noisy.sRNA.priors.matrix,sep="\t",col.names=NA,quote=F,
                file=paste("../Input_files/inf_noisy_sRNA_priors_1_n",k,"_r",u,".csv",sep=""))
  }
}