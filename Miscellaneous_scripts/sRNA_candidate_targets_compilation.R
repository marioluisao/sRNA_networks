#Compile list of candidate sRNA targets supported by experimental data
#Working directory should be the Miscellaneous_scripts folder - change to your own path
#Load the sRNA prior network used in original E. coli Inferelator run to extract manually selected sRNA priors
load("../Inferelator_output_files/Ecoli_8sRNAs/params_and_input.RData")
prior.network<-IN$priors.mat
#Save sRNAs names
sRNAs<-colnames(prior.network)[140:147]
#Create function to keep only the loci tag in the gene IDs (for example: ryhB_b4451_4)
translate.names<-function(geneNames)
{
  output<-c()
  for(i in geneNames)
  {
    temporal.name<-strsplit(as.character(i),"_")[[1]][2]
    output<-c(output,temporal.name)
  }
  output
}
#Create gold standard network (GS) table for sRNAs
sRNAs.gs<-c()
for(t in sRNAs)
{
  sRNA.priors<-names(which(prior.network[,t]!=0))
  sRNA.priors<-translate.names(sRNA.priors)
  sRNAs.gs<-rbind(sRNAs.gs,cbind(rep(t,length(sRNA.priors)),sRNA.priors))
}
#Read list of candidate targets compiled for each sRNA
sRNAs.candidate.targets.table<-c()
ryhB.targets<-read.table("../Miscellaneous_data/RyhB_putative_targets1.txt")
ryhB.targets<-as.character(ryhB.targets[,2])
sRNAs.candidate.targets.table<-rbind(sRNAs.candidate.targets.table,cbind(rep(sRNAs[1],length(ryhB.targets)),ryhB.targets))
spf.targets<-read.table("../Miscellaneous_data/Spot42_putative_targets1.txt")
spf.targets<-as.character(spf.targets[,2])
sRNAs.candidate.targets.table<-rbind(sRNAs.candidate.targets.table,cbind(rep(sRNAs[2],length(spf.targets)),spf.targets))
gcvB.targets<-read.table("../Miscellaneous_data/GcvB_putative_targets1.txt")
gcvB.targets<-as.character(gcvB.targets[,2])
sRNAs.candidate.targets.table<-rbind(sRNAs.candidate.targets.table,cbind(rep(sRNAs[3],length(gcvB.targets)),gcvB.targets))
micA.targets<-read.table("../Miscellaneous_data/MicA_putative_targets1.txt")
micA.targets<-as.character(micA.targets[,2])
sRNAs.candidate.targets.table<-rbind(sRNAs.candidate.targets.table,cbind(rep(sRNAs[4],length(micA.targets)),micA.targets))
omrA.targets<-read.table("../Miscellaneous_data/OmrA_putative_targets1.txt")
omrA.targets<-as.character(omrA.targets[,2])
sRNAs.candidate.targets.table<-rbind(sRNAs.candidate.targets.table,cbind(rep(sRNAs[5],length(omrA.targets)),omrA.targets))
cyaR.targets<-read.table("../Miscellaneous_data/CyaR_putative_targets1.txt")
cyaR.targets<-as.character(cyaR.targets[,2])
sRNAs.candidate.targets.table<-rbind(sRNAs.candidate.targets.table,cbind(rep(sRNAs[6],length(cyaR.targets)),cyaR.targets))
rybB.targets<-read.table("../Miscellaneous_data/RybB_putative_targets1.txt")
rybB.targets<-as.character(rybB.targets[,2])
sRNAs.candidate.targets.table<-rbind(sRNAs.candidate.targets.table,cbind(rep(sRNAs[7],length(rybB.targets)),rybB.targets))
fnrS.targets<-read.table("../Miscellaneous_data/FnrS_putative_targets1.txt")
fnrS.targets<-as.character(fnrS.targets[,2])
sRNAs.candidate.targets.table<-rbind(sRNAs.candidate.targets.table,cbind(rep(sRNAs[8],length(fnrS.targets)),fnrS.targets))
#Merge this table with the sRNAs GS network created above
sRNA.full.candidate.targets<-rbind(sRNAs.candidate.targets.table,sRNAs.gs)
#Remove duplicated sRNA-gene interactions
sRNA.full.candidate.targets<-unique.matrix(sRNA.full.candidate.targets)
colnames(sRNA.full.candidate.targets)<-c("sRNA","Target")
#Read Microbes online operon predictions for E. coli K12
operon.structure<-read.table("../Miscellaneous_data/microbes_online_Ecoli_K12_operons.txt",header=T)
operon.list<-list()
#Define the table rows with information about polycistronic operons
true.positions<-which(operon.structure$bOp==TRUE)
#Start the operon list - each element of the list would be an operon
operon.list[[1]]<-c(as.character(operon.structure$SysName1[true.positions[1]]),as.character(operon.structure$SysName2[true.positions[1]]))
for(s in 2:length(true.positions))
{
  
  current.position<-true.positions[s]
  previous.position<-true.positions[s-1]
  delta.position<-current.position - previous.position
  #Current gene pair (located in the same operon)
  current.line<-c(as.character(operon.structure$SysName1[current.position]),as.character(operon.structure$SysName2[current.position]))
  #Check genes in the last operon in the list
  previous.operon<-tail(operon.list,1)[[1]]
  #Add current gene pair to last entry in the operon list if there is overlap between them
  if(delta.position==1 & length(intersect(previous.operon,current.line))>0)
  {
    temporal.operon<-union(previous.operon,current.line)
    operon.list[[length(operon.list)]]<-temporal.operon
  }
  #Otherwise add current gene pair as a new element in the the operon list 
  if(delta.position!=1 | length(intersect(previous.operon,current.line))==0)
  {
    operon.list[[length(operon.list)+1]]<-current.line
  }
}
#Create function that expands list of candidate sRNA targets by including genes located in the same operons of experimentally supported targets compiled above
expand.targets<-function(input.table)
{
  output<-c()
  regulators<-unique(input.table[,1])
  for(g in regulators)
  {
    original.targets<-input.table[which(input.table[,1]==g),2]
    expanded.regulon<-c()
    for(d in original.targets)
    {
      initial.gene<-d
      expanded.operon<-initial.gene
      for(z in 1:length(operon.list))
      {
        if(length(intersect(operon.list[[z]],initial.gene))>0)
        {
          expanded.operon<-operon.list[[z]]
        }
      }
      expanded.regulon<-union(expanded.regulon,expanded.operon)
    }
    output<-rbind(output,cbind(rep(g,length(expanded.regulon)),expanded.regulon))
  }
  colnames(output)<-c("sRNA","Target")
  output
}
#Apply the function to the table of candidate sRNA-mRNA interactions compiled above
expanded.sRNA.gs<-expand.targets(sRNA.full.candidate.targets)