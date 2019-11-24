#This script reads the output of the Inferelator and creates a network model using a confidence score threshold defined by the user
#The confidence threshold is defined in terms of precision (a 0.5 precision cutoff was used in this work)
#This script assumes that the Inferelator output files has been already loaded
#The Inferelator groups regulators with similar profiles
#Save the names of the Inferelator-defined groups
columns.to.be.deleted<-grep("pred.group",colnames(comb.confs))
#Crete a matrix with the average regression coefficients (betas) of each TF-gene interaction
betas.mean<-matrix(nrow=dim(comb.confs)[1],ncol=dim(comb.confs)[2],0)
for(y in 1:length(betas))
{
  betas.mean<- betas.mean + betas[[y]]
}
combined.betas<- betas.mean/length(betas)
#Remove the groups
if(length(columns.to.be.deleted)>0)
{
  combined.betas<-combined.betas[,-1*columns.to.be.deleted]
}
#Create function to normalize confidence score matrices in the [0,1] interval
normalization<-function(confidence.matrix)
{
  score.range<-range(confidence.matrix)
  deltaScore<-score.range[2]-score.range[1]
  temporal.confidence.matrix<-(confidence.matrix-score.range[1])/deltaScore
  output<-temporal.confidence.matrix
  output
}
#Normalize the original confidence matrix
normalized.confidence.matrix<-normalization(as.matrix(comb.confs))
#Remove the groups
if(length(columns.to.be.deleted)>0)
{
  normalized.confidence.matrix<-normalized.confidence.matrix[,-1*columns.to.be.deleted]
}
#Create new gold standard (prior) matrix without the removed regulators
curated.gs.matrix<-matrix(0,nrow=dim(normalized.confidence.matrix)[1],ncol=dim(normalized.confidence.matrix)[2]
                          ,dimnames=list(rownames(normalized.confidence.matrix),colnames(normalized.confidence.matrix)))
original.gs.matrix<-IN$priors.mat
for(i in colnames(normalized.confidence.matrix))
{
  regulator.position.original.gs.matrix<-which(colnames(original.gs.matrix)==i)
  curated.gs.matrix[,i]<-original.gs.matrix[,regulator.position.original.gs.matrix]
}
#Functions to create precision-recall curves (from the Inferelator folder)
areaUnderCurve <- function(x, y)
{
  dx <- diff(x)
  my <- y[1:(length(y) - 1)] + diff(y) / 2
  return(sum(dx * my))
}
calcAupr <- function(pred, gs) {
  ord.idx <- order(abs(pred), decreasing = T)
  
  prec <- cumsum(gs[ord.idx]) / cumsum(rep(1, length(ord.idx))) #also known as positive predictive value
  rec  <- cumsum(gs[ord.idx]) / sum(gs)                     #also know as true positive rate
  fpr  <- cumsum(gs[ord.idx] == 0) / (length(gs) - sum(gs)) #false positive rate
  
  prec <- c(prec[1], prec)
  rec <- c(0, rec)
  fpr <- c(0, fpr)
  
  aupr <- areaUnderCurve(rec, prec)
  auroc <- areaUnderCurve(fpr, rec)
  
  return(list(prec=prec, rec=rec, fpr = fpr, AUPR = aupr, AUROC = auroc))
}
#Compute the Area Under the Precision Recall (AUPR) curve
precision.recall.metrics<-calcAupr(normalized.confidence.matrix,abs(curated.gs.matrix))
#Print aupr value
print(precision.recall.metrics$AUPR)
#Create final network using a 0.5 precision threshold
recovered.interactions.matrix<-curated.gs.matrix
novel.interactions.matrix<-matrix(0, nrow=nrow(normalized.confidence.matrix),ncol=ncol(normalized.confidence.matrix)
                                  ,dimnames=list(rownames(normalized.confidence.matrix),colnames(normalized.confidence.matrix)))
#Define the confidence cutoff
confidence.score.threshold<-normalized.confidence.matrix[order(normalized.confidence.matrix,decreasing=T)][max(which(precision.recall.metrics$prec==0.5))-1]
#Replace entries of TF-gene and sRNA-gene interactions in the GS network below the confidence threshold with zeroes
recovered.interactions.matrix[which(normalized.confidence.matrix < confidence.score.threshold)]<-0
#Print the number of recovered targets (from the GS network) in the inferred network
colSums(abs(recovered.interactions.matrix))
#Identify novel predictions and replace corresponding entries with the sign of the regression coefficients
novel.interactions.matrix[which(normalized.confidence.matrix >= confidence.score.threshold)]<- sign(combined.betas[which(normalized.confidence.matrix >= confidence.score.threshold)])
#Remove non-zero values that correspond to recovered priors
novel.interactions.matrix[which(recovered.interactions.matrix!=0)]<-0
#Print the number of novel targets per regulator
colSums(abs(novel.interactions.matrix))