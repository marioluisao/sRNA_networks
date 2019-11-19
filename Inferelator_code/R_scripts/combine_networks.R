library('Matrix')

get.betas <- function(beta.file) {
  load(beta.file)
  rows <- sort(rownames(betas.resc[[1]]))
  cols <- sort(colnames(betas.resc[[1]]))
  betas.resc <- lapply(betas.resc, function(br) br[rows, cols])
  betas <- lapply(betas, function(b) b[rows, cols])
  return(list(betas=betas, betas.resc=betas.resc))
}

get.mean.and.lh <- function(mat) {
  ret <- list()
  ret$mean <- apply(mat, 1, mean)
  tmp <- apply(mat, 1, quantile, probs=c(0.05, 0.95))
  ret$low <- tmp[1, ]
  ret$high <- tmp[2, ]
  return(ret)
}

to.confidence.score <- function(mat) {
  tmp <- rank(mat)
  return(1 - (tmp - 1) / (max(tmp) - 1))
}

eu.beta.file <- 'output/bsubtilis_eu_TFA_TRUE_BBSR_1.1_100_0/betas_frac_tp_100_perm_1--frac_fp_0_perm_1_1.1.RData'
us.beta.file <- 'output/bsubtilis_us_TFA_TRUE_BBSR_1.1_100_0/betas_frac_tp_100_perm_1--frac_fp_0_perm_1_1.1.RData'
gs.file <- 'input/eu/gs.mat.RData'
output.file <- 'output/combined_network.tsv'

eu.betas <- get.betas(eu.beta.file)
us.betas <- get.betas(us.beta.file)

# concatenate the lists for the two data sets
betas.resc <- c(eu.betas$betas.resc, us.betas$betas.resc)
betas <- c(eu.betas$betas, us.betas$betas)

beta.sign <- as.matrix(betas[[1]] * 0)
beta.non.zero <- as.matrix(betas[[1]] * 0)
for (n in 1:length(betas)) {
  beta.sign <- beta.sign + sign(betas[[n]])
  beta.non.zero <- beta.non.zero + (betas[[n]] != 0)
}

# we only care about interactions that are present in at least one bootstrap
non.zero <- which(beta.non.zero != 0, arr.ind=TRUE)
vec.ind <- which(beta.non.zero != 0)

# get mean beta score
betas.bs <- lapply(betas, function(x) x[vec.ind])
betas.bs <- matrix(unlist(betas.bs), length(betas.bs[[1]]))
#betas.mean <- apply(betas.bs, 1, mean)
betas.stats <- get.mean.and.lh(betas.bs)
rm(betas.bs)

# get mean rescaled beta
resc.betas <- lapply(betas.resc, function(x) x[vec.ind])
resc.betas <- matrix(unlist(resc.betas), length(resc.betas[[1]]))
resc.betas.stats <- get.mean.and.lh(resc.betas)
#resc.betas.mean <- apply(resc.betas, 1, mean)

# get mean of ranks
resc.betas.ranks <- apply(-resc.betas, 2, rank)
resc.betas.ranks.stats <- get.mean.and.lh(resc.betas.ranks)
rm(resc.betas, resc.betas.ranks)

# get confidence scores
conf.scores <- lapply(betas.resc, function(x) to.confidence.score(-as.matrix(x))[vec.ind])
conf.scores <- matrix(unlist(conf.scores), length(conf.scores[[1]]))
conf.scores.stats <- get.mean.and.lh(conf.scores)
rm(conf.scores)


# load the gs matrix (which was also the prior in this case)
gs.from.file <- as.matrix(local(get(load(gs.file))))
gs <- as.matrix(betas[[1]] * 0)
p.genes <- intersect(rownames(gs), rownames(gs.from.file))
p.tfs <- intersect(colnames(gs), colnames(gs.from.file))
#cat('p.genes', length(p.genes), 'p.tfs', length(p.tfs), '\n')
gs[p.genes, p.tfs] <- gs.from.file[p.genes, p.tfs]


# create data frame of non-zero scores 
net.df <- data.frame(regulator=colnames(beta.sign)[non.zero[, 2]], 
                     target=rownames(beta.sign)[non.zero[, 1]], 
                     beta.sign.sum=beta.sign[vec.ind],
                     beta.non.zero=beta.non.zero[vec.ind] / length(betas),
                     beta.low=betas.stats$low,
                     beta.mean=betas.stats$mean,
                     beta.high=betas.stats$high,
                     var.exp.low=resc.betas.stats$low,
                     var.exp.mean=resc.betas.stats$mean,
                     var.exp.high=resc.betas.stats$high,
                     var.exp.rank.low=resc.betas.ranks.stats$low,
                     var.exp.rank.mean=resc.betas.ranks.stats$mean,
                     var.exp.rank.high=resc.betas.ranks.stats$high,
                     conf.score.low=conf.scores.stats$low,
                     conf.score.mean=conf.scores.stats$mean,
                     conf.score.high=conf.scores.stats$high,
                     prior=gs[vec.ind])

iao <- order(-net.df$conf.score.mean, net.df$var.exp.rank.mean, -net.df$var.exp.mean, -abs(net.df$beta.non.zero), net.df$regulator, net.df$target)
#iao <- order(net.df$var.exp.rank.mean, -net.df$var.exp.mean, -abs(net.df$beta.non.zero), net.df$regulator, net.df$target)
net.df <- net.df[iao, ]

# get size of network at precision 0.5
prec <- cumsum((net.df$prior != 0)) / (1:length(net.df$prior))
net.size <- min(which(prec < 0.5)) - 1
cat(sprintf('combined network has size %d at precision 0.5\n', net.size))

net.df$precision <- prec

tmp.net.df <- net.df[1:(net.size*2), ]
for (i in 3:ncol(tmp.net.df)) {
  tmp.net.df[, i] <- round(tmp.net.df[, i], 3)
}
write.table(tmp.net.df, file=output.file, sep='\t', row.names=FALSE, quote=FALSE)

