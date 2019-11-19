## Bonneau lab
## NYU - Center for Genomics and Systems Biology

# Call this script with a job config file as arguments
# Example call: Rscript inferelator.R jobs/dream4_cfg.R


library('Matrix')

rm(list=ls())
gc()

source('R_scripts/utils.R')
source('R_scripts/design_and_response.R')
source('R_scripts/priors.R')
source('R_scripts/mi_and_clr.R')
source('R_scripts/bayesianRegression.R')
source('R_scripts/men.R')
source('R_scripts/evaluate.R')
source('R_scripts/tfa.R')
source('R_scripts/group_predictors.R')


date.time.str <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
print(date.time.str)
start.proc.time <- proc.time()

# default job parameters
PARS <- list()

PARS$input.dir <- 'input/dream4'

PARS$exp.mat.file <- 'expression.tsv'
PARS$tf.names.file <- 'tf_names.tsv'
PARS$meta.data.file <- NULL
PARS$priors.file <- NULL
PARS$gold.standard.file <- NULL
PARS$leave.out.file <- NULL
PARS$randomize.expression <- FALSE

PARS$job.seed <- 42  # set to NULL if a random seed should be used
PARS$save.to.dir <- NULL
PARS$num.boots <- 20
PARS$max.preds <- 10
PARS$mi.bins <- 10
PARS$cores <- 8

PARS$delT.max <- 110
PARS$delT.min <- 0
PARS$tau <- 45

PARS$perc.tp <- 0
PARS$perm.tp <- 1
PARS$perc.fp <- 0
PARS$perm.fp <- 1
PARS$pr.sel.mode <- 'random'  # prior selection mode: 'random' or 'tf'

PARS$eval.on.subset <- FALSE

PARS$method <- 'BBSR'  # 'BBSR' or 'MEN'
PARS$prior.weight <- 1

PARS$use.tfa <- FALSE


# some of the elastic net parameters that are essentially constants;
# only override in config script if you know what you are doing
PARS$enet.sparseModels <- TRUE    # sparser models
PARS$enet.nCv <- 10               # number of cross-validations
PARS$enet.lambda <- c(0, 1, 100)  # l2 weights
PARS$enet.verbose <- FALSE        # print progress to screen
PARS$enet.plot.it <- FALSE        # generate cross-validation plots
PARS$enet.plot.file.name <- NULL  # file name for plots


# input argument is the job config script which overrides the default parameters
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 1) {
  job.cfg <- args[1]
} else {
  #job.cfg <- 'jobs/bsubtilis_eu_201310_stfa_bbsr_22.R'
  #job.cfg <- 'jobs/bsubtilis_eu_201310_stfa_bbsr_1_tp0_fp0.R'
  #job.cfg <- '/home/ch1421/Projects/Rice/inferelator_jobs/ALL_htseq_intersection-strict_motifprior.R'
  #job.cfg <- '/home/ch1421/Projects/Emily/ILC_inferelator_job.R'
  #job.cfg <- '/home/ch1421/Projects/Kostya/inferelator_jobs/test.R'
  #job.cfg <- 'jobs/bsubtilis_eu_201310_stfa_bbsr_tf_11_tp50_fp0.R'
  #job.cfg <- 'jobs/bsubtilis_us_201310_stfa_bbsr_11_tp100_fp0.R'
  #MARIO job.cfg <- 'jobs/bsubtilits_eu_mario_20150420.R'
  job.cfg <- 'jobs/gcvB.0.R'
  #job.cfg <- 'jobs/bsubtilis_us_test.R'
}

# load job specific parameters from input config file
if (!is.null(job.cfg)) {
  source(job.cfg)
}


# read input data
IN <- read.input(PARS$input.dir, PARS$exp.mat.file, PARS$tf.names.file,
                 PARS$meta.data.file, PARS$priors.file, PARS$gold.standard.file,
                 PARS$leave.out.file, PARS$randomize.expression)


# keep only TFs that are part of the expression data
#IN$tf.names <- IN$tf.names[IN$tf.names %in% rownames(IN$exp.mat)]

# order genes so that TFs come before the other genes
#gene.order <- rownames(IN$exp.mat)
#gene.order <- c(gene.order[match(IN$tf.names, gene.order)],
#                gene.order[which(!(gene.order %in% IN$tf.names))])

#IN$exp.mat <- IN$exp.mat[gene.order, ]
#if (!is.null(IN$priors.mat)) {
#  IN$priors.mat <- IN$priors.mat[gene.order, IN$tf.names]
#}
#if (!is.null(IN$gs.mat)) {
#  IN$gs.mat <- IN$gs.mat[gene.order, IN$tf.names]
#}


# no meta data given - assume all steady state measurements
if (is.null(IN$meta.data)) {
  IN$meta.data <- trivial.meta.data(colnames(IN$exp.mat))
}


# create dummy clusterStack - a real clusterStack is only needed when inferring
# on bi-clusters
clusterStack <- trivial.cluster.stack(IN$exp.mat)


# set the random seed
if(!is.null(PARS$job.seed)) {
  set.seed(PARS$job.seed, "Mersenne-Twister", "Inversion")
  cat("RNG seed has been set to ", PARS$job.seed, "\n")
} else {
  ignore <- runif(1)
}
SEED <- .Random.seed

if(is.null(PARS$save.to.dir)) {
  PARS$save.to.dir <- file.path(PARS$input.dir, date.time.str)
}
cat("Output dir:", PARS$save.to.dir, "\n")
if (!file.exists(PARS$save.to.dir)){
  dir.create(PARS$save.to.dir, recursive=TRUE)
}


##  .-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.
# create design and response matrix
##  .-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.

cat("Creating design and response matrix ... ")
des.res <- design.and.response(IN$meta.data, IN$exp.mat, PARS$delT.min,
                               PARS$delT.max, PARS$tau)

IN$final_response_matrix <- des.res$final_response_matrix
IN$final_design_matrix <- des.res$final_design_matrix
resp.idx <- des.res$resp.idx
cat("done.\n")

if (!all(apply(resp.idx, 1, identical, resp.idx[1,]))) {
    stop('This version of the Inferelator does not support biclusters. Sorry.')
}


##  .-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.
# set up the bootstrap permutations
##  .-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.


cat("Setting up bootstrap permutations ... ")
IN$bs.pi <- matrix(0, nrow=PARS$num.boots, ncol=ncol(resp.idx))
if (PARS$num.boots == 1) {
  IN$bs.pi[1, ] <- resp.idx[1, ]
} else {
  for (bootstrap in 1:PARS$num.boots) {
    IN$bs.pi[bootstrap, ] <- resp.idx[1, sample(ncol(resp.idx), replace=TRUE)]
  }
}
cat("done.\n")


##  .-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.
# parse priors parameters and set up priors list
##  .-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.


cat("Setting up priors list ... ")
IN$priors <- getPriors(IN$exp.mat, IN$tf.names, IN$priors.mat, IN$gs.mat,
                       PARS$eval.on.subset, PARS$job.seed, PARS$perc.tp,
                       PARS$perm.tp, PARS$perc.fp, PARS$perm.fp,
                       PARS$pr.sel.mode)
cat("done.\n")


##  .-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.
# TFA specific initialization
##  .-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.

if(PARS$use.tfa) {

  IN$tf.activities <- list()
  cat("Setting up TFA specific response matrix ... ")
  des.res <- design.and.response(IN$meta.data, IN$exp.mat, PARS$delT.min,
                                 PARS$delT.max, PARS$tau/2)

  IN$final_response_matrix_halftau <- des.res$final_response_matrix
  cat("done.\n")
}


##  .-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.
# main loop
##  .-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.***.-.-.
IN$grouped.pred <- list()
for (prior.name in names(IN$priors)) {
  cat('Method:', PARS$method, '\nWeight:', PARS$prior.weight, '\nPriors:',
      prior.name, '\n')
  prior <- as.matrix(IN$priors[[prior.name]])

  # estimate transcription factor activities
  if(PARS$use.tfa) {
    IN$tf.activities[[prior.name]] <- tfa(prior, IN$final_design_matrix, IN$final_response_matrix_halftau)
  }

  des.mat <- IN$final_design_matrix[IN$tf.with.expr, ]
  if(PARS$use.tfa) {
      des.mat <- IN$tf.activities[[prior.name]][IN$tf.names, ]
  }

  cat("group predictors\n")
  gp.out <- group.predictors(des.mat, prior, IN$gs.mat, IN$bs.pi, cor.th=0.99)
  IN$grouped.pred[[prior.name]] <- gp.out

  # set the prior weights matrix
  no.pr.weight <- 1
  if (sum(prior != 0) > 0) {
    if (PARS$prior.weight == no.pr.weight) {
      warning(paste('Priors present, but they will not be used, because \
                    PARS$prior.weight is set to ', no.pr.weight, '.', sep=''),
              immediate. = TRUE)
    }
    if (PARS$method == 'BBSR') {
      no.pr.weight <- 1 / PARS$prior.weight
    }
  }
  #weights.mat <- matrix(no.pr.weight, nrow(IN$exp.mat), length(IN$tf.names))
  #weights.mat[prior != 0] <- PARS$prior.weight
  weights.mat <- IN$grouped.pred[[prior.name]]$prior.mat * 0 + no.pr.weight
  weights.mat[IN$grouped.pred[[prior.name]]$prior.mat != 0] <- PARS$prior.weight


  betas <- list()
  betas.resc <- list()
  for (bootstrap in 1:PARS$num.boots) {
    cat("Bootstrap", bootstrap, "of", PARS$num.boots, "\n")

    # set up bootstrap specific design and response
    #X <- IN$final_design_matrix[, IN$bs.pi[bootstrap, ]]
    X <- IN$grouped.pred[[prior.name]]$des.mat[, IN$bs.pi[bootstrap, ]]
    Y <- IN$final_response_matrix[, IN$bs.pi[bootstrap, ]]

    if (nrow(X) > 6000) {
      #X <- X[IN$tf.names, ]  # speeds up MI calculation for large datasets
      X <- X[IN$grouped.pred[[prior.name]]$tf.names, ]
    }

    #if(PARS$use.tfa) {
    #  X <- IN$tf.activities[[prior.name]][, IN$bs.pi[bootstrap, ]]
    #}

    # fill mutual information matrices
    cat("Calculating MI\n")
    Ms <- mi(t(Y), t(X), nbins=PARS$mi.bins, cpu.n=PARS$cores)
    diag(Ms) <- 0
    cat("Calculating Background MI\n")
    Ms_bg <- mi(t(X), t(X), nbins=PARS$mi.bins, cpu.n=PARS$cores)
    diag(Ms_bg) <- 0

    # get CLR matrix
    cat("Calculating CLR Matrix\n")
    clr.mat = mixedCLR(Ms_bg,Ms)
    dimnames(clr.mat) <- list(rownames(Y), rownames(X))
    #clr.mat <- clr.mat[, IN$tf.names]
    clr.mat <- clr.mat[, IN$grouped.pred[[prior.name]]$tf.names]

    # DREAM8 induced change:
    #for (tf1 in IN$tf.names) {
    #  for (tf2 in IN$tf.names) {
    #    if (tf1 != tf2) {
    #      #if (clr.mat[tf1, tf2] > clr.mat[tf2, tf1]) {
    #      if (Ms[tf1, tf2] > Ms[tf2, tf1]) {
    #        clr.mat[tf2, tf1] <- min(clr.mat)
    #      } else if (Ms[tf1, tf2] < Ms[tf2, tf1]) {
    #        clr.mat[tf1, tf2] <- min(clr.mat)
    #      }
    #    }
    #  }
    #}

    # get the sparse ODE models
    #X <- X[IN$tf.names, ]
    X <- X[IN$grouped.pred[[prior.name]]$tf.names, ]
    cat('Calculating sparse ODE models\n')
    if (PARS$method == 'BBSR') {
      #x <- BBSR(X, Y, clr.mat, PARS$max.preds, no.pr.weight, weights.mat,
      #          prior, PARS$cores)
      x <- BBSR(X, Y, clr.mat, PARS$max.preds, no.pr.weight, weights.mat,
                IN$grouped.pred[[prior.name]]$prior.mat, PARS$cores)
    }
    if (PARS$method == 'MEN' ) {
      stop('MEN currently not tested - remove this line and proceed at own risk')
      x <- mclapply(1:nrow(Y), callMEN, Xs=X, Y=Y,
                    clr.mat=clr.mat, nS=PARS$max.preds, nCv=PARS$enet.nCv,
                    lambda=PARS$enet.lambda, verbose=PARS$enet.verbose,
                    plot.it=PARS$enet.plot.it,
                    plot.file.name=PARS$enet.plot.file.name,
                    weights.mat=weights.mat, no.pr.val=no.pr.weight,
                    mc.cores=PARS$cores)
    }
    cat('\n')

    # our output will be a list holding two matrices: betas and betas.resc
    bs.betas <- Matrix(0, nrow(Y), nrow(X),
                       dimnames=list(rownames(Y), rownames(X)))
    bs.betas.resc <- Matrix(0, nrow(Y), nrow(X),
                            dimnames=list(rownames(Y), rownames(X)))
    for (res in x) {
      bs.betas[res$ind, res$pp] <- res$betas
      bs.betas.resc[res$ind, res$pp] <- res$betas.resc
    }
    betas[[bootstrap]] <- bs.betas
    betas.resc[[bootstrap]] <- bs.betas.resc

  }  # end bootstrap for loop

  res.file <- paste(PARS$save.to.dir, "/betas_", prior.name, "_", PARS$prior.weight, ".RData", sep="")
  save(betas, betas.resc, file = res.file)

  # rank-combine the rescaled betas (confidence scores) of the bootstraps
  confs.file <- sub('/betas_', '/combinedconf_', res.file)
  comb.confs <- Matrix(0, nrow(betas.resc[[1]]), ncol(betas.resc[[1]]),
                       dimnames=dimnames(betas.resc[[1]]))
  for (beta.resc in betas.resc) {
    comb.confs <- comb.confs + rank(as.matrix(beta.resc), ties.method='average')
  }
  save(comb.confs, file=confs.file)

}  # end prior.name loop

PROCTIME <- proc.time() - start.proc.time
save(PARS, IN, SEED, PROCTIME, file = paste(PARS$save.to.dir, "/params_and_input.RData", sep=""))

# generate network report
source('R_scripts/net_report_new.R')
Sys.sleep(2)
for (ccf in list.files(PARS$save.to.dir, pattern='combinedconf_*', full.names=TRUE)) {
  net.report(normalizePath(ccf))
}

# this part does not work for grouped predictors
#if (!is.null(IN$gs.mat)) {
#  cat('Using gold standard to evaluate results. Evaluate on subset is set to', PARS$eval.on.subset, '. \n')
#  summarizeResults(PARS$save.to.dir, PARS$eval.on.subset)
#}
