
PARS$input.dir <- '../Input_files/Ecoli_shuffled_sRNA_priors'
PARS$exp.mat.file <- 'm3_unaveraged_sRNAfiltered.RData'
PARS$tf.names.file <- 'tfs_plus_sRNAs_strong.csv'
#PARS$meta.data.file <- 'meta.data.RData'
PARS$priors.file <- 'trn_shuffled_sRNA_priors_7.csv'
PARS$gold.standard.file <- 'trn_shuffled_sRNA_priors_7.csv'

PARS$num.boots <- 20
PARS$cores <- 8

PARS$delT.max <- 60
PARS$delT.min <- 15
PARS$tau <- 15

PARS$perc.tp <- 100
PARS$perm.tp <- 1
PARS$perc.fp <- 0
PARS$perm.fp <- 1

PARS$eval.on.subset <- TRUE

PARS$method <- 'BBSR'
PARS$use.tfa <- TRUE
PARS$prior.weight <- 1.1

PARS$output.summary <- FALSE

PARS$save.to.dir <- 'output/Ecoli_shuffled_sRNA_priors7'
