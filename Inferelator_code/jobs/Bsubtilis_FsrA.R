PARS$input.dir <- '../Input_files/Bsubtilis_FsrA'
PARS$exp.mat.file <- 'bsu_exp_mat.RData'
PARS$tf.names.file <- 'bsu_tf_fsrA_names.csv'
#PARS$meta.data.file <- 'meta.data.RData'
PARS$priors.file <- 'bsu_fsrA_priors1112.csv'
PARS$gold.standard.file <- 'bsu_fsrA_priors1112.csv'

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

PARS$save.to.dir <- 'output/FsrA_Bsubtilis'
