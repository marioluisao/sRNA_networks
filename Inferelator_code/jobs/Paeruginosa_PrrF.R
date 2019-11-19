
PARS$input.dir <- '../Input_files/Paeruginosa_PrrF'
PARS$exp.mat.file <- 'pa_exp_mat.RData'
PARS$tf.names.file <- 'pa_tf_prrF_names.csv'
#PARS$meta.data.file <- 'meta.data.RData'
PARS$priors.file <- 'pa_PrrF_priors0712.csv'
PARS$gold.standard.file <- 'pa_PrrF_priors0712.csv'

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

PARS$save.to.dir <- 'output/PrrF_Paeruginosa'
