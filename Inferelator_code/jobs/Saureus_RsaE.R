PARS$input.dir <- '../Input_files/Saureus_RsaE/'
PARS$exp.mat.file <- 'staph_exp_0810_2019.RData'
PARS$tf.names.file <- 'tfs_plus_RsaE_staph.csv'
#PARS$meta.data.file <- 'meta.data.RData'
PARS$priors.file <- 'Staph_RsaE_inf_GS_0810_2019_v2.csv'
PARS$gold.standard.file <- 'Staph_RsaE_inf_GS_0810_2019_v2.csv'

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

PARS$save.to.dir <- 'output/RsaE_Saureus'
