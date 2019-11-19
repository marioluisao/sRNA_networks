
PARS$input.dir <- 'input/eu'
PARS$exp.mat.file <- 'exp.mat.RData'
PARS$tf.names.file <- 'tf.names.RData'
PARS$meta.data.file <- 'meta.data.RData'
PARS$priors.file <- 'prior.mat.RData'
PARS$gold.standard.file <- 'gs.mat.RData'

PARS$num.boots <- 5
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

PARS$output.summary <- TRUE

PARS$save.to.dir <- paste('output/bsubtilis_eu_TFA', PARS$use.tfa, PARS$method, PARS$prior.weight, PARS$perc.tp, PARS$perc.fp, sep='_')

