PARS$input.dir <- '../Input_files/Ecoli_CopraRNA_sRNA_priors'
PARS$exp.mat.file<- 'm3_unaveraged_sRNAfiltered.RData'
PARS$tf.names.file <- 'regulators_names_ryhB_06.csv'  
PARS$priors.file <- 'ryhB.copra.top15plus.csv'
PARS$gold.standard.file <- 'ryhB.copra.top15plus.csv'


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

PARS$save.to.dir <- paste('output/ryhB_silico_3', PARS$method, PARS$prior.weight, PARS$perc.tp, PARS$perc.fp, sep='_')
