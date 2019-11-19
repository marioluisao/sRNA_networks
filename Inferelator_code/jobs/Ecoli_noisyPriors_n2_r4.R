
PARS$input.dir <- '/Users/mario/Documents/NYU/Inferelator-tfa/Inferelator_2015.08.05/input_ecoli_paper/'
PARS$exp.mat.file <- 'm3_unaveraged_sRNAfiltered.RData'
PARS$tf.names.file <- 'tfs_plus_sRNAs_strong.csv'
#PARS$meta.data.file <- 'meta.data.RData'
PARS$priors.file <- 'inf_noisy_sRNA_priors_1_n2_r4.csv'
PARS$gold.standard.file <- 'inf_noisy_sRNA_priors_1_n2_r4.csv'

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

PARS$save.to.dir <- 'output/eco_sRNA_strong_m3F_unaveraged_noisyPriors_n2_r4'
