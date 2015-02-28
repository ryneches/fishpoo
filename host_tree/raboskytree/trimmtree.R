# Notes on how I generated the trimmed tree

library(ape)

rabosky <- read.tree("raboskytree/Rabosky_et_al_timetree.tre")
fishnames <- scan("raboskytree/fish_names.txt", what='character')
rabosky_trimmed <- drop.tip( rabosky, setdiff( rabosky$tip.label, fishnames ) )

rabosky_trimmed$tip.label[1] <- "Neolamprologus_buescheri"
rabosky_trimmed$tip.label[5] <- "Cyprichromis_coloratus"
rabosky_trimmed$tip.label[9] <- "Reganochromis_calliurus"
rabosky_trimmed$tip.label[11] <- "Trematochromis_benthicola"

phase1_otus <- read.table( file="fishpoo/fishpoo_open_reference_otus/fishpoo_phase1_otus_vs_host.csv", sep=',')

library('geiger')

td <- treedata( rabosky_trimmed, phase1_otus )

library('phytools')

#phase1_pca <- phyl.pca( rabosky_trimmed, phase1_otus, method="BM", mode="cov")
phase1_pca <- phyl.pca( rabosky_trimmed, phase1_otus )
