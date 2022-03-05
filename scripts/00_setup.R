################################################################################
#  _______  _____  _____ ___ ____  _   _ __  __ _____ ____  
# | ____\ \/ / _ \|  ___|_ _/ ___|| | | |  \/  | ____|  _ \ 
# |  _|  \  / | | | |_   | |\___ \| |_| | |\/| |  _| | | | |
# | |___ /  \ |_| |  _|  | | ___) |  _  | |  | | |___| |_| |
# |_____/_/\_\___/|_|   |___|____/|_| |_|_|  |_|_____|____/ 
#
################################################################################
# R script to globally set up project 
# arthur.escalas@gmail.com
################################################################################


# ==============================================================================
# cleaning memory
# ==============================================================================

rm(list = ls()) 

# ==============================================================================
# Define the directories of the project
# create objects corresponding to the directories and create directories locally
# ==============================================================================

# ------------------------------------------------------------------------------
# identify the project name, find where we are and create local directories

dir_project <- getwd()
proj_name   <- unlist(strsplit(dir_project, 
                               split = 'marbec_exofishmed', 
                               fixed = TRUE))[2]
dir_project_global <- unlist(strsplit(dir_project, 
                                      split = proj_name, 
                                      fixed = TRUE))[1]

dir_scripts   <- paste0(dir_project, "/scripts/")
dir_functions <- paste0(dir_scripts, "/functions/")

dir_data     <- paste0(dir_project, "/data/")

dir_analyses <- paste0(dir_project, "/analyses/")
dir_alpha_div <- paste0(dir_analyses, "01_alpha_diversity/")
dir_composition  <- paste0(dir_analyses, "02_microbiome_composition/")
dir_compo_res    <- paste0(dir_analyses, "03_analyze_composition_results/")
dir_dissimilarity      <- paste0(dir_analyses, "04_analyze_dissimilarity/")


dir.create(dir_scripts)

dir.create(dir_data)

dir.create(dir_analyses)

dir.create(dir_alpha_div)
dir.create(dir_composition)
dir.create(dir_compo_res)
dir.create(dir_dissimilarity)

# ------------------------------------------------------------------------------
# Source functions and libs

# source utility functions
for (f in list.files(dir_functions, full.names = T)) { source(f) }

# source package list
source(paste0(dir_scripts, "packages.R"))
install_if_not_there(cran_packages, type = "CRAN")
install_if_not_there(bioc_packages, type = "bioconductor")
sapply(c(cran_packages, bioc_packages), require, character.only = TRUE)


# Pairwise ADONIS

install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library("pairwiseAdonis")


# ------------------------------------------------------------------------------
# Objcts for data anlysis and plots

nms_region <- c("North_Red_Sea", "Levantine_Sea", "Northern_Crete")
nms_siganids <- c("Siganus_rivulatus", "Siganus_luridus")

pch_season <- c(21)

pch_reg <- c(21,23,22) %>% setNames(nms_region)
col_reg  <- c("#FF6347","#63B8FF", "#CDCD00") %>% setNames(nms_region)

# names of the alpha diversity indexes 

index_alpha <- c("taxo_q0", "taxo_q1","phylo_q0", "phylo_q1")
names(index_alpha) <- index_alpha




