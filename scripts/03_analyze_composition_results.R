################################################################################
#  _______  _____  _____ ___ ____  _   _ __  __ _____ ____  
# | ____\ \/ / _ \|  ___|_ _/ ___|| | | |  \/  | ____|  _ \ 
# |  _|  \  / | | | |_   | |\___ \| |_| | |\/| |  _| | | | |
# | |___ /  \ |_| |  _|  | | ___) |  _  | |  | | |___| |_| |
# |_____/_/\_\___/|_|   |___|____/|_| |_|_|  |_|_____|____/ 
#
#
################################################################################
# R script to summarize results of microbiome composition analyses
# arthur.escalas@gmail.com
################################################################################


############################### LOAD DATA ######################################

dir_save <- dir_compo_res

# Load the phyloseq objects with diversity estimates ----

# the list of full phyloseq at different ranks
ls_ps_full <- readRDS(paste0(dir_data, "list_glom_full_phyloseq_per_taxa_rank.rds"))
ls_ps_full$ASV <- readRDS(paste0(dir_data, "phyloseq_object_cleaned_with_tree_rarefied_with_alpha_div.rds"))


# list of dissimilarity estimates for full dataset
ls_diss_full <- readRDS(paste0(dir_data, "list_all_diversity_estimates_per_ranks_full_dataset.rds"))
ls_diss_full$ASV <- readRDS(paste0(dir_data, "beta_diversity_estimates_full_dataset_ASV.rds"))


# list of fish core phyloseq objects
ls_ps_core <- readRDS(paste0(dir_data, "list_glom_fish_phyloseq_per_taxa_rank.rds"))
tmp <- readRDS(paste0(dir_data, "list_phyloseq_object_fish_core_with_alpha_diversity.rds"))

ls_ps_core <- lapply(names(ls_ps_core), function(x) {
  ls_ps_core[[x]]$ASV <- tmp[[x]]
  ls_ps_core[[x]]
}) %>% setNames(names(ls_ps_core))


# list of dissimilarity estimates for fish core dataset
ls_diss_fish <- readRDS(paste0(dir_data, "list_all_diversity_estimates_per_ranks_fish_dataset.rds"))
tmp <- readRDS(paste0(dir_data, "list_alpha_beta_diversity_estimates_fish_ASV.rds"))

ls_diss_fish <- lapply(names(ls_diss_fish), function(x) {
  ls_diss_fish[[x]]$ASV <- tmp[[x]]
  ls_diss_fish[[x]]
}) %>% setNames(names(ls_diss_fish))


# Names of ranks

nms_ranks <- c("Phylum", "Class", "Order", "Family", "Genus", "Species", "ASV")
nms_ranks_5 <- c("Phylum", "Class", "Order", "Family", "Genus")
nms_ranks_3 <- c("Phylum", "Family", "ASV")


# change the levels of factors for plots

tab <- meta(ls_ps_full$ASV)
tab$region <- factor(tab$region, levels = nms_region)
tab$season <- factor(tab$season, levels = c("Spring", "Autumn"))
nms_reg_seas <- c("North_Red_Sea_Spring", "North_Red_Sea_Autumn", 
                  "Levantine_Sea_Spring", "Levantine_Sea_Autumn", 
                  "Northern_Crete_Spring", "Northern_Crete_Autumn")
tab$region_season <- factor(tab$region_season, levels = nms_reg_seas)
tab$sample_type <- factor(tab$sample_type, levels = c("Water","Sediment",
                                                      "Turf","Algae","Seagrass",
                                                      "Fish"))
nms_siganids <- names(ls_ps_core)

# objects for plotting ----

nms_reg_seas_plot <- c("North Red Sea\nSpring", "North Red Sea\nAutumn", 
                       "Levantine Sea\nSpring", "Levantine Sea\nAutumn", 
                       "Northern Crete\nSpring", "Northern Crete\nAutumn")
nms_reg_seas <- c("North_Red_Sea_Spring", "North_Red_Sea_Autumn", 
                  "Levantine_Sea_Spring", "Levantine_Sea_Autumn", 
                  "Northern_Crete_Spring", "Northern_Crete_Autumn")
cols_regions <- c("#FF6347", "#FFFFFF", "#63B8FF", "#FFFFFF", "#CDCD00", "#FFFFFF") %>% 
  setNames(nms_reg_seas)
cols_border <- c("#FF6347", "#FF6347", "#63B8FF", "#63B8FF", "#CDCD00", "#CDCD00") %>% 
  setNames(nms_reg_seas)
nms_index <- c("q0", "q1")



####################### SUMMARIZE THE PERMANOVA RESULTS ########################


# Load the permanova results from other folders ----

path_folds <- dir_composition
nms_folds <- c("algae", "turf", "water", "sediment", "s_luridus", 
               "s_rivulatus")

ls_permanova <- list()
for (f in nms_folds) {
  file <- list.files(paste0(path_folds, "/", f, "/stats/"),
                     pattern = "table_results_permanova.csv")
  if (length(file) != 0) {
    if (file.exists(paste0(path_folds,"/", f, "/stats/",file))) {
      ls_permanova[[f]] <- read.csv(paste0(path_folds,"/", f, "/stats/",file))
    }
  }
}


# Format the result table ----

df_summary_stats <- ls_permanova %>% reformat_as_df("compartment") %>% 
  filter(!factor %in% c("Residual","Total")) %>% 
  mutate(factor = factor(factor, levels = c("region","season","region:season"))) %>% 
  filter(rank %in% nms_ranks_3) %>% 
  arrange(compartment, data_type, rank,factor) %>% 
  dplyr::select(compartment, factor, rank, data_type, R2, F.value, p.value)

df_summary_stats$data_type[df_summary_stats$data_type == "q0"] <- "taxo_q0"
df_summary_stats$data_type[df_summary_stats$data_type == "q1"] <- "taxo_q1"

df_summary_signif <- df_summary_stats %>% 
  filter(! compartment %in% c("fish", "siganidae")) %>% 
  tidyr::separate("data_type", into = c("facet", "order"), sep = "_", remove = FALSE) %>% 
  mutate(facet = factor(facet, levels = c("taxo", "phylo")),
         compartment = factor(compartment, levels = c("algae", "sediment","turf",
                                                      "water","s_luridus","s_rivulatus"))) %>% 
  group_by(facet, compartment, factor) %>% 
  summarise(num_signif = sum(p.value < 0.05),
            avg_p = mean(p.value) %>% abs(),
            avg_f = mean(F.value[F.value > 0]) %>% abs(),
            avg_r2 = mean(R2) %>% abs()) %>% data.frame()

# export results ----

write.csv(df_summary_stats, row.names = FALSE,
          file = paste0(dir_save, "table_summary_results_permanova.csv"))
write.csv(df_summary_signif, row.names = FALSE,
          file = paste0(dir_save, "table_summary_results_significant_permanova.csv"))



# Arrange the table for supplementary material ----
# 
# tab <- df_summary_stats %>% filter(factor == "treatment") %>% data.frame()
# 
# 
# tab$compartment <- gsub("_samples", "", tab$compartment)
# tab$compartment <- gsub("sediment_", "", tab$compartment)
# tab$compartment <- gsub("_", " ", tab$compartment)
# tab$compartment <- firstup(tab$compartment)
# 
# tab$q <- gsub("q0", "Composition", tab$q)
# tab$q <- gsub("q1", "Structure", tab$q)
# tab$div_facet <- gsub("phylo", "Phylogenetic", tab$div_facet)
# tab$div_facet <- gsub("taxo", "Taxonomic", tab$div_facet)
# tab$factor <- gsub("treatment", "Treatment", tab$factor)
# tab$factor <- gsub("_", " ", tab$factor)
# tab$rank <- factor(tab$rank, levels = c("Phylum","Family","OTU"))
# 
# tab <- tab %>% 
#   arrange(desc(compartment), div_facet, q) %>% 
#   dplyr::select(compartment, div_facet, q, rank,  p.value) %>% data.frame()
# names(tab) <- c("Compartment","Diversity facet", "Diversity component", "Rank", 
#                 "p-value")
# 
# out <- tab %>% pivot_wider(names_from = c("Diversity facet", "Diversity component", "Rank"),
#                            values_from = c("p-value")) %>% data.frame()
# 
# write.csv(out, file = paste0(dir_save, "table_results_PERMANOVA_for_supplementary.csv"),
#           row.names=FALSE)
# 


####################### SUMMARIZE THE PERMDISP RESULTS ########################


# Load the permanova results from other folders ----

path_folds <- dir_composition
nms_folds <- c("algae", "turf", "water", "sediment", "s_rivulatus", "s_luridus")

ls_permdisp <- list()
for (f in nms_folds) {
  file <- list.files(paste0(path_folds, "/", f, "/stats/"),
                     pattern = "res_permdisp.rds")
  if (length(file) != 0) {
    if (file.exists(paste0(path_folds,"/", f, "/stats/",file))) {
      ls_permdisp[[f]] <- readRDS(paste0(path_folds,"/", f, "/stats/",file))
    }
  }
}


## result table of  the test -----

ls_permdisp_tab <- list()
for (f in nms_folds) {
  file <- list.files(paste0(path_folds, "/", f, "/stats/"),
                     pattern = "table_results_permdisp.csv")
  if (length(file) != 0) {
    if (file.exists(paste0(path_folds,"/", f, "/stats/",file))) {
      ls_permdisp_tab[[f]] <- read.csv(paste0(path_folds,"/", f, "/stats/",file))
    }
  }
}


# Format the result table ----

# average distance to centroid

df_dist_centroid <- lapply(ls_permdisp, function(X) {
  lapply(X, function(x) {
    lapply(x, function(xx) {
      cbind(avg = tapply(xx$distances, xx$group, mean),
            sd = tapply(xx$distances, xx$group, sd)) %>% 
        data.frame() %>% rownames_to_column("region")
    }) %>% reformat_as_df("index")
  }) %>% reformat_as_df("rank")
}) %>% reformat_as_df("compartment") %>% 
  pivot_wider(names_from = "region", values_from = c("avg", "sd")) %>% 
  mutate(merged_col = paste(rank,index,compartment, sep ="_")) %>% 
  dplyr::select(-rank, -index, -compartment)

df_summary_stats <- ls_permdisp_tab %>% reformat_as_df("compartment") %>% 
  filter(!factor %in% c("Residuals","Total")) %>% 
  dplyr::select(data_type, compartment, rank, q, F.value, p.value) %>% 
  mutate(compartment = factor(compartment, levels = c("algae", "sediment","turf",
                                                      "water","s_luridus","s_rivulatus"))) %>% 
  arrange(data_type,compartment, rank)  %>% 
  mutate(merged_col = paste(rank,data_type,compartment, sep ="_"))
  

df_summary_stats <- df_summary_stats %>% left_join(df_dist_centroid, by = "merged_col")


df_summary_stats$data_type[df_summary_stats$data_type == "q0"] <- "taxo_q0"
df_summary_stats$data_type[df_summary_stats$data_type == "q1"] <- "taxo_q1"

df_summary_signif <- df_summary_stats %>% 
  filter(! compartment %in% c("fish", "siganidae")) %>% 
  tidyr::separate("data_type", into = c("facet", "order"), sep = "_", remove = FALSE) %>% 
  mutate(facet = factor(facet, levels = c("taxo", "phylo")),
         compartment = factor(compartment, levels = c("algae", "sediment","turf",
                                                      "water","s_luridus","s_rivulatus"))) %>% 
  group_by(rank, compartment) %>% 
  summarise(num_signif = sum(p.value < 0.05),
            avg_p = mean(p.value) %>% abs(),
            avg_f = mean(F.value[F.value > 0]) %>% abs()) %>% data.frame()

# export results ----

write.csv(df_summary_stats, row.names = FALSE,
          file = paste0(dir_save, "table_summary_results_permdisp.csv"))
write.csv(df_summary_signif, row.names = FALSE,
          file = paste0(dir_save, "table_summary_results_significant_permdisp.csv"))






################ SUMMARIZE THE PERMDISP RESULTS FOR FUNCTIONS ##################


# Load the permanova results from other folders ----

path_folds <- dir_composition

ls_permdisp <-  readRDS(paste0(path_folds,"/functional/stats/res_permdisp_scfa.rds"))

nms_func_lev <- names(ls_permdisp)


## result table of  the test -----

permdisp_tab <- read.csv(paste0(path_folds,"/functional/stats/table_results_permdisp.csv"))
  


# Format the result table ----

# average distance to centroid

df_dist_centroid <- lapply(ls_permdisp, function(X) {
  lapply(X, function(x) {
    lapply(x, function(xx) {
      cbind(avg = tapply(xx$distances, xx$group, mean),
            sd = tapply(xx$distances, xx$group, sd)) %>% 
        data.frame() %>% rownames_to_column("region")
    }) %>% reformat_as_df("index")
  }) %>% reformat_as_df("rank")
}) %>% reformat_as_df("compartment") %>% 
  pivot_wider(names_from = "region", values_from = c("avg", "sd")) %>% 
  mutate(merged_col = paste(rank,index,compartment, sep ="_")) %>% 
  dplyr::select(-rank, -index, -compartment)

df_summary_stats <- permdisp_tab %>%
  filter(!factor %in% c("Residuals","Total")) %>% 
  dplyr::select(data_type, species, rank, q,div_facet, F.value, p.value) %>% 
  arrange(data_type,species, rank)  %>% 
  mutate(merged_col = paste(species, div_facet, rank, sep ="_"))


df_summary_stats <- df_summary_stats %>% left_join(df_dist_centroid, by = "merged_col")


df_summary_stats$data_type[df_summary_stats$data_type == "q0"] <- "taxo_q0"
df_summary_stats$data_type[df_summary_stats$data_type == "q1"] <- "taxo_q1"

df_summary_signif <- df_summary_stats %>% 
  group_by(species, rank) %>% 
  summarise(num_signif = sum(p.value < 0.05),
            avg_p = mean(p.value) %>% abs(),
            avg_f = mean(F.value[F.value > 0]) %>% abs()) %>% data.frame()

# export results ----

write.csv(df_summary_stats, row.names = FALSE,
          file = paste0(dir_save, "table_summary_results_permdisp_functional.csv"))
write.csv(df_summary_signif, row.names = FALSE,
          file = paste0(dir_save, "table_summary_results_significant_permdisp_functional.csv"))



###################### MAKE THE PCOA FIGURE FOR THE PAPER #############################

path_folds <- dir_composition

nms_folds <- c("algae", "sediment", "turf", "water","s_luridus","s_rivulatus")

# laod the pcoa results ----

ls_pcoa <- list()
for (f in nms_folds) {
  
  file <- paste0(path_folds,"/", f, "/pcoa/list_pcoa.rds")
  ls_pcoa[[f]] <- readRDS(file)
}
names(ls_pcoa)[5:6] <- nms_siganids[2:1]


# Make the plot for different ranks -----

for (rk in  nms_ranks_3) {
  
   jpeg(paste0(dir_save, "figure_pcoa_", rk, ".jpeg"),
      height = 15 , width = 32, unit = "cm", res = 600)
  
  par(mfrow = c(2,4), las = 1, mar = c(3,3,2,2), oma = c(3,1,0,1),
      mgp = c(2.2,.5,0), tcl = -0.3, font.lab = 2)
  cexpt <- 1
  # for ecosystem compartments ----
  
  idx <- "q1"
  
  for (X in names(ls_pcoa)[1:4]) {

      axes <- c(1,2)
      tmp <- ls_pcoa[[X]][[rk]][[idx]]
      dat <- tmp$vectors[, axes]
      eig <- tmp$values$Eigenvalues[tmp$values$Eigenvalues > 0]
      expl_var <- round(eig / sum(eig) * 100, 1)
      
      # create color and pch vector ----
      metadata <- meta(ls_ps_full$Phylum)
      nms_sples <- row.names(dat)
      metaplot <- metadata %>% filter(sample_id_fastq %in% nms_sples)
      
      col_border <- cols_border[as.character(metaplot$region_season)]
      colvec <- cols_regions[as.character(metaplot$region_season)]
      
      # do the plot ----
      xlim <- setrge(dat[,1],10) 
      ylim <- setrge(dat[,2],10)
      plot(dat, type = "n", xlim = xlim, ylim = ylim, main = "",
           xlab = paste0("Axis ", axes[1], " (", expl_var[axes[1]], "%)"),
           ylab = paste0("Axis ", axes[2], " (", expl_var[axes[2]], "%)"))
      abline(h = 0, lty = 2, col = "grey")
      abline(v = 0, lty = 2, col = "grey")
      points(dat, pch = 21, 
             bg = colvec, 
             col = col_border, cex = cexpt)
      mtext(firstup(X), 3, font = 2, cex = 1)

  }

  # for fishes -----

  idx <- "taxo_q1"

  for (X in nms_siganids[1:2]) {
    
    for (axes in list(c(1,2), c(3,4))) {
    # axes <- c(1,2)
    tmp <- ls_pcoa[[X]][[rk]][[idx]]
    dat <- tmp$vectors[, axes]
    eig <- tmp$values$Eigenvalues[tmp$values$Eigenvalues > 0]
    expl_var <- round(eig / sum(eig) * 100, 1)
    
    # create color and pch vector ----
    metadata <- meta(ls_ps_core$Siganidae$Phylum)
    nms_sples <- row.names(dat)
    metaplot <- metadata %>% filter(sample_id_fastq %in% nms_sples)
    
    col_border <- cols_border[as.character(metaplot$region_season)]
    colvec <- cols_regions[as.character(metaplot$region_season)]
    
    # do the plot ----
    xlim <- setrge(dat[,1],10) 
    ylim <- setrge(dat[,2],10)
    plot(dat, type = "n", xlim = xlim, ylim = ylim, main = "",
         xlab = paste0("Axis ", axes[1], " (", expl_var[axes[1]], "%)"),
         ylab = paste0("Axis ", axes[2], " (", expl_var[axes[2]], "%)"))
    abline(h = 0, lty = 2, col = "grey")
    abline(v = 0, lty = 2, col = "grey")
    points(dat, pch = 21,  
           bg = colvec, 
           col = col_border, cex = cexpt)
    mtext(firstup(gsub("_", " ", X)), 3, font = 4, cex = 1)
    }
  }
    dev.off()
    
}
### -----
  




#### Abundance of taxa (CLR transformed) of the core in each species ###########


dir_plot <- paste0(dir_save, "abundance_boxplots/")
dir.create(dir_plot)


ls_ps_core_clr <- lapply(ls_ps_core, function(X) {
  lapply(X, function(XX) {
    XX %>% #merge_samples(group = "region") %>% 
      microbiome::transform(transform = "clr")
  })
})

ls_ps_core_compo <- lapply(ls_ps_core, function(X) {
  lapply(X, function(XX) {
    XX %>% #merge_samples(group = "region") %>% 
      microbiome::transform(transform = "compositional")
  })
})


# objects for analyse ----

nms_rks <- nms_ranks

ls_dat <- list()
ls_kw <- list()

for (X in nms_siganids[1:2]) {
  
  for (rk in nms_ranks_5) {
    
    dat <- ls_ps_core_clr[[X]][[rk]]
    # dat <- ls_ps_core_compo[[X]][[rk]]
    
    # objects for plot and model fit ----
    
    otu <- dat@otu_table@.Data %>% as.matrix()
    row.names(otu) <- dat@tax_table[,rk]
    otu <- otu[! row.names(otu) == "Other",]
    meta_tmp <- meta(dat)
    meta_tmp$region <- factor(meta_tmp$region, levels = nms_region)
    fact <-  meta_tmp$region
    nms_taxa <- row.names(otu)
    
    # test differences between regions for each taxa ----
    

    ls_kw[[X]][[rk]] <- lapply(nms_taxa, function(x) {
      kruskal.test(otu[x,] ~ meta_tmp$region)
    }) %>% setNames(nms_taxa)
    
    # export the data ----
    
    nms_taxa <- sort(rowSums(otu), decreasing = TRUE) %>% names() %>% sort()
    otu <- otu[nms_taxa,]
    ls_dat[[X]][[rk]]$otu <- otu
    ls_dat[[X]][[rk]]$fact <- meta_tmp$region

  }
}

saveRDS(ls_dat, file = paste0(dir_plot, "list_ANOVA_data.rds"))
saveRDS(ls_kw, file = paste0(dir_plot, "list_KW_tests.rds"))


# Analysis of difference of abundance between regions ----


df_kw <- lapply(ls_kw, function(X) {
  out <- lapply(X, function(x) {
    lapply(x, function(fit) {
      data.frame(stat = fit$statistic, p_value = fit$p.value)  
    }) %>% reformat_as_df("Taxa")
  }) %>% reformat_as_df("Rank")
  # correct p-values 
  out$p_value_corr <- p.adjust(out$p_value, method = "fdr")
  out
}) %>% reformat_as_df("Host")


# round values -----


for (i in c("stat", "p_value", "p_value_corr")) {
  df_kw[,i] <- round(df_kw[,i], 4)
}


## --------------------- Check the results -------------------------------------


# Percentage of significantly changed taxa ?  very high ----

df_kw %>% group_by(Host, Rank) %>% 
  summarize(n = n(),
            n_signif = sum(p_value_corr < 0.05),
            pct_signif = round((n_signif / n)*100,0))

#   Host              Rank       n n_signif pct_signif
# 1 Siganus_luridus   Class     11        2         18
# 2 Siganus_luridus   Family    19        7         37
# 3 Siganus_luridus   Genus     29        3         10
# 4 Siganus_luridus   Order     12        2         17
# 5 Siganus_luridus   Phylum     9        2         22
# 6 Siganus_rivulatus Class     11       11        100
# 7 Siganus_rivulatus Family    19       15         79
# 8 Siganus_rivulatus Genus     29       24         83
# 9 Siganus_rivulatus Order     12       12        100
# 10 Siganus_rivulatus Phylum     9        8         89



# Add the full taxonomy to each taxa -----

nms_rks <- nms_ranks_5

ls <- split(df_kw, df_kw$Rank)

ls_taxo <- lapply(nms_ranks_5, function(rk) {
  taxo <- tax_table(ls_ps_core$Siganidae[[rk]]) %>% data.frame()
  idx <- grep(rk, names(taxo))
  taxo <- taxo[,c(1:idx)]
  taxo$Taxa <- taxo[,rk]
  
  X <- ls[[rk]]
  
  out <- left_join(taxo, X, by = "Taxa")
  out
}) %>% setNames(nms_ranks_5)

df <- bind_rows(ls_taxo) %>% 
  filter(Host != "Siganidae") %>% 
  dplyr::select(all_of(nms_ranks_5), Taxa, Rank, Host, everything(), - Domain)


# save a version of the table for the supplementary ----

# The table summarizing the ASV results for each genus/family/phylum

out <- df %>% filter(Rank == "ASV") %>% 
  group_by(Host, Phylum, Family, Genus) %>% 
  summarize(n = n(),
            n_signif = sum(p_value_corr < 0.05),
            pct_signif = round((n_signif / n)*100,0),
            avg_stat = round(mean(stat),1),
            avg_p = round(mean(p_value),3)) %>% data.frame() %>% 
  arrange(Host, Phylum, Family, dplyr::desc(pct_signif))

write.csv(out, file = paste0(dir_plot, "table_summary_KW_ASV.csv"), row.names= FALSE)

# table summarizing phylum/family and genus tests

out <- df %>% filter(Rank != "ASV") %>% 
  arrange(Host, Rank, Phylum, Class, Order, Family, Genus)

write.csv(out, file = paste0(dir_plot, "table_summary_KW_phylum_family_genus.csv"), row.names= FALSE)


# plot  for the paper ----

cols_reg  <- c("#FF6347","#63B8FF", "#CDCD00") %>% setNames(nms_region)


jpeg(paste0(dir_plot, "figure_differential_abundance.jpeg"),
     height = 18, width = 20, unit = "cm", res = 600)

par(mfrow = c(4,3), mar = c(2,2,2,1), las = 1, mgp = c(3,0.75,0),
    oma = c(3,3,1,1), tcl = -0.3)

# families   ----

nms_taxa <- c("Erysipelotrichaceae", "Lachnospiraceae", "Ruminococcaceae", 
              "Akkermansiaceae", "Arcobacteraceae", "Deferribacteraceae")
dat <- ls_dat$Siganus_rivulatus$Family
otu <- dat$otu
fact <- dat$fact

lapply(nms_taxa, function(x) {
  ylim <- setrge(otu[x,], 5)
  boxplot(otu[x,] ~ fact, outline = FALSE,main = paste0(x, " (F)"),# ylim = ylim,
          names = c("NRS","LS","NC"), border = paste0(cols_reg,"80"), ylab = "",
          col = "white", lwd = 2)
  sapply(1:length(levels(fact)), function(l) {
    f <- levels(fact)[l]
    mask <- fact == f
    mask[is.na(mask)] <- FALSE
    xval <- rep(l, sum(mask)) + sample(seq(-0.25,0.25, 0.001), sum(mask))
    points(otu[x, mask] ~ xval, pch = 21, cex = 0.6,
           col = cols_reg[f],
           bg = cols_reg[f])
  })
})

# genera   ----

nms_taxa <- c("Breznakia", "Coprobacillus", "Ruminococcaceae_UCG-014",
              "Alistipes", "Labilibacter", "Desulfovibrio")
dat <- ls_dat$Siganus_rivulatus$Genus
otu <- dat$otu
fact <- dat$fact

lapply(nms_taxa, function(x) {
  ylim <- setrge(otu[x,], 5)
  boxplot(otu[x,] ~ fact, outline = FALSE,main = gsub("_", "\n", paste0(x, " (G)")), 
          # ylim = ylim,
          names = c("NRS","LS","NC"), 
          border = paste0(cols_reg,"80"), ylab = "",
          col = "white", lwd = 2)
          
  sapply(1:length(levels(fact)), function(l) {
    f <- levels(fact)[l]
    mask <- fact == f
    mask[is.na(mask)] <- FALSE
    xval <- rep(l, sum(mask)) + sample(seq(-0.2,0.2, 0.001), sum(mask))
    points(otu[x, mask] ~ xval, pch = 21, cex = 0.6,
           col = cols_reg[f],
           bg = cols_reg[f])
  })
})

# ----
legend(-5.5,2, legend = paste0(gsub("_", " ", nms_region), c(" (NRS)"," (LS)"," (NC)")), pch = 22, ncol = 3,
       bty = "n", col = col_reg, pt.bg = col_reg, pt.cex = 2, cex = 1.2,
       text.col = col_reg, xpd = NA)
mtext(side = 2, text = "CLR transformed abundance", font = 2,
      outer = TRUE, las = 0, line = 1)

dev.off()





################################################################################
# =================== DIFFERENTIAL ABUNDANCE OF FUNCTIONS ======================

dir_plot <- paste0(dir_save, "differential_abundance_functions/")
dir.create(dir_plot)

# prepare the data -----

nms_C_num <- names(ls_ps_scfa)[1:5]

# add a column of C number in the scfa

df_tax_tables <- lapply(ls_ps_scfa[1:5], function(X) {
  X@tax_table %>% data.frame()
}) %>% reformat_as_df(new_var_name = "C_number")


# extract data

dat <- ls_ps_scfa$all_scfa
dat@tax_table <- tax_table(as.matrix(df_tax_tables))

metadata <- meta(dat)


# transform into CLR

dat <- microbiome::transform(dat, transform = "clr")

# objects for analyse ----

ls_dat <- list()
ls_kw <- list()

for (sp in nms_siganids[1:2]) {
  
  # objects for plot and model fit ----
  
  fact_sp <- metadata$taxonomy == sp
  otu <- dat@otu_table@.Data %>% as.matrix()
  otu <- otu[! row.names(otu) == "Other", fact_sp]
  meta_tmp <- metadata[fact_sp,]
  meta_tmp$region <- factor(meta_tmp$region, levels = nms_region)
  fact <-  meta_tmp$region
  nms_taxa <- row.names(otu)
  
  # test differences between regions for each taxa ----
  

  ls_kw[[sp]] <- lapply(nms_taxa, function(x) {
    kruskal.test(otu[x,] ~ fact)
  }) %>% setNames(nms_taxa)
  
  # export the data ----
  
  nms_taxa <- sort(rowSums(otu), decreasing = TRUE) %>% names() %>% sort()
  otu <- otu[nms_taxa,]
  ls_dat[[sp]]$otu <- otu
  ls_dat[[sp]]$fact <- meta_tmp$region
}

saveRDS(ls_dat, file = paste0(dir_plot, "list_ANOVA_data.rds"))
saveRDS(ls_kw, file = paste0(dir_plot, "list_KW_tests.rds"))


# Analysis of difference of abundance between regions ----


df_kw <- lapply(ls_kw, function(x) {
  out <-  lapply(x, function(fit) {
    data.frame(stat = fit$statistic, p_value = fit$p.value)  
  }) %>% reformat_as_df("lev_6")
  # correct p-values 
  out$p_value_corr <- p.adjust(out$p_value, method = "fdr")
  out
}) %>% reformat_as_df("Host_species")

# round values

for (i in c("stat", "p_value", "p_value_corr")) {
  df_kw[,i] <- round(df_kw[,i], 4)
}

# add C number and clean the dataframe

row.names(df_kw) <- 1:nrow(df_kw)
df_kw <- df_kw %>% left_join(data.frame(dat@tax_table), "lev_6")


# summarize the results ------

out <- df_kw %>% 
  group_by(Host_species, C_number) %>% 
  summarize(avg_stat = mean(`stat`), sd_stat = sd(`stat`), 
            avg_p = mean(`p_value_corr`), sd_p = sd(`p_value_corr`),
            num_KO = n(),
            num_signif = sum(`p_value_corr` < 0.05),
            prop_signif = round(num_signif / n(), 2)) %>%data.frame()

write.csv(out, file = paste0(dir_plot, "table_results_KW_summary_scfa.csv"),
          row.names = FALSE)


# final plot with selected KOs -----

sp <- nms_siganids[1]

png(paste0(dir_plot, "figure_differential_function_abundance.png"),
    height = 18, width = 20, unit = "cm", res = 600)

par(mfrow = c(4,3), mar = c(2,2,2,1), las = 1, mgp = c(3,0.75,0),
    oma = c(3,3,1,1), tcl = -0.3)

dat <- ls_dat[[sp]]
otu <- dat$otu
fact <- dat$fact
nms_taxa <- row.names(otu)
tmp <- max(apply(otu, 1, function(x) max(abs(x))))
ylim <- setrge(c(-tmp, tmp))

# identify the most responding taxa 
nms_taxa <- paste0("ko", c("00167", "18122", "00929", "00042", "21993", "22339",
                           "00187", "06718", "03932", "02688", "01912","01576")) 

# how many C for KOs ?
Cnum <- df_kw %>% filter(Host_species == sp & lev_6 %in% nms_taxa) %>% 
  dplyr::select(C_number, lev_6) %>% 
  column_to_rownames("lev_6")

lapply(nms_taxa, function(x) {
  ylim <- setrge(otu[x,], 10)
  boxplot(otu[x,] ~ fact, outline = FALSE, main = paste0(x, " (",Cnum[x,], ")"), 
          ylim = ylim, names = c("NRS","LS","NC"), 
          border = paste0(col_reg,"80"), ylab = "",
          col = "white", lwd = 2)
  sapply(1:length(levels(fact)), function(l) {
    f <- levels(fact)[l]
    mask <- fact == f
    mask[is.na(mask)] <- FALSE
    xval <- rep(l, sum(mask)) + sample(seq(-0.2,0.2, 0.001), sum(mask))
    points(otu[x, mask] ~ xval, pch = 21, cex = 0.6,
           col = col_reg[f],
           bg  = col_reg[f])
  })
})
# ----
legend(-4.5,-12, legend = paste0(gsub("_", " ", nms_region), c(" (NRS)"," (LS)"," (NC)")), 
       pch = 22, ncol = 3,
       bty = "n", col = col_reg, pt.bg = col_reg, pt.cex = 2, cex = 1.2,
       text.col = col_reg, xpd = NA)
mtext(side = 2, text = "CLR transformed function abundance", font = 2,
      outer = TRUE, las = 0, line = 1)
dev.off()




