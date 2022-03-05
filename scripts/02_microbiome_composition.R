################################################################################
#  _______  _____  _____ ___ ____  _   _ __  __ _____ ____  
# | ____\ \/ / _ \|  ___|_ _/ ___|| | | |  \/  | ____|  _ \ 
# |  _|  \  / | | | |_   | |\___ \| |_| | |\/| |  _| | | | |
# | |___ /  \ |_| |  _|  | | ___) |  _  | |  | | |___| |_| |
# |_____/_/\_\___/|_|   |___|____/|_| |_|_|  |_|_____|____/ 
#
#
################################################################################
# R script to analyse microbiome composition
# arthur.escalas@gmail.com
################################################################################


############################### LOAD DATA ######################################

#     Create a file for today' analyses ----

dir_save <- dir_composition


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




############################# COMPARE SAMPLE TYPES #############################

# ============================= PREPARE DATA ===================================

# Make a directory to store the results ----

dir_out <- paste0(dir_save, "sample_types/")
dir.create(dir_out)


# create some objects for analyses and remove unecessary samples ----

metadata_tmp <- meta(ls_ps_full$Phylum)
fact <- factor(metadata_tmp$sample_type)

ps_data_tmp <- ls_ps_full$ASV

ls_ps_data_tmp <- ls_ps_full

# list of dissimilarities ----

ls_diss_tmp <- lapply(ls_diss_full, function(rk) { 
  lapply(rk$beta_taxo, function(d) {
    as.dist(d)
  })
})


# ============================= COMPOSITION =================================

dir_composition_tables <- paste0(dir_out, "composition/")
dir.create(dir_composition_tables)

# Transform data : merge samples and tidyfy ----

ls_data_plot <- lapply(nms_ranks_5, function(rk) {
  merge_samples(ls_ps_data_tmp[[rk]], "sample_type", sum) %>% 
    microbiome::transform(transform = "compositional")
}) %>% setNames(nms_ranks_5)

# Make and save the composition tables -----

  lapply(nms_ranks_5, function(rk) {
    
    # Save the table of composition ----
    dat <- ls_data_plot[[rk]]
    otu <- dat %>% otu_table() %>% as.matrix() %>% t()
    tax <- dat %>% tax_table()
    out <- cbind(tax, round(otu * 100, 2))
    write.csv(cbind(tax, otu), row.names = FALSE,
              file = paste0(dir_composition_tables, "table_composition_", rk, ".csv"))
    
    # Make the plot with only the top 20 taxa ----
    dat <- aggregate_top_taxa(dat, top = 20, level = rk)
    otu <- dat %>% otu_table() %>% as.matrix()
    tax <- dat %>% tax_table()
    rownames(otu) <- tax[, rk]
    out <- cbind(tax, round(otu * 100, 2))
    write.csv(out, row.names = FALSE,
              file = paste0(dir_composition_tables, "table_composition_top_20_", rk, ".csv"))
  })


#============================= PCOA ORDINATION =================================

dir_pcoa <- paste0(dir_out, "pcoa/")
dir.create(dir_pcoa)

# Make the PCOA ----

ls_pcoa <- lapply(ls_diss_tmp, function(rk) { 
  lapply(rk, function(d) pcoa(D = as.dist(d)) )
})

saveRDS(ls_pcoa, file = paste0(dir_pcoa, "list_pcoa.rds"))


# =============== TEST DIFFERENCE BETWEEN GROUPS ===============================

dir_stats <- paste0(dir_out, "stats/")
dir.create(dir_stats)

# Run the tests ----

ls_res <- lapply(nms_ranks_3, function(rk) {
  lapply(ls_diss_tmp[[rk]], function(d) {
    tmp <- adonis2(d ~ fact, data = metadata_tmp, nrep = 999, by = "margin")
    out <- data.frame(tmp) %>% round(3) %>% 
      setNames(c("df","Sum_sq","R2","F-value","p-value")) %>% 
      rownames_to_column("factor")    
  }) %>% reformat_as_df(new_var_name = "data_type")
}) %>% setNames(nms_ranks_3)

saveRDS(ls_res, file = paste0(dir_stats, "res_permanova.rds"))
ls_res <- readRDS(paste0(dir_stats, "res_permanova.rds"))


#     format the results as a table ----

df_res <- reformat_as_df(ls_res, new_var_name = "rank")
df_res$q <- factor(str_split_fixed(df_res$data_type, "_", 2)[,2])
df_res$div_facet <- factor(str_split_fixed(df_res$data_type, "_", 2)[,1])

write.csv(df_res, file = paste0(dir_stats, "table_results_permanova.csv"), row.names =FALSE)




########################### COMPARE WATER SAMPLES ##############################

# ============================= PREPARE DATA ===================================

# Make a directory to store the results ----

dir_out <- paste0(dir_save, "water/")
dir.create(dir_out)

# create some objects for analyses and remove unecessary samples ----

metadata_tmp <- meta(ls_ps_full$Phylum)
metadata_tmp$region <- factor(metadata_tmp$region, levels = nms_region)

mask <- metadata_tmp$sample_type == "Water"
metadata_tmp <- metadata_tmp[mask,]%>% data.frame()

fact <- factor(metadata_tmp$sample_type)

ps_data_tmp <- subset_samples(ls_ps_full$ASV, mask) 

ls_ps_data_tmp <- lapply(ls_ps_full, function(X) {
  subset_samples(X, mask)
})

# list of dissimilarities ----

ls_diss_tmp <- lapply(ls_diss_full, function(rk) { 
  lapply(rk$beta_taxo, function(d) {
    as.dist(d[mask, mask])
  })
})

# ============================= COMPOSITION =================================

dir_composition_tables <- paste0(dir_out, "composition/")
dir.create(dir_composition_tables)

# Transform data : merge samples and tidyfy ----

ls_data_plot <- lapply(nms_ranks_5, function(rk) {
  merge_samples(ls_ps_data_tmp[[rk]], "region_season", sum) %>% 
    microbiome::transform(transform = "compositional")
}) %>% setNames(nms_ranks_5)

# Make and save the composition tables -----

lapply(nms_ranks_5, function(rk) {
  
  # Save the table of composition ----
  dat <- ls_data_plot[[rk]]
  otu <- dat %>% otu_table() %>% as.matrix() %>% t()
  tax <- dat %>% tax_table()
  out <- cbind(tax, round(otu * 100, 2))
  write.csv(cbind(tax, otu), row.names = FALSE,
            file = paste0(dir_composition_tables, "table_composition_", rk, ".csv"))
  
  # Make the plot with only the top 20 taxa ----
  dat <- aggregate_top_taxa(dat, top = 20, level = rk)
  otu <- dat %>% otu_table() %>% as.matrix()
  tax <- dat %>% tax_table()
  rownames(otu) <- tax[, rk]
  out <- cbind(tax, round(otu * 100, 2))
  write.csv(out, row.names = FALSE,
            file = paste0(dir_composition_tables, "table_composition_top_20_", rk, ".csv"))
})


#============================= PCOA ORDINATION =================================

dir_pcoa <- paste0(dir_out, "pcoa/")
dir.create(dir_pcoa)

# Make the PCOA ----

ls_pcoa <- lapply(ls_diss_tmp, function(rk) { 
  lapply(rk, function(d) pcoa(D = as.dist(d)) )
})

saveRDS(ls_pcoa, file = paste0(dir_pcoa, "list_pcoa.rds"))


# =============== TEST DIFFERENCE BETWEEN GROUPS ===============================

dir_stats <- paste0(dir_out, "stats/")
dir.create(dir_stats)

# Run the tests ----

ls_res <- lapply(nms_ranks_3, function(rk) {
  lapply(ls_diss_tmp[[rk]], function(d) {
    tmp_main <- adonis2(d ~ region + season, data = metadata_tmp, nrep = 999, by = "margin")
    tmp_inte <- adonis2(d ~ region * season, data = metadata_tmp, nrep = 999, by = "margin")
    df_main <- data.frame(tmp_main)
    df_main <- df_main[! row.names(df_main) %in% c("Residual","Total"),]
    df_inte <- data.frame(tmp_inte)
    out <- data.frame(rbind(df_main,df_inte)) %>% round(3) %>% 
      setNames(c("df","Sum_sq","R2","F-value","p-value")) %>% 
      rownames_to_column("factor")    
  }) %>% reformat_as_df(new_var_name = "data_type")
}) %>% setNames(nms_ranks_3)

saveRDS(ls_res, file = paste0(dir_stats, "res_permanova.rds"))
ls_res <- readRDS(paste0(dir_stats, "res_permanova.rds"))


#     format the results as a table ----

df_res <- reformat_as_df(ls_res, new_var_name = "rank")
df_res$q <- factor(str_split_fixed(df_res$data_type, "_", 2)[,2])
df_res$div_facet <- factor(str_split_fixed(df_res$data_type, "_", 2)[,1])

write.csv(df_res, file = paste0(dir_stats, "table_results_permanova.csv"), row.names =FALSE)


# =============== TEST DISPERSION BETWEEN GROUPS ===============================

dir_stats <- paste0(dir_out, "stats/")
dir.create(dir_stats)

# Run the tests ----

ls_res <- lapply(nms_ranks_3, function(rk) {
  lapply(ls_diss_tmp[[rk]], function(d) {
    tmp <- betadisper(d, metadata_tmp$region, type = "centroid", bias.adjust = FALSE,
               sqrt.dist = FALSE, add = FALSE)
    return(tmp)
  }) %>% setNames(c("q0","q1"))
}) %>% setNames(nms_ranks_3)

saveRDS(ls_res, file = paste0(dir_stats, "res_permdisp.rds"))


## make the results as pretty tables

ls_res_df <- lapply(ls_res, function(X) {
  lapply(X, function(x) {
    anova(x) %>% data.frame() %>% round(3) %>% 
      setNames(c("df","Sum_sq","Mean_sq","F-value","p-value")) %>% 
      rownames_to_column("factor")    
  }) %>% reformat_as_df(new_var_name = "data_type")
}) %>% setNames(nms_ranks_3)


#     format the results as a table ----

df_res <- reformat_as_df(ls_res_df, new_var_name = "rank")
df_res$q <- factor(str_split_fixed(df_res$data_type, "_", 2)[,2])
df_res$div_facet <- factor(str_split_fixed(df_res$data_type, "_", 2)[,1])

write.csv(df_res, file = paste0(dir_stats, "table_results_permdisp.csv"), row.names =FALSE)





########################### COMPARE TURF SAMPLES ##############################


# ============================= PREPARE DATA ===================================

# Make a directory to store the results ----

dir_out <- paste0(dir_save, "turf/")
dir.create(dir_out)

# create some objects for analyses and remove unecessary samples ----

metadata_tmp <- meta(ls_ps_full$Phylum)
metadata_tmp$region <- factor(metadata_tmp$region, levels = nms_region)

mask <- metadata_tmp$sample_type == "Turf"
metadata_tmp <- metadata_tmp[mask,]%>% data.frame()

fact <- factor(metadata_tmp$sample_type)

ps_data_tmp <- subset_samples(ls_ps_full$ASV, mask) 

ls_ps_data_tmp <- lapply(ls_ps_full, function(X) {
  subset_samples(X, mask)
})

# list of dissimilarities ----

ls_diss_tmp <- lapply(ls_diss_full, function(rk) { 
  lapply(rk$beta_taxo, function(d) {
    as.dist(d[mask, mask])
  })
})


# ============================= COMPOSITION =================================

dir_composition_tables <- paste0(dir_out, "composition/")
dir.create(dir_composition_tables)

# Transform data : merge samples and tidyfy ----

ls_data_plot <- lapply(nms_ranks_5, function(rk) {
  merge_samples(ls_ps_data_tmp[[rk]], "region_season", sum) %>% 
    microbiome::transform(transform = "compositional")
}) %>% setNames(nms_ranks_5)


# Make and save the composition tables -----

lapply(nms_ranks_5, function(rk) {
  
  # Save the table of composition ----
  dat <- ls_data_plot[[rk]]
  otu <- dat %>% otu_table() %>% as.matrix() %>% t()
  tax <- dat %>% tax_table()
  out <- cbind(tax, round(otu * 100, 2))
  write.csv(cbind(tax, otu), row.names = FALSE,
            file = paste0(dir_composition_tables, "table_composition_", rk, ".csv"))
  
  # Make the plot with only the top 20 taxa ----
  dat <- aggregate_top_taxa(dat, top = 20, level = rk)
  otu <- dat %>% otu_table() %>% as.matrix()
  tax <- dat %>% tax_table()
  rownames(otu) <- tax[, rk]
  out <- cbind(tax, round(otu * 100, 2))
  write.csv(out, row.names = FALSE,
            file = paste0(dir_composition_tables, "table_composition_top_20_", rk, ".csv"))
})


#============================= PCOA ORDINATION =================================

dir_pcoa <- paste0(dir_out, "pcoa/")
dir.create(dir_pcoa)

# Make the PCOA ----

ls_pcoa <- lapply(ls_diss_tmp, function(rk) { 
  lapply(rk, function(d) pcoa(D = as.dist(d)) )
})

saveRDS(ls_pcoa, file = paste0(dir_pcoa, "list_pcoa.rds"))


# =============== TEST DIFFERENCE BETWEEN GROUPS ===============================

dir_stats <- paste0(dir_out, "stats/")
dir.create(dir_stats)

# Run the tests ----

ls_res <- lapply(nms_ranks_3, function(rk) {
  lapply(ls_diss_tmp[[rk]], function(d) {
    tmp_main <- adonis2(d ~ region + season, data = metadata_tmp, nrep = 999, by = "margin")
    tmp_inte <- adonis2(d ~ region * season, data = metadata_tmp, nrep = 999, by = "margin")
    df_main <- data.frame(tmp_main)
    df_main <- df_main[! row.names(df_main) %in% c("Residual","Total"),]
    df_inte <- data.frame(tmp_inte)
    out <- data.frame(rbind(df_main,df_inte)) %>% round(3) %>% 
      setNames(c("df","Sum_sq","R2","F-value","p-value")) %>% 
      rownames_to_column("factor")   
  }) %>% reformat_as_df(new_var_name = "data_type")
}) %>% setNames(nms_ranks_3)

saveRDS(ls_res, file = paste0(dir_stats, "res_permanova.rds"))
ls_res <- readRDS(paste0(dir_stats, "res_permanova.rds"))


#     format the results as a table ----

df_res <- reformat_as_df(ls_res, new_var_name = "rank")
df_res$q <- factor(str_split_fixed(df_res$data_type, "_", 2)[,2])
df_res$div_facet <- factor(str_split_fixed(df_res$data_type, "_", 2)[,1])

write.csv(df_res, file = paste0(dir_stats, "table_results_permanova.csv"), row.names =FALSE)


# =============== TEST DISPERSION BETWEEN GROUPS ===============================

dir_stats <- paste0(dir_out, "stats/")
dir.create(dir_stats)

# Run the tests ----

ls_res <- lapply(nms_ranks_3, function(rk) {
  lapply(ls_diss_tmp[[rk]], function(d) {
    tmp <- betadisper(d, metadata_tmp$region, type = "centroid", bias.adjust = FALSE,
                      sqrt.dist = FALSE, add = FALSE)
    return(tmp)
  }) %>% setNames(c("q0","q1"))
}) %>% setNames(nms_ranks_3)

saveRDS(ls_res, file = paste0(dir_stats, "res_permdisp.rds"))


## make the results as pretty tables

ls_res_df <- lapply(ls_res, function(X) {
  lapply(X, function(x) {
    anova(x) %>% data.frame() %>% round(3) %>% 
      setNames(c("df","Sum_sq","Mean_sq","F-value","p-value")) %>% 
      rownames_to_column("factor")    
  }) %>% reformat_as_df(new_var_name = "data_type")
}) %>% setNames(nms_ranks_3)


#     format the results as a table ----

df_res <- reformat_as_df(ls_res_df, new_var_name = "rank")
df_res$q <- factor(str_split_fixed(df_res$data_type, "_", 2)[,2])
df_res$div_facet <- factor(str_split_fixed(df_res$data_type, "_", 2)[,1])

write.csv(df_res, file = paste0(dir_stats, "table_results_permdisp.csv"), row.names =FALSE)





########################### COMPARE ALGAE SAMPLES ##############################


# ============================= PREPARE DATA ===================================

# Make a directory to store the results ----

dir_out <- paste0(dir_save, "algae/")
dir.create(dir_out)


# create some objects for analyses and remove unecessary samples ----

metadata_tmp <- meta(ls_ps_full$Phylum)
metadata_tmp$region <- factor(metadata_tmp$region, levels = nms_region)

mask <- metadata_tmp$sample_type == "Algae"
metadata_tmp <- metadata_tmp[mask,]%>% data.frame()

fact <- factor(metadata_tmp$sample_type)

ps_data_tmp <- subset_samples(ls_ps_full$ASV, mask) 

ls_ps_data_tmp <- lapply(ls_ps_full, function(X) {
  subset_samples(X, mask)
})

# list of dissimilarities ----

ls_diss_tmp <- lapply(ls_diss_full, function(rk) { 
  lapply(rk$beta_taxo, function(d) {
    as.dist(d[mask, mask])
  })
})



# ============================= COMPOSITION =================================

dir_composition_tables <- paste0(dir_out, "composition/")
dir.create(dir_composition_tables)

# Transform data : merge samples and tidyfy ----

ls_data_plot <- lapply(nms_ranks_5, function(rk) {
  merge_samples(ls_ps_data_tmp[[rk]], "region_season", sum) %>% 
    microbiome::transform(transform = "compositional")
}) %>% setNames(nms_ranks_5)


# Make and save the composition tables -----

lapply(nms_ranks_5, function(rk) {
  
  # Save the table of composition ----
  dat <- ls_data_plot[[rk]]
  otu <- dat %>% otu_table() %>% as.matrix() %>% t()
  tax <- dat %>% tax_table()
  out <- cbind(tax, round(otu * 100, 2))
  write.csv(cbind(tax, otu), row.names = FALSE,
            file = paste0(dir_composition_tables, "table_composition_", rk, ".csv"))
  
  # Make the plot with only the top 20 taxa ----
  dat <- aggregate_top_taxa(dat, top = 20, level = rk)
  otu <- dat %>% otu_table() %>% as.matrix()
  tax <- dat %>% tax_table()
  rownames(otu) <- tax[, rk]
  out <- cbind(tax, round(otu * 100, 2))
  write.csv(out, row.names = FALSE,
            file = paste0(dir_composition_tables, "table_composition_top_20_", rk, ".csv"))
})


#============================= PCOA ORDINATION =================================

dir_pcoa <- paste0(dir_out, "pcoa/")
dir.create(dir_pcoa)

# Make the PCOA ----

ls_pcoa <- lapply(ls_diss_tmp, function(rk) { 
  lapply(rk, function(d) pcoa(D = as.dist(d)) )
})

saveRDS(ls_pcoa, file = paste0(dir_pcoa, "list_pcoa.rds"))


# =============== TEST DIFFERENCE BETWEEN GROUPS ===============================

dir_stats <- paste0(dir_out, "stats/")
dir.create(dir_stats)

# Run the tests ----

ls_res <- lapply(nms_ranks_3, function(rk) {
  lapply(ls_diss_tmp[[rk]], function(d) {
    tmp_main <- adonis2(d ~ region + season, data = metadata_tmp, nrep = 999, by = "margin")
    tmp_inte <- adonis2(d ~ region * season, data = metadata_tmp, nrep = 999, by = "margin")
    df_main <- data.frame(tmp_main)
    df_main <- df_main[! row.names(df_main) %in% c("Residual","Total"),]
    df_inte <- data.frame(tmp_inte)
    out <- data.frame(rbind(df_main,df_inte)) %>% round(3) %>% 
      setNames(c("df","Sum_sq","R2","F-value","p-value")) %>% 
      rownames_to_column("factor")   
  }) %>% reformat_as_df(new_var_name = "data_type")
}) %>% setNames(nms_ranks_3)

saveRDS(ls_res, file = paste0(dir_stats, "res_permanova.rds"))
ls_res <- readRDS(paste0(dir_stats, "res_permanova.rds"))


#     format the results as a table ----

df_res <- reformat_as_df(ls_res, new_var_name = "rank")
df_res$q <- factor(str_split_fixed(df_res$data_type, "_", 2)[,2])
df_res$div_facet <- factor(str_split_fixed(df_res$data_type, "_", 2)[,1])

write.csv(df_res, file = paste0(dir_stats, "table_results_permanova.csv"), row.names =FALSE)


# =============== TEST DISPERSION BETWEEN GROUPS ===============================

dir_stats <- paste0(dir_out, "stats/")
dir.create(dir_stats)

# Run the tests ----

ls_res <- lapply(nms_ranks_3, function(rk) {
  lapply(ls_diss_tmp[[rk]], function(d) {
    tmp <- betadisper(d, metadata_tmp$region, type = "centroid", bias.adjust = FALSE,
                      sqrt.dist = FALSE, add = FALSE)
    return(tmp)
  }) %>% setNames(c("q0","q1"))
}) %>% setNames(nms_ranks_3)

saveRDS(ls_res, file = paste0(dir_stats, "res_permdisp.rds"))


## make the results as pretty tables

ls_res_df <- lapply(ls_res, function(X) {
  lapply(X, function(x) {
    anova(x) %>% data.frame() %>% round(3) %>% 
      setNames(c("df","Sum_sq","Mean_sq","F-value","p-value")) %>% 
      rownames_to_column("factor")    
  }) %>% reformat_as_df(new_var_name = "data_type")
}) %>% setNames(nms_ranks_3)


#     format the results as a table ----

df_res <- reformat_as_df(ls_res_df, new_var_name = "rank")
df_res$q <- factor(str_split_fixed(df_res$data_type, "_", 2)[,2])
df_res$div_facet <- factor(str_split_fixed(df_res$data_type, "_", 2)[,1])

write.csv(df_res, file = paste0(dir_stats, "table_results_permdisp.csv"), row.names =FALSE)



########################### COMPARE SEDIMENT SAMPLES ##############################


# ============================= PREPARE DATA ===================================

# Make a directory to store the results ----

dir_out <- paste0(dir_save, "sediment/")
dir.create(dir_out)


# create some objects for analyses and remove unecessary samples ----

metadata_tmp <- meta(ls_ps_full$Phylum)
metadata_tmp$region <- factor(metadata_tmp$region, levels = nms_region)

mask <- metadata_tmp$sample_type == "Sediment"
metadata_tmp <- metadata_tmp[mask,]%>% data.frame()

fact <- factor(metadata_tmp$region)

ps_data_tmp <- subset_samples(ls_ps_full$ASV, mask) 

ls_ps_data_tmp <- lapply(ls_ps_full, function(X) {
  subset_samples(X, mask)
})

# list of dissimilarities ----

ls_diss_tmp <- lapply(ls_diss_full, function(rk) { 
  lapply(rk$beta_taxo, function(d) {
    as.dist(d[mask, mask])
  })
})


# ============================= COMPOSITION =================================

dir_composition_tables <- paste0(dir_out, "composition/")
dir.create(dir_composition_tables)

# Transform data : merge samples and tidyfy ----

ls_data_plot <- lapply(nms_ranks_5, function(rk) {
  merge_samples(ls_ps_data_tmp[[rk]], "region_season", sum) %>% 
    microbiome::transform(transform = "compositional")
}) %>% setNames(nms_ranks_5)

# Make and save the composition tables -----

lapply(nms_ranks_5, function(rk) {
  
  # Save the table of composition ----
  dat <- ls_data_plot[[rk]]
  otu <- dat %>% otu_table() %>% as.matrix() %>% t()
  tax <- dat %>% tax_table()
  out <- cbind(tax, round(otu * 100, 2))
  write.csv(cbind(tax, otu), row.names = FALSE,
            file = paste0(dir_composition_tables, "table_composition_", rk, ".csv"))
  
  # Make the plot with only the top 20 taxa ----
  dat <- aggregate_top_taxa(dat, top = 20, level = rk)
  otu <- dat %>% otu_table() %>% as.matrix()
  tax <- dat %>% tax_table()
  rownames(otu) <- tax[, rk]
  out <- cbind(tax, round(otu * 100, 2))
  write.csv(out, row.names = FALSE,
            file = paste0(dir_composition_tables, "table_composition_top_20_", rk, ".csv"))
})


#============================= PCOA ORDINATION =================================

dir_pcoa <- paste0(dir_out, "pcoa/")
dir.create(dir_pcoa)

# Make the PCOA ----

ls_pcoa <- lapply(ls_diss_tmp, function(rk) { 
  lapply(rk, function(d) pcoa(D = as.dist(d)) )
})

saveRDS(ls_pcoa, file = paste0(dir_pcoa, "list_pcoa.rds"))


# =============== TEST DIFFERENCE BETWEEN GROUPS ===============================

dir_stats <- paste0(dir_out, "stats/")
dir.create(dir_stats)

# Run the tests ----

ls_res <- lapply(nms_ranks_3, function(rk) {
  lapply(ls_diss_tmp[[rk]], function(d) {
    tmp_main <- adonis2(d ~ region + season, data = metadata_tmp, nrep = 999, by = "margin")
    tmp_inte <- adonis2(d ~ region * season, data = metadata_tmp, nrep = 999, by = "margin")
    df_main <- data.frame(tmp_main)
    df_main <- df_main[! row.names(df_main) %in% c("Residual","Total"),]
    df_inte <- data.frame(tmp_inte)
    out <- data.frame(rbind(df_main,df_inte)) %>% round(3) %>% 
      setNames(c("df","Sum_sq","R2","F-value","p-value")) %>% 
      rownames_to_column("factor")   
  }) %>% reformat_as_df(new_var_name = "data_type")
}) %>% setNames(nms_ranks_3)

saveRDS(ls_res, file = paste0(dir_stats, "res_permanova.rds"))
ls_res <- readRDS(paste0(dir_stats, "res_permanova.rds"))


#     format the results as a table ----

df_res <- reformat_as_df(ls_res, new_var_name = "rank")
df_res$q <- factor(str_split_fixed(df_res$data_type, "_", 2)[,2])
df_res$div_facet <- factor(str_split_fixed(df_res$data_type, "_", 2)[,1])

write.csv(df_res, file = paste0(dir_stats, "table_results_permanova.csv"), row.names =FALSE)


# =============== TEST DISPERSION BETWEEN GROUPS ===============================

dir_stats <- paste0(dir_out, "stats/")
dir.create(dir_stats)

# Run the tests ----

ls_res <- lapply(nms_ranks_3, function(rk) {
  lapply(ls_diss_tmp[[rk]], function(d) {
    tmp <- betadisper(d, metadata_tmp$region, type = "centroid", bias.adjust = FALSE,
                      sqrt.dist = FALSE, add = FALSE)
    return(tmp)
  }) %>% setNames(c("q0","q1"))
}) %>% setNames(nms_ranks_3)

saveRDS(ls_res, file = paste0(dir_stats, "res_permdisp.rds"))


## make the results as pretty tables

ls_res_df <- lapply(ls_res, function(X) {
  lapply(X, function(x) {
    anova(x) %>% data.frame() %>% round(3) %>% 
      setNames(c("df","Sum_sq","Mean_sq","F-value","p-value")) %>% 
      rownames_to_column("factor")    
  }) %>% reformat_as_df(new_var_name = "data_type")
}) %>% setNames(nms_ranks_3)


#     format the results as a table ----

df_res <- reformat_as_df(ls_res_df, new_var_name = "rank")
df_res$q <- factor(str_split_fixed(df_res$data_type, "_", 2)[,2])
df_res$div_facet <- factor(str_split_fixed(df_res$data_type, "_", 2)[,1])

write.csv(df_res, file = paste0(dir_stats, "table_results_permdisp.csv"), row.names =FALSE)






######################## COMPARE RIVULATUS SAMPLES #############################


# ============================= PREPARE DATA ===================================

# Make a directory to store the results ----

dir_out <- paste0(dir_save, "s_rivulatus/")
dir.create(dir_out)


# create some objects for analyses and remove unecessary samples ----

metadata_tmp <- meta(ls_ps_core$Siganus_rivulatus$Phylum)
metadata_tmp$region <- factor(metadata_tmp$region, levels = nms_region)

mask <- metadata_tmp$sample_type_2 == "Siganus_rivulatus"
metadata_tmp <- metadata_tmp[mask,]%>% data.frame()

fact <- factor(metadata_tmp$sample_type)

ps_data_tmp <- subset_samples(ls_ps_core$Siganus_rivulatus$ASV, mask) 

ls_ps_data_tmp <- lapply(ls_ps_core$Siganus_rivulatus, function(X) {
  subset_samples(X, mask)
})

# list of dissimilarities ----

ls_diss_tmp <- lapply(ls_diss_fish$Siganus_rivulatus, function(rk) { 
  out <- lapply(flatten(rk[c("beta_taxo", "beta_phylo")]), function(d) {
    as.dist(d[mask, mask])
  }) %>% setNames(c("taxo_q0","taxo_q1","phylo_q0","phylo_q1"))
  out
})


# ============================= COMPOSITION =================================

dir_composition_tables <- paste0(dir_out, "composition/")
dir.create(dir_composition_tables)

# Transform data : merge samples and tidyfy ----

ls_data_plot <- lapply(nms_ranks_5, function(rk) {
  merge_samples(ls_ps_data_tmp[[rk]], "region_season", sum) %>% 
    microbiome::transform(transform = "compositional")
}) %>% setNames(nms_ranks_5)

# Make and save the composition tables -----

lapply(nms_ranks_5, function(rk) {
  
  # Save the table of composition ----
  dat <- ls_data_plot[[rk]]
  otu <- dat %>% otu_table() %>% as.matrix() %>% t()
  tax <- dat %>% tax_table()
  out <- cbind(tax, round(otu * 100, 2))
  write.csv(out, row.names = FALSE,
            file = paste0(dir_composition_tables, "table_composition_", rk, ".csv"))
  out2 <- merge_samples(ls_ps_data_tmp[[rk]], "region", sum) %>% 
    microbiome::transform(transform = "compositional") 
  out2 <- cbind(dat %>% tax_table(), round(out2 %>% otu_table() %>% 
                                            as.matrix() %>% t() * 100, 2))
  write.csv(out2, row.names = FALSE,
            file = paste0(dir_composition_tables, "table_composition_per_region_", rk, ".csv"))
  
  # Make the plot with only the top 20 taxa ----
  dat <- aggregate_top_taxa(dat, top = 20, level = rk)
  otu <- dat %>% otu_table() %>% as.matrix()
  tax <- dat %>% tax_table()
  rownames(otu) <- tax[, rk]
  out <- cbind(tax, round(otu * 100, 2))
  write.csv(out, row.names = FALSE,
            file = paste0(dir_composition_tables, "table_composition_top_20_", rk, ".csv"))
})



#============================= PCOA ORDINATION =================================

dir_pcoa <- paste0(dir_out, "pcoa/")
dir.create(dir_pcoa)

# Make the PCOA ----

ls_pcoa <- lapply(ls_diss_tmp, function(rk) { 
  lapply(rk, function(d) pcoa(D = as.dist(d)) )
})

saveRDS(ls_pcoa, file = paste0(dir_pcoa, "list_pcoa.rds"))



# =============== TEST DIFFERENCE BETWEEN GROUPS ===============================

dir_stats <- paste0(dir_out, "stats/")
dir.create(dir_stats)

# Run the tests ----

ls_res <- lapply(nms_ranks_3, function(rk) {
  lapply(ls_diss_tmp[[rk]], function(d) {
    tmp_reg <- adonis2(d ~ region, data = metadata_tmp, nrep = 999, by = "margin")
    tmp_sea <- adonis2(d ~ season, data = metadata_tmp, nrep = 999, by = "margin")
    tmp_inte <- adonis2(d ~ region * season, data = metadata_tmp, nrep = 999, by = "margin")
    df_reg <- data.frame(tmp_reg)
    df_sea <- data.frame(tmp_sea)
    df_reg <- df_reg[! row.names(df_reg) %in% c("Residual","Total"),]
    df_sea <- df_sea[! row.names(df_sea) %in% c("Residual","Total"),]
    df_inte <- data.frame(tmp_inte)
    out <- data.frame(rbind(df_reg, df_sea ,df_inte)) %>% round(3) %>% 
      setNames(c("df","Sum_sq","R2","F-value","p-value")) %>% 
      rownames_to_column("factor")   
  }) %>% reformat_as_df(new_var_name = "data_type")
}) %>% setNames(nms_ranks_3)

saveRDS(ls_res, file = paste0(dir_stats, "res_permanova.rds"))
ls_res <- readRDS(paste0(dir_stats, "res_permanova.rds"))


#     format the results as a table ----

df_res <- reformat_as_df(ls_res, new_var_name = "rank")
df_res$q <- factor(str_split_fixed(df_res$data_type, "_", 2)[,2])
df_res$div_facet <- factor(str_split_fixed(df_res$data_type, "_", 2)[,1])

write.csv(df_res, file = paste0(dir_stats, "table_results_permanova.csv"), row.names =FALSE)


# Pairwise comparisons for revision animal µbiome -------

ls_res <- list()

for (rk in nms_ranks_3) {
  for (idx in names(ls_diss_tmp[[rk]])) {
    d <- ls_diss_tmp[[rk]][[idx]]
    res <- pairwise.adonis2(d ~ region + season, data = metadata_tmp, nrep = 999)
    out <- lapply(res[-1], function(x) data.frame(x) %>% rownames_to_column("factor")) %>% 
      reformat_as_df(new_var_name = "comparison") %>% 
      setNames(c("Factor","DF","Sum_sq","R2","F_value","p_value", "Comparison")) %>% 
      filter(! Factor %in% c("Residual", "Total"))
    ls_res[[rk]][[idx]] <- out
  }
}

df_res_pairwise <- lapply(ls_res, function(x) { 
    x %>% reformat_as_df(new_var_name = "Index")
}) %>% reformat_as_df(new_var_name = "Rank")

df_res_out <- df_res_pairwise %>% arrange(Index, Rank) %>% 
  dplyr::select(Index, Rank, Comparison, Factor, R2, F_value, p_value, everything())


write.csv(df_res_out,
          file = paste0(dir_stats, "table_results_permanova_pairwise.csv"),
          row.names =FALSE)

out <- df_res_out %>% group_by(Comparison, Factor) %>% 
  summarize(num_signif = sum(p_value < 0.05), 
            avg_F = mean(F_value), sd_F = sd(F_value), 
            avg_p = mean(p_value), sd_p = sd(p_value),
            avg_R2 = mean(R2), sd_R2 = sd(R2)) %>% data.frame()

write.csv(out,
          file = paste0(dir_stats, "table_results_permanova_pairwise_summary.csv"),
          row.names =FALSE)


# =============== TEST DISPERSION BETWEEN GROUPS ===============================

dir_stats <- paste0(dir_out, "stats/")
dir.create(dir_stats)

# Run the tests ----

ls_res <- lapply(nms_ranks_3, function(rk) {
  lapply(ls_diss_tmp[[rk]], function(d) {
    tmp <- betadisper(d, metadata_tmp$region, type = "centroid", bias.adjust = FALSE,
                      sqrt.dist = FALSE, add = FALSE)
    return(tmp)
  }) %>% setNames(names(ls_diss_tmp$Phylum))
}) %>% setNames(nms_ranks_3)

saveRDS(ls_res, file = paste0(dir_stats, "res_permdisp.rds"))


## make the results as pretty tables

ls_res_df <- lapply(ls_res, function(X) {
  lapply(X, function(x) {
    anova(x) %>% data.frame() %>% round(3) %>% 
      setNames(c("df","Sum_sq","Mean_sq","F-value","p-value")) %>% 
      rownames_to_column("factor")    
  }) %>% reformat_as_df(new_var_name = "data_type")
}) %>% setNames(nms_ranks_3)


#     format the results as a table ----

df_res <- reformat_as_df(ls_res_df, new_var_name = "rank")
df_res$q <- factor(str_split_fixed(df_res$data_type, "_", 2)[,2])
df_res$div_facet <- factor(str_split_fixed(df_res$data_type, "_", 2)[,1])

write.csv(df_res, file = paste0(dir_stats, "table_results_permdisp.csv"), row.names =FALSE)



######################## COMPARE LURIDUS SAMPLES #############################


# ============================= PREPARE DATA ===================================

# Make a directory to store the results ----

dir_out <- paste0(dir_save, "s_luridus/")
dir.create(dir_out)


# create some objects for analyses and remove unecessary samples ----

metadata_tmp <- meta(ls_ps_core$Siganus_luridus$Phylum)
metadata_tmp$region <- factor(metadata_tmp$region, levels = nms_region)

mask <- metadata_tmp$sample_type_2 == "Siganus_luridus"
metadata_tmp <- metadata_tmp[mask,]%>% data.frame()

fact <- factor(metadata_tmp$sample_type)

ps_data_tmp <- subset_samples(ls_ps_core$Siganus_luridus$ASV, mask) 

ls_ps_data_tmp <- lapply(ls_ps_core$Siganus_luridus, function(X) {
  subset_samples(X, mask)
})


# list of dissimilarities ----

ls_diss_tmp <- lapply(ls_diss_fish$Siganus_luridus, function(rk) { 
  out <- lapply(flatten(rk[c("beta_taxo", "beta_phylo")]), function(d) {
    as.dist(d[mask, mask])
  }) %>% setNames(c("taxo_q0","taxo_q1","phylo_q0","phylo_q1"))
  out
})


# ============================= COMPOSITION =================================

dir_composition_tables <- paste0(dir_out, "composition/")
dir.create(dir_composition_tables)

# Transform data : merge samples and tidyfy ----

ls_data_plot <- lapply(nms_ranks_5, function(rk) {
  merge_samples(ls_ps_data_tmp[[rk]], "region_season", sum) %>% 
    microbiome::transform(transform = "compositional")
}) %>% setNames(nms_ranks_5)

# Make and save the composition tables -----

lapply(nms_ranks_5, function(rk) {
  
  # Save the table of composition ----
  dat <- ls_data_plot[[rk]]
  otu <- dat %>% otu_table() %>% as.matrix() %>% t()
  tax <- dat %>% tax_table()
  out <- cbind(tax, round(otu * 100, 2))
  write.csv(out, row.names = FALSE,
            file = paste0(dir_composition_tables, "table_composition_", rk, ".csv"))
  out2 <- merge_samples(ls_ps_data_tmp[[rk]], "region", sum) %>% 
    microbiome::transform(transform = "compositional") 
  out2 <- cbind(dat %>% tax_table(), round(out2 %>% otu_table() %>% 
                                             as.matrix() %>% t() * 100, 2))
  write.csv(out2, row.names = FALSE,
            file = paste0(dir_composition_tables, "table_composition_per_region_", rk, ".csv"))
  
  # Make the plot with only the top 20 taxa ----
  dat <- aggregate_top_taxa(dat, top = 20, level = rk)
  otu <- dat %>% otu_table() %>% as.matrix()
  tax <- dat %>% tax_table()
  rownames(otu) <- tax[, rk]
  out <- cbind(tax, round(otu * 100, 2))
  write.csv(out, row.names = FALSE,
            file = paste0(dir_composition_tables, "table_composition_top_20_", rk, ".csv"))
})



#============================= PCOA ORDINATION =================================

dir_pcoa <- paste0(dir_out, "pcoa/")
dir.create(dir_pcoa)

# Make the PCOA ----

ls_pcoa <- lapply(ls_diss_tmp, function(rk) { 
  lapply(rk, function(d) pcoa(D = as.dist(d)) )
})

saveRDS(ls_pcoa, file = paste0(dir_pcoa, "list_pcoa.rds"))


## =============== TEST DIFFERENCE BETWEEN GROUPS ===============================

dir_stats <- paste0(dir_out, "stats/")
dir.create(dir_stats)

# Run the tests ----

ls_res <- lapply(nms_ranks_3, function(rk) {
  lapply(ls_diss_tmp[[rk]], function(d) {
    tmp_reg <- adonis2(d ~ region, data = metadata_tmp, nrep = 999, by = "margin")
    tmp_sea <- adonis2(d ~ season, data = metadata_tmp, nrep = 999, by = "margin")
    tmp_inte <- adonis2(d ~ region * season, data = metadata_tmp, nrep = 999, by = "margin")
    df_reg <- data.frame(tmp_reg)
    df_sea <- data.frame(tmp_sea)
    df_reg <- df_reg[! row.names(df_reg) %in% c("Residual","Total"),]
    df_sea <- df_sea[! row.names(df_sea) %in% c("Residual","Total"),]
    df_inte <- data.frame(tmp_inte)
    out <- data.frame(rbind(df_reg, df_sea ,df_inte)) %>% round(3) %>% 
      setNames(c("df","Sum_sq","R2","F-value","p-value")) %>% 
      rownames_to_column("factor")   
  }) %>% reformat_as_df(new_var_name = "data_type")
}) %>% setNames(nms_ranks_3)

saveRDS(ls_res, file = paste0(dir_stats, "res_permanova.rds"))
ls_res <- readRDS(paste0(dir_stats, "res_permanova.rds"))


#     format the results as a table ----

df_res <- reformat_as_df(ls_res, new_var_name = "rank")
df_res$q <- factor(str_split_fixed(df_res$data_type, "_", 2)[,2])
df_res$div_facet <- factor(str_split_fixed(df_res$data_type, "_", 2)[,1])

write.csv(df_res, file = paste0(dir_stats, "table_results_permanova.csv"), row.names =FALSE)


# Pairwise comparisons for revision animal µbiome -------

ls_res <- list()

for (rk in nms_ranks_3) {
  for (idx in names(ls_diss_tmp[[rk]])) {
    d <- ls_diss_tmp[[rk]][[idx]]
    res <- pairwise.adonis2(d ~ region + season, data = metadata_tmp, nrep = 999)
    out <- lapply(res[-1], function(x) data.frame(x) %>% rownames_to_column("factor")) %>% 
      reformat_as_df(new_var_name = "comparison") %>% 
      setNames(c("Factor","DF","Sum_sq","R2","F_value","p_value", "Comparison")) %>% 
      filter(! Factor %in% c("Residual", "Total"))
    ls_res[[rk]][[idx]] <- out
  }
}

df_res_pairwise <- lapply(ls_res, function(x) { 
  x %>% reformat_as_df(new_var_name = "Index")
}) %>% reformat_as_df(new_var_name = "Rank")

df_res_out <- df_res_pairwise %>% arrange(Index, Rank) %>% 
  dplyr::select(Index, Rank, Comparison, Factor, R2, F_value, p_value, everything())


write.csv(df_res_out,
          file = paste0(dir_stats, "table_results_permanova_pairwise.csv"),
          row.names =FALSE)

out <- df_res_out %>% group_by(Comparison, Factor) %>% 
  summarize(num_signif = sum(p_value < 0.05), 
            avg_F = mean(F_value), sd_F = sd(F_value), 
            avg_p = mean(p_value), sd_p = sd(p_value),
            avg_R2 = mean(R2), sd_R2 = sd(R2)) %>% data.frame()

write.csv(out,
          file = paste0(dir_stats, "table_results_permanova_pairwise_summary.csv"),
          row.names =FALSE)


# =============== TEST DISPERSION BETWEEN GROUPS ===============================

dir_stats <- paste0(dir_out, "stats/")
dir.create(dir_stats)

# Run the tests ----

ls_res <- lapply(nms_ranks_3, function(rk) {
  lapply(ls_diss_tmp[[rk]], function(d) {
    tmp <- betadisper(d, metadata_tmp$region, type = "centroid", bias.adjust = FALSE,
                      sqrt.dist = FALSE, add = FALSE)
    return(tmp)
  }) %>% setNames(names(ls_diss_tmp$Phylum))
}) %>% setNames(nms_ranks_3)

saveRDS(ls_res, file = paste0(dir_stats, "res_permdisp.rds"))


## make the results as pretty tables

ls_res_df <- lapply(ls_res, function(X) {
  lapply(X, function(x) {
    anova(x) %>% data.frame() %>% round(3) %>% 
      setNames(c("df","Sum_sq","Mean_sq","F-value","p-value")) %>% 
      rownames_to_column("factor")    
  }) %>% reformat_as_df(new_var_name = "data_type")
}) %>% setNames(nms_ranks_3)


#     format the results as a table ----

df_res <- reformat_as_df(ls_res_df, new_var_name = "rank")
df_res$q <- factor(str_split_fixed(df_res$data_type, "_", 2)[,2])
df_res$div_facet <- factor(str_split_fixed(df_res$data_type, "_", 2)[,1])

write.csv(df_res, file = paste0(dir_stats, "table_results_permdisp.csv"), row.names =FALSE)






################################################################################
## ANALYZE FUNCTIONAL PREDICTIONS DATA =======================================


# Load the functional data ----

ls_ps_scfa <- readRDS(paste0(dir_data, "list_phyloseq_scfa.rds"))

# estimated dissimilarity
ls_diss_scfa <- readRDS(paste0(dir_data, "list_functional_dissimilarity_scfa.rds"))


dir_out <- paste0(dir_save, "functional/")
dir.create(dir_out)


#### test differences between group =======================================


dir_stats <- paste0(dir_out, "stats/")
dir.create(dir_stats)

metadata_tmp <- meta(ls_ps_scfa$all_scfa)
metadata_tmp$region <- factor(metadata_tmp$region, levels = nms_region)


#### Test all variables at once for each species  --------

nms_func_lev <- paste0("lev_", 1:3)

df_res_scfa <- lapply(names(ls_diss_scfa), function(rk) {
  
  lapply(nms_siganids[1:2], function(sp) {
    
    lapply(nms_index, function(idx) {
      
      meta_sp <- metadata_tmp %>% filter(taxonomy == sp)
      nms_sples <- meta_sp %>% dplyr::select(sample_id_fastq) %>% unlist() %>% as.vector()
      
      d <- as.matrix(ls_diss_scfa[[rk]][[idx]])
      d <- d[row.names(d) %in% nms_sples, colnames(d) %in% nms_sples]
      d <- as.dist(d)
      tmp_reg <- adonis2(d ~ region, data = meta_sp, nrep = 999, by = "margin")
      tmp_sea <- adonis2(d ~ season, data = meta_sp, nrep = 999, by = "margin")
      tmp_inte <- adonis2(d ~ region : season, data = meta_sp, nrep = 999, by = "margin")
      df_reg <- data.frame(tmp_reg)
      df_sea <- data.frame(tmp_sea)
      df_reg <- df_reg[! row.names(df_reg) %in% c("Residual","Total"),]
      df_sea <- df_sea[! row.names(df_sea) %in% c("Residual","Total"),]
      df_inte <- data.frame(tmp_inte)
      out <- data.frame(rbind(df_reg, df_sea ,df_inte)) %>% round(3) %>% 
        setNames(c("df","Sum_sq","R2","F-value","p-value")) %>% 
        rownames_to_column("factor") %>% 
        filter(! factor %in% c("Residual", "Total"))
      return(out)
    }) %>% setNames(nms_index) %>% reformat_as_df(new_var_name = "index")
  }) %>% setNames(nms_siganids[1:2]) %>% reformat_as_df(new_var_name = "species")
}) %>% setNames(names(ls_diss_scfa)) %>% reformat_as_df(new_var_name = "C_number")

write.csv(df_res_scfa, 
          file = paste0(dir_stats, "table_results_permanova_per_species.csv"), 
          row.names =FALSE)


## summarize the results ----


df_res_scfa <- read.csv(file = paste0(dir_stats, "table_results_permanova_per_species.csv"))

out <- df_res_scfa %>% 
  group_by(index, species, factor) %>% 
  summarize(avg_F = mean(`F.value`), sd_F = sd(`F.value`), 
            avg_p = mean(`p.value`), sd_p = sd(`p.value`),
            avg_R2 = mean(`R2`), sd_R2 = sd(`R2`), 
            num_signif = sum(`p.value` < 0.05)) %>% 
  arrange(species) %>% data.frame()

write.csv(out, file = paste0(dir_stats, "table_results_permanova_summary_per_species.csv"), 
          row.names =FALSE)


# Pairwise comparisons for revision animal µbiome -------

ls_res <- list()

for (rk in names(ls_diss_scfa)) {
  for (sp in nms_siganids[1:2]) {
    for (idx in nms_index[2]) {
      meta_sp <- metadata_tmp %>% filter(taxonomy == sp)
      nms_sples <- meta_sp %>% dplyr::select(sample_id_fastq) %>% unlist() %>% as.vector()
      
      d <- as.matrix(ls_diss_scfa[[rk]][[idx]])
      d <- d[row.names(d) %in% nms_sples, colnames(d) %in% nms_sples]
      d <- as.dist(d)
      res <- pairwise.adonis2(d ~ region + season, data = meta_sp, nrep = 999)[-1]
      out <- lapply(res, function(x) data.frame(x) %>% rownames_to_column("factor")) %>% 
        reformat_as_df(new_var_name = "comparison") %>% 
        setNames(c("Factor","DF","Sum_sq","R2","F_value","p_value", "Comparison")) %>% 
        filter(! Factor %in% c("Residual", "Total"))
      ls_res[[rk]][[sp]][[idx]] <- out
    }
  }
}

df_res_pairwise <- lapply(ls_res, function(X) { 
  lapply(X, function(x) {
    x %>% reformat_as_df(new_var_name = "Index")
  }) %>% reformat_as_df(new_var_name = "Species")
}) %>% reformat_as_df(new_var_name = "C_number")

df_res_scfa_out <- df_res_pairwise %>% arrange(Index, Species, C_number) %>% 
  dplyr::select(Index, Species, C_number, Comparison, Factor, R2, F_value, p_value, everything())


write.csv(df_res_scfa_out,
          file = paste0(dir_stats, "table_results_permanova_pairwise.csv"),
          row.names =FALSE)

out <- df_res_scfa_out %>% group_by(Species, Comparison, Factor) %>% 
  summarize(num_signif = sum(p_value < 0.05), 
            avg_F = mean(F_value), sd_F = sd(F_value), 
            avg_p = mean(p_value), sd_p = sd(p_value),
            avg_R2 = mean(R2), sd_R2 = sd(R2)) %>% 
  arrange(Species) %>% data.frame()

write.csv(out,
          file = paste0(dir_stats, "table_results_permanova_pairwise_summary.csv"),
          row.names =FALSE)


# Compare the species in each region separately ----

ls_res <- lapply(names(ls_diss_scfa), function(rk) {
  lapply(nms_region[-2], function(reg) {
    lapply(ls_diss_scfa[[rk]], function(d) {
      meta_reg <- metadata_tmp %>% filter(region == reg)
      nms_sples <- meta_reg %>% dplyr::select(sample_id_fastq) %>% unlist() %>% as.vector()
      
      d <- as.matrix(d)
      d <- d[row.names(d) %in% nms_sples, colnames(d) %in% nms_sples]
      d <- as.dist(d)
      tmp <- adonis2(d ~ taxonomy, data = meta_reg, nrep = 999, by = "margin")
      out <- data.frame(tmp) %>% round(3) %>% 
        setNames(c("df","Sum_sq","R2","F-value","p-value")) %>% 
        rownames_to_column("factor")    
      return(out)
    }) %>% reformat_as_df(new_var_name = "data_type")
  }) %>% setNames(nms_region[-2])
}) %>% setNames(names(ls_diss_scfa))

saveRDS(ls_res, file = paste0(dir_stats, "res_permanova_per_region.rds"))
ls_res <- readRDS(paste0(dir_stats, "res_permanova_per_region.rds"))


#     format the results as a table ----

df_res <- lapply(ls_res, function(X) {
  X %>% reformat_as_df(new_var_name = "region")
}) %>% reformat_as_df(new_var_name = "rank") %>% 
  filter(!factor %in% c("Residual","Total")) 
df_res$q <- factor(str_split_fixed(df_res$data_type, "_", 2)[,2])
df_res$div_facet <- factor(str_split_fixed(df_res$data_type, "_", 2)[,1])

write.csv(df_res, file = paste0(dir_stats, "table_results_permanova_per_region.csv"), row.names =FALSE)



# =============== TEST DISPERSION BETWEEN GROUPS ===============================

# Run the tests ----

nms_func_lev <- names(ls_diss_scfa)

ls_res <- lapply(nms_func_lev, function(rk) {
  lapply(nms_siganids[1:2], function(sp) {
    lapply(ls_diss_scfa[[rk]], function(d) {      
      meta_sp <- metadata_tmp %>% filter(taxonomy == sp)
      nms_sples <- meta_sp %>% dplyr::select(sample_id_fastq) %>% unlist() %>% as.vector()
      d <- as.matrix(d)
      d <- d[row.names(d) %in% nms_sples, colnames(d) %in% nms_sples]
      d <- as.dist(d)
      tmp <- betadisper(d, meta_sp$region, type = "centroid", bias.adjust = FALSE,
                        sqrt.dist = FALSE, add = FALSE)
      return(tmp)
    }) %>% setNames(nms_index)
  }) %>% setNames(nms_siganids[1:2])
}) %>% setNames(nms_func_lev)

saveRDS(ls_res, file = paste0(dir_stats, "res_permdisp_scfa.rds"))


## make the results as pretty tables

df_res <- lapply(ls_res, function(X) {
  lapply(X, function(x) {
    lapply(x, function(xx) {
    anova(xx) %>% data.frame() %>% round(3) %>% 
      setNames(c("df","Sum_sq","Mean_sq","F-value","p-value")) %>% 
      rownames_to_column("factor")    
    })%>% reformat_as_df(new_var_name = "data_type")
  }) %>% reformat_as_df(new_var_name = "species")
})  %>% reformat_as_df(new_var_name = "rank") %>% 
  filter(factor != "Residuals")


#     format the results as a table ----

df_res$q <- factor(str_split_fixed(df_res$data_type, "_", 2)[,2])
df_res$div_facet <- factor(str_split_fixed(df_res$data_type, "_", 2)[,1])

write.csv(df_res, file = paste0(dir_stats, "table_results_permdisp.csv"), row.names =FALSE)


#============================= PCOA ORDINATION =================================

dir_pcoa <- paste0(dir_out, "pcoa/")
dir.create(dir_pcoa)

# Make the PCOA ----


# for all the fishes simultaneously

ls_pcoa_scfa <- lapply(ls_diss_scfa, function(rk) { 
  lapply(rk, function(d) pcoa(D = as.dist(d)) )
})

saveRDS(ls_pcoa_scfa, file = paste0(dir_pcoa, "ls_pcoa_scfa.rds"))


# for each species separately 

nms_func_lev <- paste0("lev_", 1:3)

nms_scfa_C <- names(ls_diss_scfa)

ls_pcoa_scfa_per_species <- lapply(nms_siganids[1:2], function(sp) {
  
  lapply(nms_scfa_C, function(rk) {
    lapply(nms_index, function(idx) {
    
    meta_sp <- metadata_tmp %>% filter(taxonomy == sp)
    nms_sples <- meta_sp %>% dplyr::select(sample_id_fastq) %>% unlist() %>% as.vector()
    
    d <- as.matrix(ls_diss_scfa[[rk]][[idx]])
    d <- d[row.names(d) %in% nms_sples, colnames(d) %in% nms_sples]
    pcoa(D = as.dist(d))
    }) %>% setNames(nms_index) 
  }) %>% setNames(nms_scfa_C)
}) %>% setNames(nms_siganids[1:2])

saveRDS(ls_pcoa_scfa_per_species, file = paste0(dir_pcoa, "ls_pcoa_scfa_per_species.rds"))





# ============ TEST DIFFERENCES BETWEEN SITES WITHIN REGIONS ===================
# for the revision in animal µbiome

nms_idx <- c("beta_taxo", "beta_phylo")
nms_q <- c("q0", "q1")
nms_seasons <- c("Spring","Autumn")

metadata_tmp <- meta(ls_ps_core$Siganidae$ASV)
metadata_tmp$region <- factor(metadata_tmp$region, levels = nms_region)

sp <- "Siganus_rivulatus"  

df_res_rivu <-   lapply(nms_region, function(reg) {
  lapply(nms_seasons, function(season) {
    
    meta_sp <- metadata_tmp %>% filter(taxonomy == sp & region == reg & season == season)
    nms_sples <- meta_sp %>% dplyr::select(sample_id_fastq) %>% unlist() %>% as.vector()
    
    lapply(nms_ranks_3, function(rk) {
      lapply(nms_idx, function(idx) {
        lapply(nms_q, function(q) {
    
          if (nrow(meta_sp) != 0) {
            d <- as.matrix(ls_diss_fish[[sp]][[rk]][[idx]][[q]])
            d <- d[row.names(d) %in% nms_sples, colnames(d) %in% nms_sples]
            d <- as.dist(d)
            tmp <- adonis2(d ~ zone, data = meta_sp, nrep = 999, by = "margin")
            out <- data.frame(tmp) %>% round(3) %>% 
              setNames(c("df","Sum_sq","R2","F_value","p_value")) %>% 
              rownames_to_column("factor") %>% 
              filter(! factor %in% c("Residual", "Total"))
            return(out)
          }
      }) %>% setNames(nms_q) %>% reformat_as_df(new_var_name = "q_value")
    }) %>% setNames(nms_idx) %>% reformat_as_df(new_var_name = "Index")
  }) %>% setNames(nms_ranks_3) %>% reformat_as_df(new_var_name = "Rank")
  }) %>% setNames(nms_seasons) %>% reformat_as_df(new_var_name = "Season")
}) %>% setNames(nms_region) %>% reformat_as_df(new_var_name = "Region")

sp <- "Siganus_luridus"  

df_res_luri <-   lapply(nms_region[-2], function(reg) {
  
  lapply(nms_seasons, function(season) {
    
    meta_sp <- metadata_tmp %>% filter(taxonomy == sp & region == reg & season == season)
    nms_sples <- meta_sp %>% dplyr::select(sample_id_fastq) %>% unlist() %>% as.vector()
  
  lapply(nms_ranks_3, function(rk) {
    lapply(nms_idx, function(idx) {
      lapply(nms_q, function(q) {
        
        if (nrow(meta_sp) != 0) {
          d <- as.matrix(ls_diss_fish[[sp]][[rk]][[idx]][[q]])
          d <- d[row.names(d) %in% nms_sples, colnames(d) %in% nms_sples]
          d <- as.dist(d)
          tmp <- adonis2(d ~ zone, data = meta_sp, nrep = 999, by = "margin")
          out <- data.frame(tmp) %>% round(3) %>% 
            setNames(c("df","Sum_sq","R2","F_value","p_value")) %>% 
            rownames_to_column("factor") %>% 
            filter(! factor %in% c("Residual", "Total"))
          return(out)
        }
      }) %>% setNames(nms_q) %>% reformat_as_df(new_var_name = "q_value")
    }) %>% setNames(nms_idx) %>% reformat_as_df(new_var_name = "Index")
  }) %>% setNames(nms_ranks_3) %>% reformat_as_df(new_var_name = "Rank")
  }) %>% setNames(nms_seasons) %>% reformat_as_df(new_var_name = "Season")
}) %>% setNames(nms_region[-2]) %>% reformat_as_df(new_var_name = "Region")

df_res_rivu$Species <- rep("Siganus_rivulatus", nrow(df_res_rivu))
df_res_luri$Species <- rep("Siganus_luridus", nrow(df_res_luri))


df_res <- rbind(df_res_rivu, df_res_luri) %>% data.frame()
out <- df_res %>% group_by(Species, Region, Season) %>% 
  summarize(avg_p = mean(p_value),
            avg_f = mean(F_value),
            mean_r = mean(R2),
            num_sign = sum(p_value < 0.05))


write.csv(df_res, file = paste0(dir_save, "table_test_of_site_differences.csv"), row.names = FALSE)
write.csv(out, file = paste0(dir_save, "table_test_of_site_differences_summary.csv"), row.names = FALSE)


