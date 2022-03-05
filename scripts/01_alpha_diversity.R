################################################################################
#  _______  _____  _____ ___ ____  _   _ __  __ _____ ____  
# | ____\ \/ / _ \|  ___|_ _/ ___|| | | |  \/  | ____|  _ \ 
# |  _|  \  / | | | |_   | |\___ \| |_| | |\/| |  _| | | | |
# | |___ /  \ |_| |  _|  | | ___) |  _  | |  | | |___| |_| |
# |_____/_/\_\___/|_|   |___|____/|_| |_|_|  |_|_____|____/ 
#
################################################################################
# R script to analyse alpha diversity
# arthur.escalas@gmail.com
################################################################################


############################### LOAD DATA ######################################

dir_save <- dir_alpha_div

# Load the phyloseq objects with diversity estimates ----

# the full rarefied dataset to compare sample types
ps_full <- readRDS(paste0(dir_data, "phyloseq_object_cleaned_with_tree_rarefied_with_alpha_div.rds"))

# The fish core dataset to compare fishes
ls_ps_core <- readRDS(paste0(dir_data, "list_phyloseq_object_fish_core_with_alpha_diversity.rds"))



############### TEST DIFFERENCES BETWEEN ECOSYSTEM COMPARTMENTS ################


# ============== Compute summary statistics for alpha diversity ================

tab <- meta(ps_full)

#     for each sample type globally ----

res_summary_alpha_spletype <- lapply(index_alpha[1:2], function(idx) {
  tmp <- tab %>% mutate(datavec = tab[, idx])
  res <- tmp %>% 
    dplyr::group_by(sample_type) %>% 
    dplyr::summarize(avg = mean(datavec),
                     sd  = sd(datavec),
                     min = min(datavec),
                     max = max(datavec)) %>% 
    mutate(index = rep(idx, length(avg)))
return(res)
}) %>% bind_rows()

write.csv(res_summary_alpha_spletype, 
          paste0(dir_save, "res_table_summary_alphadiv_per_habitat.csv"))
  
  
#     for each sample_type and region ----

res_summary_alpha_spletype_region <- lapply(index_alpha[1:2], function(idx) {
  tmp <- tab %>% mutate(datavec = tab[, idx])
  res <- tmp %>% 
    dplyr::group_by(sample_type, region) %>% 
    dplyr::summarize(avg = mean(datavec),
                     sd  = sd(datavec),
                     min = min(datavec),
                     max = max(datavec)) %>% 
    mutate(index = rep(idx, length(avg)))
  return(res)
}) %>% bind_rows()

write.csv(res_summary_alpha_spletype_region, 
          paste0(dir_save, "res_table_summary_alphadiv_per_habitat_&_region.csv"))

res_summary_alpha_spletype_regionseason <- lapply(index_alpha[1:2], function(idx) {
  tmp <- tab %>% mutate(datavec = tab[, idx])
  res <- tmp %>% 
    dplyr::group_by(sample_type, region, season) %>% 
    dplyr::summarize(avg = mean(datavec),
                     sd  = sd(datavec),
                     min = min(datavec),
                     max = max(datavec)) %>% 
    mutate(index = rep(idx, length(avg)))
  return(res)
}) %>% bind_rows()

write.csv(res_summary_alpha_spletype_region, 
          paste0(dir_save, "res_table_summary_alphadiv_per_habitat_&_region_&_season.csv"))




# ================= Plot alpha diversity in each compartment ===================

# objects for plotting ----

nm_cols <- c(Water = "deepskyblue2", Sediment = "darkgoldenrod4", 
             Turf = "gold2", Algae = "chartreuse4", Seagrass = "limegreen",
             Fish = "firebrick3")
nm_cols <- c(Water = "#1874CD", Sediment = "#8B6508", 
             Turf = "#EEC900", Algae = "#458B00", Seagrass = "#32CD32",
             Fish = "#EE7600")
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



# make the plot ----

png(paste0(dir_save, "figure_alpha_diversity.png"),
    height = 25, width = 20, unit = "cm", res = 400)

tab <- meta(ps_full)
tab$region <- factor(tab$region, levels = nms_region)
tab$season <- factor(tab$season, levels = c("Spring", "Autumn"))
tab$region_season <- factor(tab$region_season, levels = nms_reg_seas)
tab$sample_type <- factor(tab$sample_type, levels = c("Water","Sediment",
                                                      "Turf","Algae","Seagrass",
                                                      "Fish"))

layout(matrix(c(1,1,2,3,4,5,6,7,8,9), 5,2, TRUE))
par(mar = c(3,3,0,2), oma = c(3,3,2,1), las = 1, 
    mgp = c(2.5,0.5,0), xpd = NA, font.lab = 2, tcl = -0.3, yaxs = "i")

idx <- "taxo_q1"
ylim = c(0, 800)
nm_div_idx <- "Taxonomic entropy (q = 1)"

# sample type ----

fact <- tab$sample_type

boxplot(tab[,idx] ~ tab$sample_type, col = paste0(nm_cols, 70), 
        main = "", names = rep("", length(levels(fact))),
        outline = FALSE, ylim = ylim, xlab = "", ylab = nm_div_idx)
sapply(1:length(levels(fact)), function(l) {
  f <- levels(fact)[l]
  mask <- fact == f
  mask[is.na(mask)] <- FALSE
  xval <- rep(l, sum(mask)) + sample(seq(-0.2,0.2, 0.001), sum(mask))
  points(tab[mask, idx] ~ xval, pch = 21, cex = 0.6,
         col = paste0(nm_cols[f]), bg = paste0(nm_cols[f]))  
})
text(1:length(levels(fact)), rep(min(ylim) - 0.1 * diff(ylim), 3),
     labels = levels(fact), font = 2, cex = 1.2)
mtext("Sample type", 3, font = 2)

# separately for each compartment ----

par(mar = c(2,3,1,2))
for (i in levels(tab$sample_type)[c(1:4)]) {
  
  tmp <- tab %>% filter(sample_type == i)
  fact <- tmp$region_season
  
  boxplot(tmp[,idx] ~ tmp$region_season, col = paste0(cols_regions, 70), 
          border = cols_border, names = rep("", length(levels(fact))),
          main = "", outline = FALSE, ylim = ylim, xlab = "", ylab = nm_div_idx)
  sapply(1:length(levels(fact)), function(l) {
    f <- levels(fact)[l]
    mask <- fact == f
    mask[is.na(mask)] <- FALSE
    xval <- rep(l, sum(mask)) + sample(seq(-0.2,0.2, 0.001), sum(mask))
    points(tmp[mask, idx] ~ xval, pch = 21, cex = 0.6,
           col = paste0(cols_border[f]), bg = paste0(cols_border[f]))  
  })
  mtext(i, 3, font = 2)
}


# separately for each fish species ----

tab <- meta(ls_ps_core$Siganidae)
tab$region <- factor(tab$region, levels = nms_region)
tab$season <- factor(tab$season, levels = c("Spring", "Autumn"))
nms_reg_seas <- c("North_Red_Sea_Spring", "North_Red_Sea_Autumn", 
                  "Levantine_Sea_Spring", "Levantine_Sea_Autumn", 
                  "Northern_Crete_Spring", "Northern_Crete_Autumn")
tab$region_season <- factor(tab$region_season, levels = nms_reg_seas)

nm_div_idx <- "Taxonomic entropy\nof the core microbiome (q = 1)"
ylim = c(0, 300)

for (i in nms_siganids[1:2]) {
  tmp <- tab %>% filter(taxonomy == i)
  fact <- tmp$region_season
  
  boxplot(tmp[,idx] ~ tmp$region_season, col = paste0(cols_regions, 70), 
          border = cols_border, names = rep("", length(levels(fact))),
          main = "", outline = FALSE, ylim = ylim, xlab = "", ylab = nm_div_idx)
  sapply(1:length(levels(fact)), function(l) {
    f <- levels(fact)[l]
    mask <- fact == f
    mask[is.na(mask)] <- FALSE
    xval <- rep(l, sum(mask)) + sample(seq(-0.2,0.2, 0.001), sum(mask))
    points(tmp[mask, idx] ~ xval, pch = 21, cex = 0.6,
           col = paste0(cols_border[f]), bg = paste0(cols_border[f]))  
  })
  mtext(gsub("_", " ", i), 3, font = 4)
}

idx <- "phylo_q1"
nm_div_idx <- "Phylogenetic entropy\nof the core microbiome (q = 1)"

ylim = c(0, 15)
  
for (i in nms_siganids[1:2]) {
  tmp <- tab %>% filter(taxonomy == i)
  fact <- tmp$region_season
  
  boxplot(tmp[,idx] ~ tmp$region_season, col = paste0(cols_regions, 70), 
          border = cols_border, names = rep("", length(levels(fact))),
          main = "", outline = FALSE, ylim = ylim, xlab = "", ylab = nm_div_idx)
  sapply(1:length(levels(fact)), function(l) {
    f <- levels(fact)[l]
    mask <- fact == f
    mask[is.na(mask)] <- FALSE
    xval <- rep(l, sum(mask)) + sample(seq(-0.2,0.2, 0.001), sum(mask))
    points(tmp[mask, idx] ~ xval, pch = 21, cex = 0.6,
           col = paste0(cols_border[f]), bg = paste0(cols_border[f]))  
  })
  mtext(gsub("_", " ", i), 3, font = 4)
  
  text(seq(1.5, by = 2, length.out = 3), rep(min(ylim) - 0.25 * diff(ylim), 3),
       labels = gsub("_", " ", nms_region), font = 2, 
       col = col_reg, cex = 1.1)
  text(1:length(levels(fact)), rep(min(ylim) - 0.1 * diff(ylim), 3),
       labels = c("Spring"," Autumn"), font = 1)
}

dev.off()

# ----



# ================ Test differences between compartments =======================
# PB: as we have different number of samples per groups so we bootstrap by
# randomly taking a similar number of samples per group and repeating the analysis
# the average statistic and p-values will be reported


#     Test the differences in alpha diversity between COMPARTMENTS ------

meta(ps_full) %>% group_by(sample_type) %>% tally()

res_kw_compartment <- lapply(index_alpha[1:2], function(q) {
  get_bootstraped_KW_or_W_test(meta(ps_full), "sample_type", 100, q)  
})

res_dunn_compartment <- lapply(index_alpha[1:2], function(q) {
  get_bootstraped_Dunn_test(meta(ps_full), "sample_type", 100, q)  
})


#   Test the differences in alpha dversity between REGIONS for each compartment ----

meta(ps_full) %>% group_by(sample_type, region) %>% tally()

ls_data <- split(meta(ps_full), meta(ps_full)$sample_type)
ls_data <- ls_data[c("Algae", "Fish", "Sediment", "Turf", "Water")]

# global effect
res_kw_reg_per_comp <- lapply(ls_data, function(X) {
  lapply(index_alpha[1:2], function(q) {
    get_bootstraped_KW_or_W_test(X, "region", 100, q)  
  })
})

# pairwise test
# we don't test for the habitat for which there is only 2 regions
mask <- names(ls_data)[-c(1,3)]

res_dunn_reg_per_comp <- lapply(ls_data[mask], function(X) {
  lapply(index_alpha[1:2], function(q) {
    get_bootstraped_Dunn_test(X, "region", 100, q)  
  })
})


#     Test the differences in alpha dversity between SEASONS for each compartment ----

meta(ps_full) %>% group_by(sample_type, season) %>% tally()

# global effect
res_kw_seas_per_comp <- lapply(ls_data, function(X) {
  lapply(index_alpha[1:2], function(q) {
    get_bootstraped_KW_or_W_test(X, "season", 100, q)  
  })
})


#     Test the differences in alpha diversity between REGION_SEASONS for each compartment ----

meta(ps_full) %>% group_by(sample_type, region, season) %>% tally() %>% data.frame()

# global effect
res_kw_regseas_per_comp <- lapply(ls_data[-c(1,3)], function(X) {
  lapply(index_alpha[1:2], function(q) {
    get_bootstraped_KW_or_W_test(X, "region_season", 100, q)  
  })
})

# pairwise test
# we don't test for the habitat for which there is only 2 regions

res_dunn_regseas_per_comp <- lapply(ls_data[-c(1,3)], function(X) {
  lapply(index_alpha[1:2], function(q) {
    get_bootstraped_Dunn_test(X, "region_season", 100, q)  
  })
})


res_aov_regseas_per_comp <- lapply(ls_data[-1], function(X) {
  lapply(index_alpha[1:2], function(q) {
    get_bootstraped_2way_ANOVA(X, "region", "season", 100, q)  
  })
})



# ----------------------- SUMMARIZE THE RESULTS --------------------------------

dir_stats <- paste0(dir_save, "stats/")
dir.create(dir_stats)

# between habitats
tab_kw <- reformat_as_df(res_kw_compartment, new_var_name = "div_index") %>% 
  group_by(div_index) %>% 
  dplyr::summarize(avg_stat = mean(as.numeric(test_stat)),
                   sd_stat  = sd(test_stat),
                   avg_pval = mean(as.numeric(p_value)),
                   sd_pval  = sd(p_value),
                   num_sign = sum(as.numeric(p_value) < 0.05))

write.csv(tab_kw, paste0(dir_stats, "table_kw_sample_type.csv"))

tab_dunn <- lapply(res_dunn_compartment, function(X) {
  X %>% group_by(comparisons) %>% 
  dplyr::summarize(avg_z = mean(as.numeric(Z)),
                   sd_z  = sd(Z),
                   avg_pval = mean(as.numeric(P)),
                   sd_pval  = sd(P),
                   num_sign = sum(as.numeric(P) < 0.05)) %>% 
    arrange(dplyr::desc(num_sign), avg_pval)
}) %>% reformat_as_df(new_var_name = "div_index") %>% data.frame()

write.csv(tab_dunn, paste0(dir_stats, "table_dunn_sample_type.csv"))

# seasons effect

tab_kw <- lapply(res_kw_seas_per_comp, function(X) {
  reformat_as_df(X, new_var_name = "div_index") %>%
    group_by(div_index) %>% 
    dplyr::summarize(avg_stat = mean(as.numeric(test_stat)),
                     sd_stat  = sd(test_stat),
                     avg_pval = mean(as.numeric(p_value)),
                     sd_pval  = sd(p_value),
                     num_sign = sum(as.numeric(p_value) < 0.05)) %>% 
    arrange(dplyr::desc(num_sign), avg_pval)
}) %>% reformat_as_df(new_var_name = "habitat") %>% 
  arrange(div_index, habitat) %>% 
  dplyr::select(div_index, habitat, num_sign, everything()) %>% data.frame()

write.csv(tab_kw, paste0(dir_stats, "table_kw_sample_type_season.csv"))


# region effect


tab_kw <- lapply(res_kw_reg_per_comp, function(X) {
  reformat_as_df(X, new_var_name = "div_index") %>% 
    group_by(div_index) %>% 
    dplyr::summarize(avg_stat = mean(as.numeric(test_stat)),
                     sd_stat  = sd(test_stat),
                     avg_pval = mean(as.numeric(p_value)),
                     sd_pval  = sd(p_value),
                     num_sign = sum(as.numeric(p_value) < 0.05))  
}) %>% reformat_as_df(new_var_name = "habitat") %>% 
  arrange(div_index, habitat, dplyr::desc(num_sign)) %>% 
  dplyr::select(div_index, habitat, num_sign, everything()) %>% data.frame()

write.csv(tab_kw, paste0(dir_stats, "table_kw_sample_type_region.csv"))

tab_dunn <- lapply(res_dunn_reg_per_comp, function(X) {
  reformat_as_df(X, new_var_name = "div_index") %>% 
    group_by(div_index, comparisons) %>% 
    dplyr::summarize(avg_stat = mean(as.numeric(Z)),
                     sd_stat  = sd(Z),
                     avg_pval = mean(as.numeric(P)),
                     sd_pval  = sd(P),
                     num_sign = sum(as.numeric(P) < 0.05))  %>% ungroup()
}) %>% reformat_as_df(new_var_name = "habitat") %>% 
  arrange(div_index, habitat, dplyr::desc(num_sign)) %>% 
  dplyr::select(div_index, habitat, num_sign, everything()) %>% data.frame()

write.csv(tab_dunn, paste0(dir_stats, "table_dunn_sample_type_region.csv"))


# region and season effect

tab_kw <- lapply(res_kw_regseas_per_comp, function(X) {
  reformat_as_df(X, new_var_name = "div_index") %>% 
  group_by(div_index) %>% 
    dplyr::summarize(avg_stat = mean(as.numeric(test_stat)),
                     sd_stat  = sd(test_stat),
                     avg_pval = mean(as.numeric(p_value)),
                     sd_pval  = sd(p_value),
                     num_sign = sum(as.numeric(p_value) < 0.05))  
}) %>% reformat_as_df(new_var_name = "habitat") %>% 
  arrange(div_index, habitat, dplyr::desc(num_sign)) %>% 
  dplyr::select(div_index, habitat, num_sign, everything()) %>% data.frame()

write.csv(tab_kw, paste0(dir_stats, "table_kw_sample_type_regionXseason.csv"))


tab_dunn <- lapply(res_dunn_regseas_per_comp, function(X) {
  reformat_as_df(X, new_var_name = "div_index") %>% 
    group_by(div_index, comparisons) %>% 
    dplyr::summarize(avg_stat = mean(as.numeric(Z)),
                     sd_stat  = sd(Z),
                     avg_pval = mean(as.numeric(P)),
                     sd_pval  = sd(P),
                     num_sign = sum(as.numeric(P) < 0.05))  %>% ungroup()
}) %>% reformat_as_df(new_var_name = "habitat") %>% 
  arrange(div_index, habitat) %>% 
  dplyr::select(div_index, habitat, num_sign, everything()) %>% data.frame()

write.csv(tab_dunn, paste0(dir_stats, "table_dunn_sample_type_regionXseason.csv"))


tab_aov_region_x_season <- lapply(res_aov_regseas_per_comp, function(X) {
  lapply(X, function(x) {
   x %>% group_by(factor) %>% 
    dplyr::summarize(avg_f = mean(as.numeric(f_value)),
                     sd_f  = sd(f_value),
                     avg_pval = round(mean(as.numeric(p_value)),4),
                     sd_pval  = sd(p_value),
                     avg_r2 = mean(as.numeric(r_squared)),
                     sd_r2  = sd(r_squared),
                     num_sign = sum(as.numeric(p_value) < 0.05))  
  }) %>% reformat_as_df(new_var_name = "div_index")
}) %>% reformat_as_df(new_var_name = "habitat") %>% 
  mutate(factor = factor(factor, levels = c("region","season","region:season"))) %>% 
  arrange(div_index, habitat, factor) %>% 
  dplyr::select(div_index, habitat, num_sign, everything()) %>% data.frame()

write.csv(tab_aov_region_x_season, paste0(dir_stats, "table_anova_sample_type_regionXseason.csv"))





##################### TEST DIFFERENCES BETWEEN FISHES ##########################


# ============== Compute summary statistics for alpha diversity ================

tab <- do.call(rbind, lapply(nms_siganids, function(X) meta(ls_ps_core[[X]])))

#     for each species ----

res_summary_alpha_fish_region <- lapply(index_alpha, function(idx) {
  tmp <- tab %>% mutate(datavec = tab[, idx])
  res <- tmp %>% 
    dplyr::group_by(sample_type_3, region) %>% 
    dplyr::summarize(avg = mean(datavec),
                     sd  = sd(datavec),
                     min = min(datavec),
                     max = max(datavec)) %>% 
    mutate(index = rep(idx, length(avg)))
  return(res)
}) %>% bind_rows() %>% data.frame() %>% 
  mutate(region = factor(region, levels = c("North_Red_Sea", "Levantine_Sea",
                                            "Northern_Crete"))) %>% 
  arrange(sample_type_3, index, region)

write.csv(res_summary_alpha_fish_region, 
          paste0(dir_save, "res_table_summary_alphadiv_fish_per_region.csv"))


#     for each sample_type and region ----

res_summary_alpha_fish_region_season <- lapply(index_alpha, function(idx) {
  tmp <- tab %>% mutate(datavec = tab[, idx])
  res <- tmp %>% 
    dplyr::group_by(sample_type_3, region, season) %>% 
    dplyr::summarize(avg = mean(datavec),
                     sd  = sd(datavec),
                     min = min(datavec),
                     max = max(datavec)) %>% 
    mutate(index = rep(idx, length(avg)))
  return(res)
}) %>% bind_rows() %>% data.frame() %>% 
  mutate(region = factor(region, levels = c("North_Red_Sea", "Levantine_Sea",
                                            "Northern_Crete"))) %>% 
  arrange(sample_type_3, index, region)

write.csv(res_summary_alpha_fish_region_season, 
          paste0(dir_save, "res_table_summary_alphadiv_fish_per_region_season.csv.csv"))



# ============== Test differences between sites within regions ==================

ls_data <- lapply(ls_ps_core, microbiome::meta)


# global effect
res_kw_site_effect <- lapply(ls_data[-c(3)], function(X) {
  lapply(index_alpha, function(q) {
    lapply(split(X, X$region), function(x) {
      get_bootstraped_KW_or_W_test(x, "zone", 100, q)  
    })
  })
})


tab_kw <- lapply(res_kw_site_effect, function(X) {
  lapply(X, function(x) {
    reformat_as_df(x, new_var_name = "Region") %>% 
      group_by(Region) %>% 
      dplyr::summarize(avg_stat = mean(as.numeric(test_stat)),
                       sd_stat  = sd(test_stat),
                       avg_pval = mean(as.numeric(p_value)),
                       sd_pval  = sd(p_value),
                       num_sign = sum(as.numeric(p_value) < 0.05))  
    }) %>% reformat_as_df(new_var_name = "div_index") 
  }) %>% reformat_as_df(new_var_name = "Species")  %>% 
  arrange(Species, Region, div_index, dplyr::desc(num_sign)) %>% 
  dplyr::select(Species, Region, div_index, num_sign, everything()) %>% data.frame()

write.csv(tab_kw, paste0(dir_stats, "table_kw_wilcox_site_effect.csv"))



# ============== Test differences between regions and seasons ==================

ls_data <- lapply(ls_ps_core, microbiome::meta)

res_aov_regseas_per_comp <- lapply(ls_data, function(X) {
  lapply(index_alpha, function(q) {
    get_bootstraped_2way_ANOVA(X, "region", "season", 100, q)  
  })
})

tab_aov_region_x_season <- lapply(res_aov_regseas_per_comp, function(X) {
  lapply(X, function(x) {
    x %>% group_by(factor) %>% 
      dplyr::summarize(avg_f = mean(as.numeric(f_value)),
                       sd_f  = sd(f_value),
                       avg_pval = round(mean(as.numeric(p_value)),4),
                       sd_pval  = sd(p_value),
                       avg_r2 = mean(as.numeric(r_squared)),
                       sd_r2  = sd(r_squared),
                       num_sign = sum(as.numeric(p_value) < 0.05))  
  }) %>% reformat_as_df(new_var_name = "div_index")
}) %>% reformat_as_df(new_var_name = "species") %>% 
  mutate(factor = factor(factor, levels = c("region","season","region:season"))) %>% 
  arrange(species, div_index, factor) %>% 
  dplyr::select(div_index, species, num_sign, everything()) %>% data.frame()

write.csv(tab_aov_region_x_season, paste0(dir_stats, "table_anova_fish_regionXseason.csv"))



