################################################################################
#  _______  _____  _____ ___ ____  _   _ __  __ _____ ____  
# | ____\ \/ / _ \|  ___|_ _/ ___|| | | |  \/  | ____|  _ \ 
# |  _|  \  / | | | |_   | |\___ \| |_| | |\/| |  _| | | | |
# | |___ /  \ |_| |  _|  | | ___) |  _  | |  | | |___| |_| |
# |_____/_/\_\___/|_|   |___|____/|_| |_|_|  |_|_____|____/ 
#
################################################################################
#
# R script to anlyze dissimilarity between fishes based on
# - taxonomy/phylogeny of the microbiome
# - functions
# arthur.escalas@gmail.com
################################################################################


#     Create a file for today' analyses ----

dir_save <- dir_dissimilarity


########################### LOAD AND PREPARE DATA ################################

# =========================  microbiome data ===================================

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


# reformat the list of dissimilarity for easier analysis ----

ls_diss_tmp <- lapply(ls_diss_fish$Siganidae, function(rk) {
  out <- list()
  out$taxo_q0 <- rk$beta_taxo$q0
  out$taxo_q1 <- rk$beta_taxo$q1
  out$phylo_q0 <- rk$beta_phylo$q0
  out$phylo_q1 <- rk$beta_phylo$q1
  out
})

# µbiome metadata
metadata_µ <- meta(ls_ps_core$Siganidae$Phylum)

df_alldiss <- readRDS(paste(dir_data, "df_all_dissimilarities.rds"))
ls_df <- readRDS(paste(dir_data, "list_dissimilarity_df.rds"))
ls_diss <- readRDS(paste(dir_data, "list_dissimilarity.rds"))

nms_ranks_3 <- c("Phylum", "Family", "ASV")




######################## ANALYSIS OF DISSIMILARITY #############################


# ========== create the labels of the different combinations of factors ========

# dataframe with samples combinations ----

df <- ls_diss_tmp$Phylum$taxo_q0 %>% 
  as.dist() %>% dendextend::dist_long() %>% 
  dplyr::select(- distance) %>% 
  setNames(c("sample_a","sample_b"))

# add metadata to this dataframe ----

df$spl_comb <- sapply(1:nrow(df), function(x) { 
  tmp <- sort(c(as.character(df$sample_a)[x], 
                as.character(df$sample_b)[x]))
  paste0(tmp[1], "__", tmp[2])
})
df$sp_comb <- sapply(df$spl_comb, function(x) {
  sapply(str_split_fixed(x, "__", 2), function(xx) {
    metadata_µ[metadata_µ$sample_id_fastq == xx, "taxonomy"]
  }) %>% sort() %>% paste0(collapse = "__")
})
df$reg_comb <- sapply(df$spl_comb, function(x) {
  sapply(str_split_fixed(x, "__", 2), function(xx) {
    metadata_µ[metadata_µ$sample_id_fastq == xx, "region"]
  }) %>% sort() %>% paste0(collapse = "__")
})

# reformat and clean the df ----

df <- df %>% 
  separate(sp_comb, c("species_a", "species_b"), sep = "__", remove = FALSE) %>%
  separate(reg_comb, c("region_a", "region_b"), sep = "__", remove = FALSE) %>%
  dplyr::select(sample_a, sample_b, species_a, species_b, region_a, region_b, 
                spl_comb, sp_comb, reg_comb,
                everything())

df$sp_lev <- rep("intra_species", nrow(df))
df$sp_lev[df$species_a != df$species_b] <- "inter_species"
df$reg_lev <- rep("intra_region", nrow(df))
df$reg_lev[df$region_a != df$region_b] <- "inter_region"

# change factor levels
for (i in names(df)[names(df) != "distance"]) {
  df[, i] <- factor(df[,i], levels = sort(unique(df[,i])))
}

# clean the keys used to split data -----

key <- group_by(df, sp_lev, reg_lev, sp_comb, reg_comb) %>% 
  group_keys() %>% data.frame()
key <- apply(key, 2, as.character)

for (i in 1:nrow(key)) {  
  x <- key[i,]
  if (length(grep("intra", x["sp_lev"])) == 1) {
    x["sp_comb"] <- strsplit(x[ "sp_comb"], "__")[[1]][1]
  }
  if (length(grep("inter", x["sp_lev"])) == 1) {
    x["sp_comb"] <- NA
  }
  if (length(grep("intra", x["reg_lev"])) == 1) {
    x["reg_comb"] <- strsplit(x[ "reg_comb"], "__")[[1]][1]
  }
  key[i,] <- x
}

# create data names
nms <- apply(key, 1, function(x) {
  x <- paste0(x[! is.na(x)], collapse = "_|_")
  gsub("__", "-", x)
})


# =======================  Summarize the trends ================================


# average  diss for each combination of factors ----

df_summ_alldiss <- df_alldiss %>% 
  group_by(sp_lev, reg_lev, sp_comb, reg_comb) %>% 
  summarize_at(names(df_alldiss)[12:32], 
               list(avg = mean, sd = sd), na.rm = TRUE) %>% data.frame()


# make the dataframe tidy ----------------

df_alldiss <- df_alldiss %>% 
  pivot_longer(cols = where(is.numeric), names_to = c("index", "rank"),
               names_sep = "__") %>% data.frame() %>% 
  separate("index", into = c("facet","order"), sep = "_", remove = FALSE) %>% 
  rename(value = "dissimilarity")


# summarize dissimilarity by combinations of factors ----

df_summ <- df_alldiss %>% 
  group_by(sp_lev, reg_lev, sp_comb, reg_comb, index, facet, order, rank) %>% 
  summarize(avg_diss = mean(dissimilarity, na.rm = TRUE),
            sd_diss = sd(dissimilarity, na.rm = TRUE)) %>% data.frame()



# estimate the drop in dissimilarity in the non native range -----

# intra species intra regions ----

df_summ %>% filter(reg_lev == "intra_region" & sp_lev == "intra_species") %>% 
  dplyr::select(-sd_diss) %>% 
  pivot_wider(names_from = reg_comb, values_from = avg_diss) %>% 
  mutate(red_levant = (Levantine_Sea__Levantine_Sea - North_Red_Sea__North_Red_Sea) / North_Red_Sea__North_Red_Sea * 100,
         red_crete = (Northern_Crete__Northern_Crete - North_Red_Sea__North_Red_Sea) / North_Red_Sea__North_Red_Sea * 100) %>% 
  data.frame() %>% group_by(sp_comb, index) %>% 
  summarize(avg_red_levant = mean(red_levant, na.rm = TRUE), 
            sd_red_levant = sd(red_levant, na.rm = TRUE), 
            avg_red_crete = mean(red_crete, na.rm = TRUE),
            sd_red_crete = sd(red_crete, na.rm = TRUE))
out <- df_summ %>% filter(reg_lev == "intra_region" & sp_lev == "intra_species") %>% 
  dplyr::select(-sd_diss) %>% 
  pivot_wider(names_from = reg_comb, values_from = avg_diss) %>% 
  mutate(red_levant = (Levantine_Sea__Levantine_Sea - North_Red_Sea__North_Red_Sea) / North_Red_Sea__North_Red_Sea * 100,
         red_crete = (Northern_Crete__Northern_Crete - North_Red_Sea__North_Red_Sea) / North_Red_Sea__North_Red_Sea * 100) %>% 
  data.frame() %>%  dplyr::select(facet, rank, order, sp_comb,
                                  North_Red_Sea__North_Red_Sea,
                                  Levantine_Sea__Levantine_Sea,
                                  Northern_Crete__Northern_Crete,
                                  red_levant, red_crete)

names(out) <- c("Facet","Rank","Order", "Species", 
                "North Red Sea", "Levantine Sea", "Northern Crete", 
                "Change in dissimilarity red_levant",
                "Change in dissimilarity red_crete")
write.csv(out, file = paste0(dir_save, "table_loss_dissimiliarity_intra_species.csv"),
          row.names = FALSE)


# inter species intra regions ----


df_summ %>% filter(reg_lev == "intra_region" & sp_lev == "inter_species") %>% 
  dplyr::select(-sd_diss) %>% 
  pivot_wider(names_from = reg_comb, values_from = avg_diss) %>% 
  mutate(red_crete = (Northern_Crete__Northern_Crete - North_Red_Sea__North_Red_Sea) / North_Red_Sea__North_Red_Sea * 100) %>% 
  data.frame() %>% 
  summarize(mean(red_crete),sd(red_crete))

out <- df_summ %>% filter(reg_lev == "intra_region" & sp_lev == "inter_species") %>% 
  dplyr::select(-sd_diss) %>% 
  pivot_wider(names_from = reg_comb, values_from = avg_diss) %>% 
  mutate(red_crete = (Northern_Crete__Northern_Crete - North_Red_Sea__North_Red_Sea) / North_Red_Sea__North_Red_Sea * 100) %>% 
  data.frame() %>% dplyr::select(rank, facet, order, North_Red_Sea__North_Red_Sea,
                                 Northern_Crete__Northern_Crete,
                                 red_crete)

names(out) <- c("Rank","Facet","Order","North Red Sea","Northern Crete", "Change in dissimilarity red_crete")

write.csv(out, file = paste0(dir_save, "table_loss_dissimiliarity_inter_species.csv"),
          row.names = FALSE)





################################################################################
# --------------------- Figure Distance to centroid ----------------------

dir_plot_centroid <- paste0(dir_save, "plots_centroid/")
dir.create(dir_plot_centroid)

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

# plot the distances to centroids ------
cols_reg <- col_reg
nms_index <- c("q0", "q1")

for (idx in nms_index) {
  
  for (rk in nms_ranks_3) {
    
    jpeg(paste0(dir_plot_centroid, "boxplot_distance_to_centroids_environment_", idx,"_per_", rk, ".jpeg"),
         height = 15, width = 20, unit = "cm", res = 600)
    
    par(mfrow = c(2,2), mar = c(2,2,1,1), las = 1, mgp = c(2.5,0.5,0),
        oma = c(2,3,1,1))
    for (x in names(ls_permdisp)[1:4]) {
      tmp <- ls_permdisp[[x]][[rk]][[idx]]
      ylim <- setrge(tmp$distances, 5)
      fact <- tmp$group
      boxplot(tmp$distances ~ fact, outline = FALSE, main = firstup(x), ylim = ylim,
              names = rep("", length(levels(fact))), 
              border = paste0(cols_reg[levels(fact)],"80"), ylab = "",
              col = "white", lwd = 2)
      sapply(1:length(levels(fact)), function(l) {
        f <- levels(fact)[l]
        mask <- fact == f
        mask[is.na(mask)] <- FALSE
        xval <- rep(l, sum(mask)) + sample(seq(-0.25,0.25, 0.001), sum(mask))
        points(tmp$distances[mask] ~ xval, pch = 21, cex = 0.6,
               col = cols_reg[f],
               bg = cols_reg[f])
      })
    }
    mtext(side = 2, text = "Distance to centroid", font = 2,
          outer = TRUE, las = 0, line = 1)    
    legend(-1.5,min(ylim - 0.1*diff(ylim)), legend = gsub("_", " ", nms_region), pch = 22, ncol = 3,
           bty = "n", col = col_reg, pt.bg = col_reg, pt.cex = 2, cex = 1.2,
           text.col = col_reg, xpd = NA)
    dev.off()
    
  }
}


for (rk in nms_ranks_3) {
  
  jpeg(paste0(dir_plot_centroid, "boxplot_distance_to_centroids_fishes_all_idx_per_", rk, ".jpeg"),
       height = 25, width = 20, unit = "cm", res = 400)
  
  par(mfrow = c(4,2), mar = c(2,2,1,1), las = 1, mgp = c(2.5,0.5,0),
      oma = c(2,3,1,1))
  for (idx in names(ls_permdisp$s_rivulatus$Phylum)) {
    
    for (x in names(ls_permdisp)[5:6]) {
      tmp <- ls_permdisp[[x]][[rk]][[idx]]
      ylim <- setrge(tmp$distances, 5)
      fact <- tmp$group
      boxplot(tmp$distances ~ fact, outline = FALSE, main = paste0(firstup(x)," (",idx,")"), 
              ylim = ylim,
              names = rep("", length(levels(fact))), 
              border = paste0(cols_reg[levels(fact)],"80"), ylab = "",
              col = "white", lwd = 2)
      sapply(1:length(levels(fact)), function(l) {
        f <- levels(fact)[l]
        mask <- fact == f
        mask[is.na(mask)] <- FALSE
        xval <- rep(l, sum(mask)) + sample(seq(-0.25,0.25, 0.001), sum(mask))
        points(tmp$distances[mask] ~ xval, pch = 21, cex = 0.6,
               col = cols_reg[f],
               bg = cols_reg[f])
      })
    }
  }
  mtext(side = 2, text = "Distance to centroid", font = 2,
        outer = TRUE, las = 0, line = 1)    
  legend(-.5,min(ylim - 0.1*diff(ylim)), legend = gsub("_", " ", nms_region), pch = 22, ncol = 3,
         bty = "n", col = col_reg, pt.bg = col_reg, pt.cex = 2, cex = 1.2,
         text.col = col_reg, xpd = NA)
  
  dev.off()
  
}




################################################################################
# --------------------- Figure for the paper FONCTIONNEL ----------------------

nms_index <- c("q0", "q1")

dir_plot_functio <- paste0(dir_save, "plots_functio/")
dir.create(dir_plot_functio)

for (idx in nms_index) {

tmp <- df_alldiss %>% filter(facet == "func" & order == idx)
ls_df <- split(tmp, tmp$rank)

nms <- apply(key, 1, function(x) {
  x <- paste0(x[! is.na(x)], collapse = "_|_")
  gsub("__", "-", x)
})

ls_diss <- lapply(ls_df, function(x) {
    x <- x %>% 
      group_split(sp_lev, reg_lev, sp_comb, reg_comb) %>% setNames(nms)
    ls <- lapply(x, function(xx) xx$dissimilarity)
    ls[rev(c(13,10,4,12,14,11,5,7,8,9,6))]
})

saveRDS(ls_diss, file = paste(dir_save, "list_dissimilarity_fonctio_", idx, ".rds"))


# Make dissimilarity plots for each index and rank ----

ls_restest <- list()

  for (rk in names(ls_df)) {
    
    # define axes labels ----  
    ylab <- "Functional dissimilarity" 
    
    nms <- gsub("_", " ", nms_region)
    nms_pairs <- c(paste0(nms[1], "\nvs.\n", nms[2]),
                   paste0(nms[1], "\nvs.\n", nms[3]), 
                   paste0(nms[2], "\nvs.\n", nms[3]))
    
    # make the plot ----
    jpeg(paste0(dir_plot_functio, "figure_dissimilarity_functional_",idx, "_", rk, ".jpeg"),
         width = 20, height = 22, unit = "cm", res = 600)
    
    # layout(matrix(c(1,3,2,4,5,6), 3, 2, TRUE))
    layout(matrix(c(1,1,3,3,2,2,4,4,6,5,5,6), 3, 4, TRUE))
    par(mar = c(2,2,3,1), las = 1, oma = c(3,4,2,1), xpd = NA, 
        bty = "n", tcl = -0.3, font.lab = 2, cex.lab = 1.2,mgp = c(2.5,0.5,0))
    
    dat <- ls_diss[[rk]]
    
    # SIGANUS RIVULATUS ----
    
    # intra region ----
    at <- c(1,2,3)
    sel <- c(11,8,7)
    
    # prepare data for test
    tmp <- dat[sel]
    tab <- do.call(c, tmp)
    g   <- do.call(c, lapply(1:length(tmp), function(x) rep(x, length(tmp[[x]]))))
    if (length(sel[!is.na(sel)]) == 2) {
      restest <- wilcox.test(tmp[[1]], tmp[[2]])
    } else {
      restest <- kruskal.test(tab, g)
    }
    ls_restest$intra_species$Siganus_rivulatus[[rk]] <- restest
    
    col_bord <- c("#FF6347","#63B8FF","#CDCD00")
    col_fill <- c("#FFFFFF","#FFFFFF","#FFFFFF")
    boxplot(dat[sel], at = at, xlim = c(0.5,3.5),
            ylab = paste0("Within regions\n", ylab), ylim = c(0,1), xaxt = "n",
            col = col_fill,
            border = col_bord, lwd = 2,
            main = "", add = FALSE, outline = FALSE,
            names = rep("", length(sel)))
    text(at, rep(-0.15,3), labels = c("North\nRed Sea", "Levantine\nSea", "Northern\nCrete"), 
         font = 2, cex = 1)
    text(0.5, 1.05, labels = "A", font = 1, cex = 2)
    text(4, 1.2, labels = "Intra-species dissimilarity", font = 2, cex = 1.5)
    text(2,1.05, "Siganus rivulatus", font = 3, cex = 1.5)
    lapply(1:length(sel), function(l) {
      if (! is.na(sel[l])) {
        xval <- rep(l, length(dat[[sel[l]]])) + sample(seq(-0.2,0.2, 0.0001), length(dat[[sel[l]]]))
        points(dat[[sel[l]]] ~ xval, pch = 21, cex = 0.3,
               col = paste0(col_reg[l], 20),
               bg = paste0(col_reg[l], 20))
      }
    })
    
    # inter region ----
    at <- c(1,2,3)
    sel <- c(4,2,3)
    
    # prepare data for test
    tmp <- dat[sel]
    tab <- do.call(c, tmp)
    g   <- do.call(c, lapply(1:length(tmp), function(x) rep(x, length(tmp[[x]]))))
    
    col_bord <- c("#FF6347","#63B8FF","#FF6347")
    col_fill <- c("#FFFFFF","#FFFFFF","#FFFFFF")
    boxplot(dat[sel], at = at, xlim = c(0.5,3.5),
            ylab = paste0("Between regions\n", ylab), ylim = c(0,1), xaxt = "n",
            col = col_fill,
            border = col_bord, lwd = 2,
            main = "", add = FALSE, outline = FALSE,
            names = rep("", length(sel)))
    col_bord <- c("#FF6347","#FF6347","#63B8FF")
    col_fill <- c("#63B8FF","#CDCD00", "#CDCD00")
    boxplot(dat[sel], at = at, xlim = c(0.5,3.5),
            ylab = "", ylim = c(0,1), xaxt = "n",
            col = col_fill,
            border = col_bord,lwd = 2,
            main = "", add = TRUE, outline = FALSE,
            names = rep("", length(sel)))
    text(0.5, 1.05, labels = "C", font = 1, cex = 2)
    text(at, rep(-0.15,3), labels = nms_pairs, font = 2, cex = 1)
    lapply(1:length(sel), function(l) {
      if (! is.na(sel[l])) {
        xval <- rep(at[l], length(dat[[sel[l]]])) + sample(seq(-0.2,0.2, 0.00001), 
                                                           length(dat[[sel[l]]]))
        points(dat[[sel[l]]] ~ xval, pch = 21, cex = 0.3,
               col = paste0(col_bord[l], 20),
               bg = paste0(col_bord[l], 20))
      }
    })
    
    
    # ----
    
    # SIGANUS LURIDUS ----
    
    # intra region ----
    at <- c(1,2,3)
    sel <- c(10, NA, 6)
    
    # prepare data for test
    tmp <- dat[sel]
    tab <- do.call(c, tmp)
    g   <- do.call(c, lapply(1:length(tmp), function(x) rep(x, length(tmp[[x]]))))
    if (length(sel[!is.na(sel)]) == 2) {
      restest <- wilcox.test(tmp[[1]], tmp[[2]])
    } else {
      restest <- kruskal.test(tab, g)
    }
    ls_restest$intra_species$Siganus_luridus[[rk]] <- restest
    
    col_bord <- c("#FF6347","","#CDCD00")
    col_fill <- c("#FFFFFF","","#FFFFFF")
    boxplot(dat[sel], at = at, xlim = c(0.5,3.5),
            ylab = "", ylim = c(0,1), xaxt = "n",
            col = col_fill, 
            border = col_bord, lwd = 2,
            main = "", add = FALSE, outline = FALSE,
            names = rep("", length(sel)))
    text(0.5, 1.05, labels = "B", font = 1, cex = 2)
    text(at, rep(-0.15,3), labels = c("North\nRed Sea", "", "Northern\nCrete"), 
         font = 2, cex = 1)    
    text(2,1.05, "Siganus luridus", font = 3, cex = 1.5)
    lapply(1:length(sel), function(l) {
      if (! is.na(sel[l])) {
        xval <- rep(l, length(dat[[sel[l]]])) + sample(seq(-0.2,0.2, 0.00001), length(dat[[sel[l]]]))
        points(dat[[sel[l]]] ~ xval, pch = 21, cex = 0.3,
               col = paste0(col_reg[l], 20),
               bg = paste0(col_reg[l], 20))
      }
    })
    
    # inter region ----
    at <- c(2)
    sel <- c(1)
    
    # prepare data for test
    tmp <- dat[sel]
    tab <- do.call(c, tmp)
    g   <- do.call(c, lapply(1:length(tmp), function(x) rep(x, length(tmp[[x]]))))
    
    col_bord <- c("#FF6347","#63B8FF","#FF6347")
    col_fill <- c("#FFFFFF","#FFFFFF","#FFFFFF")
    boxplot(dat[sel], at = at, xlim = c(0.5,3.5),
            ylab = "", ylim = c(0,1), xaxt = "n",
            col = col_fill, width = 1,
            border = col_bord, lwd = 2,
            main = "", add = FALSE, outline = FALSE,
            names = rep("", length(sel)))
    col_bord <- c("#FF6347","#FF6347","#63B8FF")
    col_fill <- c("#CDCD00","#CDCD00", "#CDCD00")
    boxplot(dat[sel], at = at, xlim = c(0.5,3.5),
            ylab = "", ylim = c(0,1), xaxt = "n",
            col = col_fill,width = 1,
            border = col_bord,lwd = 2,
            main = "", add = TRUE, outline = FALSE,
            names = rep("", length(sel)))
    text(0.5, 1.05, labels = "D", font = 1, cex = 2)
    text(at, rep(-0.15,3), labels = c("", nms_pairs[2], ""),
         font = 2, cex = 1)
    lapply(1:length(sel), function(l) {
      if (! is.na(sel[l])) {
        xval <- rep(at[l], length(dat[[sel[l]]])) + sample(seq(-0.2,0.2, 0.00001), length(dat[[sel[l]]]))
        points(dat[[sel[l]]] ~ xval, pch = 21, cex = 0.3,
               col = paste0(col_bord[l], 20),
               bg = paste0(col_bord[l], 20))
      }
    })
    
    
    # ----
    
    # inter species HOMOGENIZATION ----
    at <- c(1,2,3)
    sel <- c(9, NA, 5)
    par(mar = c(1,3,4,1))
    # prepare data for test
    tmp <- dat[sel]
    tab <- do.call(c, tmp)
    g   <- do.call(c, lapply(1:length(tmp), function(x) rep(x, length(tmp[[x]]))))
    if (length(sel[!is.na(sel)]) == 2) {
      restest <- wilcox.test(tmp[[1]], tmp[[2]])
    } else {
      restest <- kruskal.test(tab, g)
    }
    ls_restest$inter_species[[rk]] <- restest
    
    col_bord <- c("#FF6347","#63B8FF","#CDCD00")
    col_fill <- c("#FFFFFF","#FFFFFF","#FFFFFF")
    boxplot(dat[sel], at = at, xlim = c(0.5,3.5),
            ylab = paste0("Within region\n", ylab), ylim = c(0,1), xaxt = "n",
            col = col_fill,
            border = col_bord, lwd = 2,
            main = "", add = FALSE, outline = FALSE,
            names = rep("", length(sel)))
    text(at, rep(-0.15,3), labels = c("North\nRed Sea", "", "Northern\nCrete"), 
         font = 2, cex = 1)   
    text(0.5, 1.05, labels = "E", font = 1, cex = 2)
    text(2,1.05, "Inter-species dissimilarity", font = 2, cex = 1.5)
    
    lapply(1:length(sel), function(l) {
      if (! is.na(sel[l])) {
        xval <- rep(l, length(dat[[sel[l]]])) + sample(seq(-0.2,0.2, 0.00001), length(dat[[sel[l]]]))
        points(dat[[sel[l]]] ~ xval, pch = 21, cex = 0.3,
               col = paste0(col_reg[l], 20),
               bg = paste0(col_reg[l], 20))
      }
    }) 
    

    dev.off()
    
}
}

# Analyse the results of tests ----

df_restest <- lapply(ls_restest$intra_species, function(X) {
  lapply(X, function(x) {
      data.frame(stat = x$statistic, p_value = x$p.value, method = x$method)  
    }) %>% reformat_as_df(new_var_name = "rank")
}) %>% reformat_as_df(new_var_name = "species") %>% 
  dplyr::select(species, rank, everything()) %>% 
  arrange(species, rank)

df_restest$method <- gsub("Wilcoxon signed rank exact test",
                          "Wilcoxon test", df_restest$method)
df_restest$method <- gsub("Wilcoxon signed rank test with continuity correction",
                          "Wilcoxon test", df_restest$method)
df_restest$method <- gsub("Kruskal-Wallis rank sum test",
                          "Kruskal test", df_restest$method)

write.csv(df_restest, file = paste0(dir_save, "table_test_intra_sp_functio.csv"),
          row.names = FALSE)

df_restest <- lapply(ls_restest$inter_species, function(x) {
    data.frame(stat = x$statistic, p_value = x$p.value, method = x$method)  
  }) %>% reformat_as_df(new_var_name = "rank") %>% 
  dplyr::select(rank, everything()) %>% 
  arrange(rank)

df_restest$method <- gsub("Wilcoxon signed rank test with continuity correction",
                          "Wilcoxon test", df_restest$method)

write.csv(df_restest, file = paste0(dir_save, "table_test_inter_sp_functio.csv"),
          row.names = FALSE)




################################################################################
# ----------------------- Figure for the paper TAXO/PHYLO ------------------------------

dir_plot_taxo <- paste0(dir_save, "plots_taxo/")
dir.create(dir_plot_taxo)

ls_diss <- readRDS(file = paste(dir_data, "list_dissimilarity.rds"))

# Make dissimilarity plots for each index and rank ----

ls_restest <- list()

for (idx in index_alpha) {
  
  for (rk in c("Phylum","Family","ASV")) {
    
    # define axes labels ----  
    if (length(grep("phylo", idx)) == 1) { 
      a <- "Phylogenetic dissimilarity" 
    } else { 
      a <- "Taxonomic dissimilarity" 
    }
    ylab <- a
    
    nms <- gsub("_", " ", nms_region)
    nms_pairs <- c(paste0(nms[1], "\nvs.\n", nms[2]),
                   paste0(nms[1], "\nvs.\n", nms[3]), 
                   paste0(nms[2], "\nvs.\n", nms[3]))
    
    # make the plot ----
    jpeg(paste0(dir_plot_taxo, "figure_dissimilarity_", rk, "_", idx,".jpeg"),
         width = 20, height = 22, unit = "cm", res = 600)
    
    # layout(matrix(c(1,3,2,4,5,6), 3, 2, TRUE))
    layout(matrix(c(1,1,3,3,2,2,4,4,6,5,5,6), 3, 4, TRUE))
    par(mar = c(2,2,3,1), las = 1, oma = c(3,4,2,1), xpd = NA, 
        bty = "n", tcl = -0.3, font.lab = 2, cex.lab = 1.2,mgp = c(2.5,0.5,0))
    
    dat <- ls_diss[[idx]][[rk]]
    
    # SIGANUS RIVULATUS ----
    
    # intra region ----
    at <- c(1,2,3)
    sel <- c(11,8,7)
    
    # prepare data for test
    tmp <- dat[sel]
    tab <- do.call(c, tmp)
    g   <- do.call(c, lapply(1:length(tmp), function(x) rep(x, length(tmp[[x]]))))
    if (length(sel[!is.na(sel)]) == 2) {
      restest <- wilcox.test(tmp[[1]], tmp[[2]])
    } else {
      restest <- kruskal.test(tab, g)
    }
    ls_restest$intra_species$Siganus_rivulatus[[idx]][[rk]] <- restest
    
    col_bord <- c("#FF6347","#63B8FF","#CDCD00")
    col_fill <- c("#FFFFFF","#FFFFFF","#FFFFFF")
    boxplot(dat[sel], at = at, xlim = c(0.5,3.5),
            ylab = paste0("Within regions\n", ylab), ylim = c(0,1), xaxt = "n",
            col = col_fill,
            border = col_bord, lwd = 2,
            main = "", add = FALSE, outline = FALSE,
            names = rep("", length(sel)))
    text(at, rep(-0.15,3), labels = c("North\nRed Sea", "Levantine\nSea", "Northern\nCrete"), 
         font = 2, cex = 1)
    text(0.5, 1.05, labels = "A", font = 1, cex = 2)
    text(4, 1.2, labels = "Intra-species dissimilarity", font = 2, cex = 1.5)
    text(2,1.05, "Siganus rivulatus", font = 3, cex = 1.5)
    lapply(1:length(sel), function(l) {
      if (! is.na(sel[l])) {
        xval <- rep(l, length(dat[[sel[l]]])) + sample(seq(-0.2,0.2, 0.0001), length(dat[[sel[l]]]))
        points(dat[[sel[l]]] ~ xval, pch = 21, cex = 0.3,
               col = paste0(col_reg[l], 20),
               bg = paste0(col_reg[l], 20))
      }
    })
    
    # inter region ----
    at <- c(1,2,3)
    sel <- c(4,2,3)
    
    # prepare data for test
    tmp <- dat[sel]
    tab <- do.call(c, tmp)
    g   <- do.call(c, lapply(1:length(tmp), function(x) rep(x, length(tmp[[x]]))))
    
    col_bord <- c("#FF6347","#63B8FF","#FF6347")
    col_fill <- c("#FFFFFF","#FFFFFF","#FFFFFF")
    boxplot(dat[sel], at = at, xlim = c(0.5,3.5),
            ylab = paste0("Between regions\n", ylab), ylim = c(0,1), xaxt = "n",
            col = col_fill,
            border = col_bord, lwd = 2,
            main = "", add = FALSE, outline = FALSE,
            names = rep("", length(sel)))
    col_bord <- c("#FF6347","#FF6347","#63B8FF")
    col_fill <- c("#63B8FF","#CDCD00", "#CDCD00")
    boxplot(dat[sel], at = at, xlim = c(0.5,3.5),
            ylab = "", ylim = c(0,1), xaxt = "n",
            col = col_fill,
            border = col_bord,lwd = 2,
            main = "", add = TRUE, outline = FALSE,
            names = rep("", length(sel)))
    text(0.5, 1.05, labels = "C", font = 1, cex = 2)
    text(at, rep(-0.15,3), labels = nms_pairs, font = 2, cex = 1)
    lapply(1:length(sel), function(l) {
      if (! is.na(sel[l])) {
        xval <- rep(at[l], length(dat[[sel[l]]])) + sample(seq(-0.2,0.2, 0.0001), length(dat[[sel[l]]]))
        points(dat[[sel[l]]] ~ xval, pch = 21, cex = 0.3,
               col = paste0(col_bord[l], 20),
               bg = paste0(col_bord[l], 20))
      }
    })
    
    
    # ----
    
    # SIGANUS LURIDUS ----
    
    # intra region ----
    at <- c(1,2,3)
    sel <- c(10, NA, 6)
    
    # prepare data for test
    tmp <- dat[sel]
    tab <- do.call(c, tmp)
    g   <- do.call(c, lapply(1:length(tmp), function(x) rep(x, length(tmp[[x]]))))
    if (length(sel[!is.na(sel)]) == 2) {
      restest <- wilcox.test(tmp[[1]], tmp[[2]])
    } else {
      restest <- kruskal.test(tab, g)
    }
    ls_restest$intra_species$Siganus_luridus[[idx]][[rk]] <- restest
    
    col_bord <- c("#FF6347","","#CDCD00")
    col_fill <- c("#FFFFFF","","#FFFFFF")
    boxplot(dat[sel], at = at, xlim = c(0.5,3.5),
            ylab = "", ylim = c(0,1), xaxt = "n",
            col = col_fill, 
            border = col_bord, lwd = 2,
            main = "", add = FALSE, outline = FALSE,
            names = rep("", length(sel)))
    text(0.5, 1.05, labels = "B", font = 1, cex = 2)
    text(at, rep(-0.15,3), labels = c("North\nRed Sea", "", "Northern\nCrete"), 
         font = 2, cex = 1)    
    text(2,1.05, "Siganus luridus", font = 3, cex = 1.5)
    lapply(1:length(sel), function(l) {
      if (! is.na(sel[l])) {
        xval <- rep(l, length(dat[[sel[l]]])) + sample(seq(-0.2,0.2, 0.0001), length(dat[[sel[l]]]))
        points(dat[[sel[l]]] ~ xval, pch = 21, cex = 0.3,
               col = paste0(col_reg[l], 20),
               bg = paste0(col_reg[l], 20))
      }
    })
    
    # inter region ----
    at <- c(2)
    sel <- c(1)
    
    # prepare data for test
    tmp <- dat[sel]
    tab <- do.call(c, tmp)
    g   <- do.call(c, lapply(1:length(tmp), function(x) rep(x, length(tmp[[x]]))))
    
    col_bord <- c("#FF6347","#63B8FF","#FF6347")
    col_fill <- c("#FFFFFF","#FFFFFF","#FFFFFF")
    boxplot(dat[sel], at = at, xlim = c(0.5,3.5),
            ylab = "", ylim = c(0,1), xaxt = "n",
            col = col_fill, width = 1,
            border = col_bord, lwd = 2,
            main = "", add = FALSE, outline = FALSE,
            names = rep("", length(sel)))
    col_bord <- c("#FF6347","#FF6347","#63B8FF")
    col_fill <- c("#CDCD00","#CDCD00", "#CDCD00")
    boxplot(dat[sel], at = at, xlim = c(0.5,3.5),
            ylab = "", ylim = c(0,1), xaxt = "n",
            col = col_fill,width = 1,
            border = col_bord,lwd = 2,
            main = "", add = TRUE, outline = FALSE,
            names = rep("", length(sel)))
    text(0.5, 1.05, labels = "D", font = 1, cex = 2)
    text(at, rep(-0.15,3), labels = c("", nms_pairs[2], ""),
         font = 2, cex = 1)
    lapply(1:length(sel), function(l) {
      if (! is.na(sel[l])) {
        xval <- rep(at[l], length(dat[[sel[l]]])) + sample(seq(-0.2,0.2, 0.0001), length(dat[[sel[l]]]))
        points(dat[[sel[l]]] ~ xval, pch = 21, cex = 0.3,
               col = paste0(col_bord[l], 20),
               bg = paste0(col_bord[l], 20))
      }
    })
    
    
    # ----
    
    # inter species HOMOGENIZATION ----
    at <- c(1,2,3)
    sel <- c(9, NA, 5)
    par(mar = c(1,3,4,1))
    # prepare data for test
    tmp <- dat[sel]
    tab <- do.call(c, tmp)
    g   <- do.call(c, lapply(1:length(tmp), function(x) rep(x, length(tmp[[x]]))))
    if (length(sel[!is.na(sel)]) == 2) {
      restest <- wilcox.test(tmp[[1]], tmp[[2]])
    } else {
      restest <- kruskal.test(tab, g)
    }
    ls_restest$inter_species[[idx]][[rk]] <- restest
    
    col_bord <- c("#FF6347","#63B8FF","#CDCD00")
    col_fill <- c("#FFFFFF","#FFFFFF","#FFFFFF")
    boxplot(dat[sel], at = at, xlim = c(0.5,3.5),
            ylab = paste0("Within region\n", ylab), ylim = c(0,1), xaxt = "n",
            col = col_fill,
            border = col_bord, lwd = 2,
            main = "", add = FALSE, outline = FALSE,
            names = rep("", length(sel)))
    text(at, rep(-0.15,3), labels = c("North\nRed Sea", "", "Northern\nCrete"), 
         font = 2, cex = 1)   
    text(0.5, 1.05, labels = "E", font = 1, cex = 2)
    text(2,1.05, "Inter-species dissimilarity", font = 2, cex = 1.5)
    
    lapply(1:length(sel), function(l) {
      if (! is.na(sel[l])) {
        xval <- rep(l, length(dat[[sel[l]]])) + sample(seq(-0.2,0.2, 0.0001), length(dat[[sel[l]]]))
        points(dat[[sel[l]]] ~ xval, pch = 21, cex = 0.3,
               col = paste0(col_reg[l], 20),
               bg = paste0(col_reg[l], 20))
      }
    }) 
    
    dev.off()
    
  }
}


# Analyse the results of tests ----

df_restest <- lapply(ls_restest$intra_species, function(X) {
  lapply(X, function(XX) {
    lapply(XX, function(x) {
      data.frame(stat = x$statistic, p_value = x$p.value, method = x$method)  
    }) %>% reformat_as_df(new_var_name = "rank")
  }) %>% reformat_as_df(new_var_name = "index")
}) %>% reformat_as_df(new_var_name = "species") %>% 
  dplyr::select(species, rank, index, everything()) %>% 
  arrange(species, rank, index)

df_restest$method <- gsub("Wilcoxon signed rank exact test",
                          "Wilcoxon test", df_restest$method)
df_restest$method <- gsub("Wilcoxon signed rank test with continuity correction",
                          "Wilcoxon test", df_restest$method)
df_restest$method <- gsub("Kruskal-Wallis rank sum test",
                          "Kruskal test", df_restest$method)

write.csv(df_restest, file = paste0(dir_save, "table_test_intra_sp.csv"),
          row.names = FALSE)

df_restest <- lapply(ls_restest$inter_species, function(XX) {
  lapply(XX, function(x) {
    data.frame(stat = x$statistic, p_value = x$p.value, method = x$method)  
  }) %>% reformat_as_df(new_var_name = "rank")
}) %>% reformat_as_df(new_var_name = "index") %>% 
  dplyr::select(rank, index, everything()) %>% 
  arrange(rank, index)

df_restest$method <- gsub("Wilcoxon signed rank test with continuity correction",
                          "Wilcoxon test", df_restest$method)

write.csv(df_restest, file = paste0(dir_save, "table_test_inter_sp.csv"),
          row.names = FALSE)


