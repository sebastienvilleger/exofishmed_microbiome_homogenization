################################################################################
#  _______  _____  _____ ___ ____  _   _ __  __ _____ ____  
# | ____\ \/ / _ \|  ___|_ _/ ___|| | | |  \/  | ____|  _ \ 
# |  _|  \  / | | | |_   | |\___ \| |_| | |\/| |  _| | | | |
# | |___ /  \ |_| |  _|  | | ___) |  _  | |  | | |___| |_| |
# |_____/_/\_\___/|_|   |___|____/|_| |_|_|  |_|_____|____/ 
#
################################################################################
#
# Functions to manipulate and analyze data
#
# Arthur Escalas January 2019
# arthur.escalas@gmail.com
################################################################################


#################################################### CCA ANALYSIS

# Schlaeppi et al., PNAS, 2013

# variability_table generates a table with the relative variability explained 
# by each transformation given a CCA object

variability_table <- function(cca){
  
  chi <- c(cca$tot.chi, cca$CCA$tot.chi, cca$CA$tot.chi)
  variability_table <- cbind(chi, chi / chi[1])
  colnames(variability_table) <- c("inertia", "proportion")
  rownames(variability_table) <- c("total", "constrained", "unconstrained") 
  return(variability_table)
  
}
# cap_var_props returns the variance explained by the CPCoA transformation 
# for a given CCA object

cap_var_props <- function(cca){
  
  eig_tot <- sum(cca$CCA$eig)
  var_propdf <- cca$CCA$eig / eig_tot
  return(round(var_propdf,3))
}

# cca_ci returns the confidence interval for the variance explained by a 
# CPCoA transformation for a given CCA object using a ANOVA-based permutation 
# test with 5000 permutations by default

cca_ci <- function(cca, permutations = 5000) {
  
  var_tbl <- variability_table(cca)
  p <- permutest(cca, permutations=permutations)
  ci <- quantile(p$F.perm, 
                 c(.05, .95)) * p$chi[1] / var_tbl["total", "inertia"]
  return(ci)
  
}



################################################  BOOTSTRAPPED STATISTICAL TESTS


# functions to make analyses on bootrstraped data ===============================

get_bootstraped_KW_or_W_test <- function(dat, fact, n_perm = 100, y_var = "hn_q0") {
  
  min_num <- min(table(dat[, fact]))
  
  ls_gp <- split(dat, dat[, fact])
  
  ls_sple_combs <- lapply(1:n_perm, function(x) {
    do.call(c, lapply(ls_gp, function(xx) {
      sample(xx$sample_id_fastq, min_num)
    }))
  })
  
  df_res_test <- do.call(rbind, lapply(ls_sple_combs, function(nm_sples) {
    tab <- dat %>% filter(sample_id_fastq %in% nm_sples)
    if (length(ls_gp) > 2) {
      tmp <- kruskal.test(tab[, y_var] ~ tab[, fact])
      out <- c("KW", tmp$statistic, tmp$parameter, tmp$p.value)
    } else{
      tmp <- wilcox.test(tab[, y_var] ~ tab[, fact])
      out <- c("W", tmp$statistic, NA, tmp$p.value)
    }
    out
  })) %>% data.frame() 
  names(df_res_test) <- c("test_type", "test_stat", "df", "p_value")
  
  return(df_res_test)
}

get_bootstraped_Dunn_test <- function(dat, fact, n_perm = 100, y_var = "hn_q0") {
  
  min_num <- min(table(dat[, fact]))
  
  ls_gp <- split(dat, dat[, fact])
  
  ls_sple_combs <- lapply(1:n_perm, function(x) {
    do.call(c, lapply(ls_gp, function(xx) {
      sample(xx$sample_id_fastq, min_num)
    }))
  })
  
  df_res_test <- lapply(ls_sple_combs, function(nm_sples) {
    tab <- dat %>% filter(sample_id_fastq %in% nm_sples)
    tmp <- dunn.test::dunn.test(tab[, y_var], tab[, fact], method = "bh")
    do.call(cbind, tmp[c(5,2,3)]) %>% data.frame()
  }) %>% setNames(paste0("rdm_", 1:n_perm)) %>% 
    reformat_as_df(new_var_name = "permutation")
  
  return(df_res_test)
}

get_bootstraped_1way_ANOVA <- function(dat, fact, n_perm = 100, y_var = "hn_q0") {
  
  min_num <- min(table(dat[, fact]))
  
  ls_gp <- split(dat, dat[, fact])
  
  ls_sple_combs <- lapply(1:n_perm, function(x) {
    do.call(c, lapply(ls_gp, function(xx) {
      sample(xx$sample_id_fastq, min_num)
    }))
  })
  
  df_res_test <- lapply(ls_sple_combs, function(nm_sples) {
    tab <- dat %>% filter(sample_id_fastq %in% nm_sples)
    fit <- lm(formula(paste0(y_var, "~", fact)), data = tab)
    tmp <- summary(fit)
    tmp2 <- anova(fit)
    out <- do.call(cbind, list("1-way-ANOVA", tmp$adj.r.squared,  
                               tmp2[1, c("Df", "F value", "Pr(>F)")]))
    names(out) <- c("test_type", "r_squared", "df", "f_value", "p_value")
    out
  }) %>% setNames(paste0("rdm_", 1:n_perm)) %>% 
    reformat_as_df(new_var_name = "permutation")

  return(df_res_test)
}

get_bootstraped_2way_ANOVA <- function(dat, fact1, fact2, n_perm = 100, y_var = "hn_q0") {
  
  fact_comb <- apply(dat[, c(fact1, fact2)], 1, function(x) {
    paste(x, collapse = "_")
  })
  min_num <- min(table(fact_comb))
  
  ls_gp <- split(dat, fact_comb)
  
  ls_sple_combs <- lapply(1:n_perm, function(x) {
    do.call(c, lapply(ls_gp, function(xx) {
      sample(xx$sample_id_fastq, min_num)
    }))
  })
  
  df_res_test <- lapply(ls_sple_combs, function(nm_sples) {
    tab <- dat %>% filter(sample_id_fastq %in% nm_sples)
    fit <- lm(formula(paste0(y_var, "~", fact1, "*", fact2)), data = tab) 
    tmp <- summary(fit)
    tmp2 <- anova(fit) %>% data.frame() %>% 
      rownames_to_column(var = "factor")%>% 
      filter(factor != "Residuals")
    
    out <- cbind(rep("2-way-ANOVA", 3), rep(tmp$adj.r.squared,3), tmp2)
    names(out) <- c("test_type", "r_squared", "factor",  "df", "sum_sq", 
                    "mean_sq",  "f_value", "p_value")
    out
  }) %>% setNames(paste0("rdm_", 1:n_perm)) %>% 
    reformat_as_df(new_var_name = "permutation")
  
  return(df_res_test)
}







