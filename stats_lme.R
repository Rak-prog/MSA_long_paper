library(car)
library(effectsize)
library(emmeans)
library(MuMIn)
library(ggplot2)
library(lme4) 
library(lmerTest)   # p-values + df methods :contentReference[oaicite:2]{index=2}
library(performance)

doCT  <- FALSE
doVOL <- TRUE
OUTLIER_REMOVAL <- TRUE
ATLAS <- "Schaefer" # "DKT" or "Schaefer"
NPARC <- 400
NNET <- 17
method_fdr = "BH"

fmt_p <- function(p, digits = 5, min_p = 0.001) {
  # returns character vector like "0.048", "0.320", "<0.001"
  out <- ifelse(is.na(p), NA_character_,
                ifelse(p < min_p, paste0("<", format(min_p, nsmall = digits, scientific = FALSE)),
                       formatC(round(p, digits), format = "f", digits = digits)))
  out
}

if (doCT && !doVOL) {
  
  if (ATLAS == "DKT") {
    areas <- c(
      "superiorfrontal","caudalmiddlefrontal","rostralmiddlefrontal","parstriangularis","parsopercularis",
      "parsorbitalis","lateralorbitofrontal","medialorbitofrontal","precentral","paracentral",
      "superiorparietal","inferiorparietal","supramarginal","postcentral","precuneus",
      "superiortemporal","middletemporal","inferiortemporal","fusiform","transversetemporal",
      "parahippocampal","entorhinal","lateraloccipital","lingual","cuneus","pericalcarine",
      "rostralanteriorcingulate","caudalanteriorcingulate","posteriorcingulate","isthmuscingulate","insula"
    )
    p_bonf = 0.05/length(areas)
    
    if (OUTLIER_REMOVAL) {
      df <- read.csv( "/home/riccardo/codici_progetti/MSA_paper/MSA_long_paper/matrices/CT_MSAlong_data_clean_outlier_removed.csv", stringsAsFactors = FALSE)
      savename <- "/home/riccardo/codici_progetti/MSA_paper/MSA_long_paper/matrices/fsquare_MSAlong_CT_outlier_removed.csv"
    }
    else if (!OUTLIER_REMOVAL) {
      df <- read.csv( "/home/riccardo/codici_progetti/MSA_paper/MSA_long_paper/matrices/CT_MSAlong.csv", stringsAsFactors = FALSE)
      savename <- "/home/riccardo/codici_progetti/MSA_paper/MSA_long_paper/matrices/fsquare_MSAlong_CT.csv"
    }
    
  } else if (ATLAS == "Schaefer") {
    
    areas <- c(
      "VisCent","VisPeri","SomMotA","SomMotB","DorsAttnA","DorsAttnB","SalVentAttnA","SalVentAttnB",
      "LimbicA","LimbicB","ContA","ContB","ContC","DefaultA","DefaultB","DefaultC","TempPar"
    )
    p_bonf = 0.05/length(areas)
    
    if (OUTLIER_REMOVAL) {
      df <- read.csv(paste0("/home/riccardo/codici_progetti/MSA_paper/MSA_long_paper/matrices/CT_MSAlong_", ATLAS, "_Parc", NPARC, "_NET", NNET, "_data_clean_outlier_removed.csv"),
                     stringsAsFactors = FALSE)
      savename <- paste0("/home/riccardo/codici_progetti/MSA_paper/MSA_long_paper/matrices/fsquare_MSAlong_CT_",ATLAS, "_Parc", NPARC, "_NET", NNET, "_outlier_removed.csv")
    }
    else if (!OUTLIER_REMOVAL) {
      df <- read.csv(paste0("/home/riccardo/codici_progetti/MSA_paper/MSA_long_paper/matrices/CT_MSAlong_", ATLAS, "_Parc", NPARC, "_NET", NNET, ".csv"),
                     stringsAsFactors = FALSE)
      savename <- paste0("/home/riccardo/codici_progetti/MSA_paper/MSA_long_paper/matrices/fsquare_MSAlong_CT_",ATLAS, "_Parc", NPARC, "_NET", NNET, ".csv")
    }
  } else {
    stop("ATLAS must be 'DKT' or 'SCHAEFER' when doCT=TRUE.")
  }
  
} else if (!doCT && doVOL) {
  areas <- c("Putamen","CerebellumWhiteMatter","CerebellumCortex","Pons")
  p_bonf = 0.05/length(areas)
  res_pheno <- data.frame(
    ROI = areas,
    Estimate = NA_real_,
    CI95_Lower = NA_real_,
    CI95_Upper = NA_real_,
    tStat = NA_real_,
    pValue = NA_real_,
    PctDiff_vs_Ref = NA_real_,
    f2_Phenotype = NA_real_,
    stringsAsFactors = FALSE,
    Nobs = NA_integer_,
    Nsubj = NA_integer_
  )
  
  if (OUTLIER_REMOVAL) { 
    df <- read.csv( "/home/riccardo/codici_progetti/MSA_paper/MSA_long_paper/matrices/VOL_MSAlong_data_clean_outlier_removed.csv", stringsAsFactors = FALSE)
    savename <- "/home/riccardo/codici_progetti/MSA_paper/MSA_long_paper/matrices/fsquare_MSAlong_VOL_outlier_removed.csv"
    savename_phenotype <- "/home/riccardo/codici_progetti/MSA_paper/MSA_long_paper/matrices/fsquare_MSAlong_VOL_Phenotype_outlier_removed.csv"
  }
  else if (!OUTLIER_REMOVAL) {
    df <- read.csv( "/home/riccardo/codici_progetti/MSA_paper/MSA_long_paper/matrices/VOL_MSAlong.csv", stringsAsFactors = FALSE)
    savename <- "/home/riccardo/codici_progetti/MSA_paper/MSA_long_paper/matrices/fsquare_MSAlong_VOL.csv"
    savename_phenotype <- "/home/riccardo/codici_progetti/MSA_paper/MSA_long_paper/matrices/fsquare_MSAlong_VOL_Phenotype.csv"
  }  
  
} else {
  
  stop("Choose either doCT=TRUE, doVOL=FALSE OR doCT=FALSE, doVOL=TRUE")
  
}

# make sure factors are factors (as you already do)
df$Subj_ID <- as.factor(df$Subj_ID)
df$Time <- as.factor(df$Time)
df$MCI_s__1_no_0_ <- as.factor(df$MCI_s__1_no_0_)
df$SottotipoMotorio <- as.factor(df$SottotipoMotorio)
df$Sesso <- as.factor(df$Sesso)
options(contrasts = c("contr.treatment", "contr.poly"))  # to encode the difference between FUP and BL as in matlab (this is important for the parameter estimates) - using BL as reference

nROI <- length(areas)
term_regex <- "^Time"  # will match TimeT1

res <- data.frame(
  ROI = areas,
  Term = NA_character_,
  Estimate = NA_real_,
  SE = NA_real_,
  tStat = NA_real_,
  DF = NA_real_,
  pValue = NA_real_,
  CI95_Lower = NA_real_,
  CI95_Upper = NA_real_,
  Mean_T0 = NA_real_,
  Mean_T1 = NA_real_,
  Delta_T1_T0 = NA_real_,
  PctChange_vs_T0 = NA_real_,
  R2m_full = NA_real_,
  R2m_reduced = NA_real_,
  f2_Time = NA_real_,
  stringsAsFactors = FALSE
)

res$Nobs <- NA_integer_ 
res$Nsubj <- NA_integer_

for (roi in seq_along(areas)) {
  y <- areas[roi]
  df_roi <- df
  bad_subj <- unique(df_roi$Subj_ID[is.na(df_roi[[y]])])
  
  # if any, set BOTH timepoints to NA for that ROI
  if (length(bad_subj) > 0) {
    df_roi[df_roi$Subj_ID %in% bad_subj, y] <- NA
  }
  
  df_fit <- df_roi[complete.cases(df_roi[, c(y, "Time","MCI_s__1_no_0_","SottotipoMotorio","Et_","Sesso","Subj_ID")]), ] # this stay the same when remove_outliers is FALSE
  n_obs  <- nrow(df_fit)
  n_subj <- length(unique(df_fit$Subj_ID))
  
  f <- as.formula(
    paste0(y, " ~ Time + MCI_s__1_no_0_ + SottotipoMotorio + Et_ + Sesso + (1|Subj_ID)")
  )
  
  lme_fit <- lmer(f, data = df_fit, REML = FALSE)
  
  cs <- summary(lme_fit)$coefficients
  rn <- rownames(cs)
  
  time_rows <- grep(term_regex, rn, value = TRUE)
  if (length(time_rows) == 0) {
    warning("No Time term found for ROI: ", y)
    next
  }
  term <- time_rows[1]  # for 2 levels, this should be the single TimeT1 term
  
  est <- cs[term, "Estimate"]      # THIS is (T1 - T0)
  se  <- cs[term, "Std. Error"]
  tval <- cs[term, "t value"]
  
  df_matlab <- nobs(lme_fit) - length(fixef(lme_fit))
  p_matlab_like <- 2 * pt(abs(tval), df = df_matlab, lower.tail = FALSE)
  
  tcrit <- qt(0.975, df = df_matlab)
  ci_lo <- est - tcrit * se
  ci_hi <- est + tcrit * se
  
  # model-based marginal means at T0 and T1 (averaged over covariates)
  emt <- emmeans(lme_fit, ~ Time)
  em_df <- as.data.frame(emt)[, c("Time", "emmean")]
  
  lv <- levels(df_fit$Time)
  mean_t0 <- em_df$emmean[em_df$Time == lv[1]] # "T0", BL
  mean_t1 <- em_df$emmean[em_df$Time == lv[2]] # "T1", FUP
  
  delta <- est  # since est is T1 - T0 under treatment coding
  se_delta <- se 
  
  pct_change <- 100 * (delta / mean_t0)
  
  # f2 for Time
  lme_reduced_fit <- update(lme_fit, . ~ . - Time)
  r2_full <- r.squaredGLMM(lme_fit)[1, "R2m"]
  r2_red  <- r.squaredGLMM(lme_reduced_fit)[1, "R2m"]
  f2_time <- (r2_full - r2_red) / (1 - r2_full)
  
  res[roi, c("ROI","Term","Estimate","SE","tStat","DF","pValue","CI95_Lower","CI95_Upper",
             "Mean_T0","Mean_T1","Delta_T1_T0","PctChange_vs_T0","R2m_full","R2m_reduced","f2_Time")] <-
    list(y, term, est, se, tval, df_matlab, p_matlab_like, ci_lo, ci_hi,
         mean_t0, mean_t1, delta, pct_change, r2_full, r2_red, f2_time)
  
  res[roi, "Nobs"]  <- n_obs
  res[roi, "Nsubj"] <- n_subj
  
  # create the table for the phenotype effect as well when doVOL is True
  if (doVOL) {
    
    # ---- identify phenotype coefficient ----
    cs <- summary(lme_fit)$coefficients
    rn <- rownames(cs)
    
    pheno_rows <- grep("^SottotipoMotorio", rn, value = TRUE)
    if (length(pheno_rows) == 0) {
      warning("No phenotype term found for ROI: ", y)
    } else {
      
      term_ph <- pheno_rows[1]   # for 2 levels, single coefficient
      
      est_ph <- cs[term_ph, "Estimate"]
      se_ph  <- cs[term_ph, "Std. Error"]
      t_ph   <- cs[term_ph, "t value"]
      
      df_matlab <- nobs(lme_fit) - length(fixef(lme_fit))
      p_ph <- 2 * pt(abs(t_ph), df = df_matlab, lower.tail = FALSE)
      
      tcrit <- qt(0.975, df = df_matlab)
      ci_lo_ph <- est_ph - tcrit * se_ph
      ci_hi_ph <- est_ph + tcrit * se_ph
      
      # ---- model-adjusted marginal means by phenotype ----
      emp <- emmeans(lme_fit, ~ SottotipoMotorio)
      emp_df <- as.data.frame(emp)[, c("SottotipoMotorio", "emmean")]
      
      ref_lvl <- levels(df_fit$SottotipoMotorio)[1]
      mean_ref <- emp_df$emmean[emp_df$SottotipoMotorio == ref_lvl]
      
      pct_diff <- 100 * (est_ph / mean_ref)
      
      # ---- f2 for phenotype ----
      lme_red_ph <- update(lme_fit, . ~ . - SottotipoMotorio)
      r2_red_ph <- r.squaredGLMM(lme_red_ph)[1, "R2m"]
      f2_ph <- (r2_full - r2_red_ph) / (1 - r2_full)
      
      
      res_pheno[roi, c("ROI","Estimate","CI95_Lower","CI95_Upper","tStat","pValue","PctDiff_vs_Ref","f2_Phenotype")] <-
        list(y, est_ph, ci_lo_ph, ci_hi_ph, t_ph, p_ph, pct_diff, f2_ph)
      
      res_pheno[roi, "Nobs"]  <- n_obs
      res_pheno[roi, "Nsubj"] <- n_subj
    }
  }
}

res_final <- res[, c("ROI","Nsubj","Estimate","CI95_Lower","CI95_Upper","tStat","pValue","f2_Time","PctChange_vs_T0")]
res_final$P_txt <- fmt_p(res_final$pValue, digits = 3, min_p = 0.001)
# adding flags for reliability
res_final$Sig_FDR  <- p.adjust(res_final$pValue, method = method_fdr) < 0.05
res_final$Sig_Bonf <- res_final$pValue < p_bonf

res_final$Estimate         <- round(res_final$Estimate, 3)
res_final$pValue         <- round(res_final$pValue, 5)
res_final$CI95_Lower       <- round(res_final$CI95_Lower, 3)
res_final$CI95_Upper       <- round(res_final$CI95_Upper, 3)
res_final$tStat            <- round(res_final$tStat, 2)
res_final$PctChange_vs_T0  <- round(res_final$PctChange_vs_T0, 2)
res_final$f2_Time          <- round(res_final$f2_Time, 3)

print(res_final)
res_final_sorted <- res_final[order(res_final$PctChange_vs_T0), ]
print(res_final_sorted)
write.csv(res_final, file = savename, row.names = FALSE)
savename_sorted <- sub("\\.csv$", "_sorted.csv", savename)
write.csv(res_final_sorted, file = savename_sorted, row.names = FALSE)

print(savename)
if (doVOL) {
  
  # round numeric columns
  num_cols <- sapply(res_pheno, is.numeric)
  res_pheno[num_cols] <- round(res_pheno[num_cols], 2)
  
  # formatted p-values (reuse your helper)
  res_pheno$P_txt <- fmt_p(res_pheno$pValue, digits = 3, min_p = 0.001)
  
  # significance flags
  res_pheno$Sig_FDR  <- p.adjust(res_pheno$pValue, method = method_fdr) < 0.05
  res_pheno$Sig_Bonf <- res_pheno$pValue < p_bonf
  res_pheno$pValue <- res_pheno$P_txt
  res_pheno$P_txt
  res_pheno[ , c("P_txt")]
  write.csv(res_pheno, savename_phenotype, row.names = FALSE)
  print(res_pheno)
  print(ref_lvl)
}

 # sort the final tables according to % change 


 