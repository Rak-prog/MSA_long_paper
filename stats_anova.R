library(car)
library(effectsize)
library(ggplot2)

# the partial eta of this script matches the MatLab results (HC script plues do_between group anova) both p values and partial value correspond 

doCT  <- TRUE
doVOL <- FALSE 
ATLAS <- "Schaefer" #"DKT" or Schaefer
NPARC <- 400
NNET <- 17

fmt_p <- function(p, digits = 3, min_p = 0.001) {
  ifelse(is.na(p), NA_character_,
         ifelse(p < min_p,
                paste0("<", format(min_p, nsmall = digits, scientific = FALSE)),
                formatC(round(p, digits), format = "f", digits = digits)))
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
    
    df <- read.csv(
      "/home/riccardo/codici_progetti/MSA_paper/MSA_long_paper/matrices/CT_cross_HC_vs_crossMSA.csv",
      stringsAsFactors = FALSE
    )
    
    savename <- "/home/riccardo/codici_progetti/MSA_paper/MSA_long_paper/matrices/fsquare_ancova_CT.csv"
    
  } else if (ATLAS == "Schaefer") {
    
    areas <- c(
      "VisCent","VisPeri","SomMotA","SomMotB","DorsAttnA","DorsAttnB","SalVentAttnA","SalVentAttnB",
      "LimbicA","LimbicB","ContA","ContB","ContC","DefaultA","DefaultB","DefaultC","TempPar"
    )
    
    df <- read.csv(
      paste0(
        "/home/riccardo/codici_progetti/MSA_paper/MSA_long_paper/matrices/CT_cross_HC_vscrossMSA_",
        ATLAS, "_Parc", NPARC, "_NET", NNET, ".csv"  ),
      stringsAsFactors = FALSE
    )
    savename <- paste0(
      "/home/riccardo/codici_progetti/MSA_paper/MSA_long_paper/matrices/fsquare_ancova_CT_",
      ATLAS, "_Parc", NPARC, "_NET", NNET, ".csv" )
    
  } else {
    stop("ATLAS must be 'DKT' or 'SCHAEFER' when doCT=TRUE.")
  }
  
} else if (!doCT && doVOL) {
  
  areas <- c("Putamen","CerebellumWhiteMatter","CerebellumCortex","Pons")
  
  df <- read.csv(
    "/home/riccardo/codici_progetti/MSA_paper/MSA_long_paper/matrices/VOL_cross_HC_vs_crossMSA.csv",
    stringsAsFactors = FALSE
  )
  
  savename <- "/home/riccardo/codici_progetti/MSA_paper/MSA_long_paper/matrices/fsquare_ancova_VOL.csv"
  
} else {
  
  stop("Choose either doCT=TRUE, doVOL=FALSE OR doCT=FALSE, doVOL=TRUE")
  
}

# ---- make sure factors are factors ----
df$GROUP <- as.factor(df$GROUP)
df$sex_final <- as.factor(df$sex_final)

options(contrasts = c("contr.treatment", "contr.poly"))
df$GROUP <- relevel(df$GROUP, ref = "HC")

nROI <- length(areas)

# ---- storage ----
CI_low   <- vector("list", nROI)
CI_high  <- vector("list", nROI)
R2_full  <- numeric(nROI)
R2_reduced <- numeric(nROI)
f2_GROUP_ancova <-numeric(nROI)

anova_tbl <- vector("list", nROI)
eta_tbl  <- vector("list", nROI)
eta_GROUP   <- numeric(nROI)
p_GROUP     <- numeric(nROI)



res_group <- data.frame(
  ROI = areas,
  Estimate = NA_real_,
  CI95_Lower = NA_real_,
  CI95_Upper = NA_real_,
  tStat = NA_real_,
  p_txt = NA_character_,
  f2_GROUP = NA_real_,
  stringsAsFactors = FALSE
)

for (roi in seq_len(nROI)) {
  y <- areas[roi]
  cat("\n=============================\nROI:", y, "\n")
  
  f <- as.formula(paste(y, "~ GROUP + age_final + sex_final"))
  ancova_fit <- lm(f, data = df)
  ancova_reduced_fit <- update(ancova_fit, . ~ . - GROUP)
  
  cs <- summary(ancova_fit)$coefficients
  grp_row <- grep("^GROUP", rownames(cs), value = TRUE)
  stopifnot(length(grp_row) == 1)
  est <- cs[grp_row, "Estimate"]
  tval <- cs[grp_row, "t value"]
  pval <- cs[grp_row, "Pr(>|t|)"]
  ci <- confint(ancova_fit)[grp_row, ]
  lo <- ci[1]; hi <- ci[2]
  
  # --- f2 from R2 drop ---
  R2_full <- summary(ancova_fit)$r.squared
  R2_red  <- summary(ancova_reduced_fit)$r.squared
  f2g <- (R2_full - R2_red) / (1 - R2_full)
  
  res_group[roi, ] <- list(y, est, lo, hi, tval, fmt_p(pval), f2g)
}

# round numeric columns (leave p_txt as text)
num_cols <- sapply(res_group, is.numeric)
res_group[num_cols] <- round(res_group[num_cols], 3)

write.csv(res_group, savename, row.names = FALSE)
print(res_group)
print(df$GROUP)

print(savename)
#eta_GROUP <- sapply(eta_tbl, function(x) {
# effectsize output has a "Parameter" column naming terms
#  x$Eta2_partial[x$Parameter == "GROUP"]
#})

#names(eta_GROUP) <- areas
#names(f2_GROUP_ancova) <- areas
#names(p_GROUP)   <- areas

#print(p_GROUP)
#print(f2_GROUP_ancova)

#df_f2 <- data.frame(
#  ROI = names(f2_GROUP_ancova),
#  Value = as.numeric(f2_GROUP_ancova)
#)
#df_f2$ROI <- factor(df_f2$ROI, levels = df_f2$ROI[order(df_f2$Value)])

#write.csv(df_f2, savename, row.names = FALSE)
#p<-ggplot(df_f2, aes(x = Value, y = ROI)) +
#  geom_col(fill = "gray60") +
#  labs(x = "Effect sizes", y = NULL) +
#  theme_minimal(base_size = 14) +
#  theme(
#    axis.text.y = element_text(size = 12),
#    panel.grid.major.y = element_blank()
#  )

#print(p)