Expected input files (simple, rectangular .csv; one row per SampleID)

data/
  labels.csv                # SampleID, CES_score, PBS_score, PopCS_score, CES_group(high/low), ZIP(or Tract)
  microbiome.csv            # SampleID, <species/taxa...>  (already variance-stabilized or CLR/log CPM)
  brain_fc.csv              # SampleID, <edge_1 ... edge_P> (Fisher z)
  cytokines.csv             # SampleID, <GMCSF, IL17F, IL13, TNFA, VEGFA, ...>
  brainstem_volumes.csv     # SampleID, <AAN_LC, AAN_PAG, AAN_VTA, AAN_PnO, ...> (mm^3)
  scfa_hormones.csv         # SampleID, <acetate, propionate, butyrate, lactate, GLP1, PYY, insulin, Cpeptide, ghrelin>

# ---------------------------
# Minimal code to replicate core findings
# by Dr. Gabriel A. Vignolle
# UCLA Health
# ---------------------------
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr)
  library(mixOmics)     # sPLS-DA / DIABLO
  library(pROC)         # ROC/AUC
})

# ---------- 1) Load data ----------
lbl   <- read_csv("data/labels.csv")          # SampleID, CES_score, PBS_score, PopCS_score, CES_group, ZIP
mb    <- read_csv("data/microbiome.csv")      # SampleID + taxa
fc    <- read_csv("data/brain_fc.csv")        # SampleID + edges
cyt   <- read_csv("data/cytokines.csv")       # SampleID + cytokines
aan   <- read_csv("data/brainstem_volumes.csv") # SampleID + AAN nuclei
sh    <- read_csv("data/scfa_hormones.csv")   # SampleID + SCFAs + hormones

# Harmonize IDs
all_ids <- Reduce(intersect, list(lbl$SampleID, mb$SampleID, fc$SampleID))
lbl <- lbl %>% filter(SampleID %in% all_ids)
mb  <- mb  %>% filter(SampleID %in% all_ids) %>% arrange(SampleID)
fc  <- fc  %>% filter(SampleID %in% all_ids) %>% arrange(SampleID)

# Response for integration
Y <- factor(lbl$CES_group[match(mb$SampleID, lbl$SampleID)], levels = c("low_CES","high_CES"))

# Drop ID columns and coerce to matrices
X_micro <- mb %>% select(-SampleID) %>% as.matrix()
X_brain <- fc %>% select(-SampleID) %>% as.matrix()

# ---------- 2) Dose–response (univariate) ----------
# 2a) PBS ~ cytokines (expect predominantly negative slopes)
message("Dose-response: PBS ~ cytokines")
xy <- cyt %>% inner_join(lbl, by = "SampleID")
res_cyt <- lapply(names(cyt)[-1], function(k) {
  fit <- lm(scale(xy[[k]]) ~ scale(xy$PBS_score))  # age/sex/BMI can be added if present
  c(feature = k,
    beta   = coef(fit)[2],
    pval   = summary(fit)$coefficients[2,4])
}) %>% bind_rows()
res_cyt$FDR <- p.adjust(res_cyt$pval, "BH")
print(res_cyt %>% arrange(pval) %>% head(10))

# 2b) CES/PopCS ~ AAN volumes (expect negative slopes)
message("Dose-response: CES/PopCS ~ AAN volumes")
bv <- aan %>% inner_join(lbl, by = "SampleID")
aan_cols <- setdiff(names(aan), "SampleID")
fit_AAN <- lapply(aan_cols, function(k) {
  fit1 <- lm(scale(bv[[k]]) ~ scale(bv$CES_score))
  fit2 <- lm(scale(bv[[k]]) ~ scale(bv$PopCS_score))
  tibble(region = k,
         beta_CES = coef(fit1)[2], p_CES = summary(fit1)$coefficients[2,4],
         beta_Pop = coef(fit2)[2], p_Pop = summary(fit2)$coefficients[2,4])
}) %>% bind_rows()
fit_AAN <- fit_AAN %>% mutate(FDR_CES = p.adjust(p_CES, "BH"),
                              FDR_Pop  = p.adjust(p_Pop,  "BH"))
print(fit_AAN %>% arrange(p_CES) %>% head(10))

# ---------- 3) sPLS-DA (microbiome + brain FC), nested CV with spatial group-K-fold ----------
# Helper: create spatial group K-fold indices (by ZIP/tract)
make_group_folds <- function(groups, K = 5, seed = 1) {
  set.seed(seed)
  glev <- sample(unique(groups))
  # split unique groups into K folds
  split(glev, cut(seq_along(glev), breaks = K, labels = FALSE))
}
# Inner tuner for keepX per block
tune_inner <- function(Xlist, Y, groups, K_outer = 5, K_inner = 4,
                       keepX_grid_micro = c(30, 60, 90),
                       keepX_grid_brain = c(30, 60, 90),
                       ncomp = 3, seed = 1) {
  set.seed(seed)
  outer_folds <- make_group_folds(groups, K = K_outer, seed = seed)
  auc_outer <- c()
  sel_keepX <- list()
  for (f in seq_along(outer_folds)) {
    test_groups <- outer_folds[[f]]
    test_idx <- which(groups %in% test_groups)
    train_idx <- setdiff(seq_along(Y), test_idx)

    Xtr <- lapply(Xlist, function(x) x[train_idx, , drop = FALSE])
    Ytr <- Y[train_idx]

    # mixOmics tuner (inner CV, Mfold by group labels to avoid leakage)
    tune <- tune.block.splsda(
      X = Xtr, Y = Ytr, ncomp = ncomp,
      test.keepX = list(microbiome = keepX_grid_micro,
                        brain      = keepX_grid_brain),
      design = matrix(0.1, nrow = 2, ncol = 2, dimnames = list(names(Xlist), names(Xlist))),
      validation = 'Mfold',
      folds = K_inner,  # inner-CV folds
      progressBar = FALSE, dist = "max.dist"
    )
    keepX <- tune$choice.keepX
    sel_keepX[[f]] <- keepX

    # Fit on train, predict test
    final <- block.splsda(X = Xtr, Y = Ytr, ncomp = ncomp, keepX = keepX,
                          design = matrix(0.1, 2, 2, dimnames = list(names(Xlist), names(Xlist))))
    Xte <- lapply(Xlist, function(x) x[test_idx, , drop = FALSE])
    pred <- predict(final, Xte)
    # Use majority vote; also compute block‑wise linear predictor for ROC
    prob_micro <- pred$predict$microbiome[, , 1][, "high_CES", drop = TRUE]
    prob_brain <- pred$predict$brain[, , 1][, "high_CES", drop = TRUE]
    yte <- Y[test_idx]
    auc_micro <- pROC::roc(yte, prob_micro, quiet = TRUE)$auc
    auc_brain <- pROC::roc(yte, prob_brain, quiet = TRUE)$auc
    auc_outer <- rbind(auc_outer, c(microbiome = auc_micro, brain = auc_brain))
  }
  list(outer_auc = as.data.frame(auc_outer),
       selected_keepX = sel_keepX)
}

# Assemble blocks and groups
ord <- match(mb$SampleID, lbl$SampleID)
groups_zip <- lbl$ZIP[ord]
Xlist <- list(microbiome = X_micro, brain = X_brain)

tuned <- tune_inner(Xlist, Y, groups = groups_zip,
                    K_outer = 5, K_inner = 4,
                    keepX_grid_micro = c(50, 100, 150),
                    keepX_grid_brain = c(30, 50, 70),
                    ncomp = 3, seed = 123)

print(tuned$outer_auc)
cat("Mean AUC (microbiome):", round(mean(tuned$outer_auc$microbiome), 3), "\n")
cat("Mean AUC (brain):",      round(mean(tuned$outer_auc$brain), 3), "\n")

# Fit final model on full data with median keepX across outer folds
median_keep <- function(lst, block) {
  vals <- sapply(lst, function(kx) kx[[block]])
  # if values are vectors per comp, take round of colwise median
  round(apply(vals, 1, median))
}
keepX_final <- list(
  microbiome = median_keep(tuned$selected_keepX, "microbiome"),
  brain      = median_keep(tuned$selected_keepX, "brain")
)
final_model <- block.splsda(X = Xlist, Y = Y, ncomp = 3,
                            keepX = keepX_final,
                            design = matrix(0.1, 2, 2, dimnames = list(names(Xlist), names(Xlist))))

# Block‑wise ROC on full data (resubstitution; report for reference)
pp <- predict(final_model, Xlist)
prob_m <- pp$predict$microbiome[, , 1][, "high_CES"]
prob_b <- pp$predict$brain[, , 1][, "high_CES"]
cat("Resubstitution AUC (microbiome):", round(pROC::roc(Y, prob_m, quiet = TRUE)$auc, 4), "\n")
cat("Resubstitution AUC (brain):",      round(pROC::roc(Y, prob_b, quiet = TRUE)$auc, 4), "\n")

# Export selected features
sel_micro <- selectVar(final_model, block = "microbiome", comp = 1)$name
sel_brain <- selectVar(final_model, block = "brain",       comp = 1)$name
writeLines(sel_micro, "outputs/selected_microbiome_comp1.txt")
writeLines(sel_brain, "outputs/selected_brain_edges_comp1.txt")

# ---------- 4) SCFA ↔ enteroendocrine hub screen (Spearman) ----------
message("SCFA–enteroendocrine hub screen")
sh2 <- sh %>% inner_join(lbl %>% select(SampleID), by = "SampleID")
scfas   <- c("acetate","propionate","butyrate","lactate")
horms   <- c("GLP1","PYY","insulin","Cpeptide","ghrelin")
M <- sh2 %>% select(all_of(scfas), all_of(horms))
cors <- cor(M[, scfas], M[, horms], method = "spearman", use = "pairwise.complete.obs")
hub <- which(abs(cors) >= 0.55, arr.ind = TRUE)  # threshold as in manuscript
hub_df <- tibble(SCFA = rownames(cors)[hub[,1]],
                 Hormone = colnames(cors)[hub[,2]],
                 rho = cors[hub])
print(hub_df %>% arrange(desc(abs(rho))))
write.csv(hub_df, "outputs/scfa_hormone_hubs.csv", row.names = FALSE)

# ---------- 5) Session info ----------
dir.create("outputs", showWarnings = FALSE)
sink("outputs/sessionInfo.txt"); print(sessionInfo()); sink()
