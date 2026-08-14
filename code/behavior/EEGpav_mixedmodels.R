rm(list = ls())
clc <- function() cat("\014")
clc()

rootdir <- "H:/mt_behavalluse"
setwd(rootdir)

suppressPackageStartupMessages({
  library(lme4)
  library(lmerTest)
  library(car)
  library(emmeans)
  library(dplyr)
  library(ggplot2)
  library(writexl)
})

outdir <- file.path(rootdir, "behavmixmodel")
figdir <- file.path(rootdir, "group主效应")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
if (!dir.exists(figdir)) dir.create(figdir, recursive = TRUE)

rawdata <- read.csv("EEG4mixedmodelR_withGroup.csv", header = TRUE)


# ============================================================
# 1. Data preparation
# ============================================================

prepare_data <- function(dat){

  dat$valence[dat$valence == 0] <- -1
  dat$action[dat$action == 0] <- -1

  dat$accuracy <- dat$DVacc
  dat$accuracy[dat$accuracy == 0] <- -1

  dat$DVcorrectGo <- dat$DVgo * dat$DVacc
  dat$DVincorrectGo <- dat$DVgo * (1 - dat$DVacc)

  dat$DVrt[dat$DVrt < 0.1] <- NA
  dat$DVlnrt <- log(dat$DVrt + 0.9)

  dat
}

data_all <- prepare_data(rawdata)


# ============================================================
# 2. Helper functions
# ============================================================

extract_fixed <- function(mod, model_name){

  sm <- coef(summary(mod))

  if ("z value" %in% colnames(sm)) {

    data.frame(
      model = model_name,
      term = rownames(sm),
      beta = sm[, "Estimate"],
      se = sm[, "Std. Error"],
      df = NA,
      statistic = sm[, "z value"],
      statistic_type = "z",
      p = sm[, "Pr(>|z|)"],
      row.names = NULL
    )

  } else {

    data.frame(
      model = model_name,
      term = rownames(sm),
      beta = sm[, "Estimate"],
      se = sm[, "Std. Error"],
      df = if ("df" %in% colnames(sm)) sm[, "df"] else NA,
      statistic = sm[, "t value"],
      statistic_type = "t",
      p = if ("Pr(>|t|)" %in% colnames(sm)) sm[, "Pr(>|t|)"] else NA,
      row.names = NULL
    )
  }
}


extract_anova <- function(mod, model_name){

  a <- as.data.frame(car::Anova(mod, type = 2))
  a$effect <- rownames(a)
  rownames(a) <- NULL
  a$model <- model_name

  a %>% select(model, effect, everything())
}


extract_model_results <- function(models){

  fixed_list <- list()
  anova_list <- list()

  for (nm in names(models)) {
    fixed_list[[nm]] <- extract_fixed(models[[nm]], nm)
    anova_list[[nm]] <- extract_anova(models[[nm]], nm)
  }

  list(
    fixed = bind_rows(fixed_list),
    anova = bind_rows(anova_list)
  )
}


check_models <- function(models){

  data.frame(
    model = names(models),
    singular = sapply(models, function(x) isSingular(x, tol = 1e-4)),
    row.names = NULL
  )
}


fmt3 <- function(x){
  ifelse(is.na(x), NA, formatC(x, digits = 3, format = "f"))
}


p_to_stars <- function(p){
  ifelse(
    p < 0.001, "***",
    ifelse(
      p < 0.01, "**",
      ifelse(p < 0.05, "*", "ns")
    )
  )
}


# ============================================================
# 3. HC / OCD separate models
# ============================================================

fit_group_models <- function(dat, group_value, group_name){

  d <- dat[dat$group == group_value, ]

  cat("\n================ ", group_name, " ================\n")
  cat("Subjects:", length(unique(d$sID)), "\n")
  cat("Trials:", nrow(d), "\n")
  cat("Valid RT trials:", sum(!is.na(d$DVlnrt)), "\n")

  Mgo.basic <- glmer(
    DVgo ~ valence * action + (1 + valence * action | sID),
    family = binomial,
    data = d,
    control = glmerControl(
      optimizer = "bobyqa",
      optCtrl = list(maxfun = 10000)
    )
  )

  Mcorrectgo <- glmer(
    DVcorrectGo ~ valence + (1 + valence | sID),
    family = binomial,
    data = d[d$DVincorrectGo != 1 & d$action == 1, ],
    control = glmerControl(
      optimizer = "bobyqa",
      optCtrl = list(maxfun = 10000)
    )
  )

  Mincorrectgo <- glmer(
    DVincorrectGo ~ valence + (1 + valence | sID),
    family = binomial,
    data = d[d$DVcorrectGo != 1 & d$action == 1, ],
    control = glmerControl(
      optimizer = "bobyqa",
      optCtrl = list(maxfun = 10000)
    )
  )

  Maccuracyofgo <- glmer(
    DVacc ~ valence + (1 + valence | sID),
    family = binomial,
    data = d[d$DVgo == 1 & d$action == 1, ],
    control = glmerControl(
      optimizer = "bobyqa",
      optCtrl = list(maxfun = 10000)
    )
  )

  Mrt.basic <- lmerTest::lmer(
    DVlnrt ~ valence * action + (1 + valence * action | sID),
    data = d[!is.na(d$DVlnrt), ],
    control = lmerControl(
      optimizer = "bobyqa",
      optCtrl = list(maxfun = 10000)
    )
  )

  Mrt.accuracy <- lmerTest::lmer(
    DVlnrt ~ valence * accuracy + (1 + valence * accuracy | sID),
    data = d[d$action == 1 & !is.na(d$DVlnrt), ],
    control = lmerControl(
      optimizer = "bobyqa",
      optCtrl = list(maxfun = 10000)
    )
  )

  models <- list(
    Mgo_basic = Mgo.basic,
    McorrectGo = Mcorrectgo,
    MincorrectGo = Mincorrectgo,
    MaccuracyOfGo = Maccuracyofgo,
    Mrt_basic = Mrt.basic,
    Mrt_accuracy = Mrt.accuracy
  )

  results <- extract_model_results(models)
  singular <- check_models(models)

  emm_pgo <- as.data.frame(
    emmeans(
      Mgo.basic,
      ~ valence * action,
      at = list(
        valence = c(-1, 1),
        action = c(-1, 1)
      ),
      type = "response"
    )
  )

  emm_valence_by_action <- as.data.frame(
    pairs(
      emmeans(
        Mgo.basic,
        ~ valence | action,
        at = list(
          valence = c(-1, 1),
          action = c(-1, 1)
        )
      )
    )
  )

  emm_action_by_valence <- as.data.frame(
    pairs(
      emmeans(
        Mgo.basic,
        ~ action | valence,
        at = list(
          valence = c(-1, 1),
          action = c(-1, 1)
        )
      )
    )
  )

  write_xlsx(
    list(
      Fixed_effects = results$fixed,
      TypeII_tests = results$anova,
      Singularity = singular,
      PGo_EMMeans = emm_pgo,
      PGo_valence_by_action = emm_valence_by_action,
      PGo_action_by_valence = emm_action_by_valence
    ),
    path = file.path(
      outdir,
      paste0("EEG_", group_name, "_mixedmodel_effects.xlsx")
    )
  )

  ranefGo <- ranef(Mgo.basic)$sID
  ranefRT <- ranef(Mrt.basic)$sID

  cor_pearson <- cor.test(
    ranefGo$valence,
    ranefRT$valence,
    method = "pearson"
  )

  cor_spearman <- cor.test(
    ranefGo$valence,
    ranefRT$valence,
    method = "spearman",
    exact = FALSE
  )

  cor_table <- data.frame(
    method = c("Pearson", "Spearman"),
    estimate = c(
      unname(cor_pearson$estimate),
      unname(cor_spearman$estimate)
    ),
    statistic = c(
      unname(cor_pearson$statistic),
      unname(cor_spearman$statistic)
    ),
    p = c(
      cor_pearson$p.value,
      cor_spearman$p.value
    )
  )

  write.csv(
    cor_table,
    file.path(
      outdir,
      paste0(group_name, "_PGo_RT_random_slope_correlations.csv")
    ),
    row.names = FALSE
  )

  ranef_df <- data.frame(
    PGo_valence = ranefGo$valence,
    RT_valence = ranefRT$valence
  )

  p_cor <- ggplot(
    ranef_df,
    aes(x = PGo_valence, y = RT_valence)
  ) +
    geom_point(size = 2) +
    geom_smooth(
      method = "lm",
      se = FALSE,
      linewidth = 0.7
    ) +
    labs(
      x = "P(Go) valence random slope",
      y = "RT valence random slope"
    ) +
    theme_classic(base_size = 12)

  ggsave(
    file.path(
      outdir,
      paste0("ranefpGo_RT_", group_name, ".png")
    ),
    p_cor,
    width = 4,
    height = 4,
    dpi = 600
  )

  saveRDS(
    list(
      data = d,
      models = models
    ),
    file.path(
      outdir,
      paste0("EEGpav_mm_", group_name, ".rds")
    )
  )

  invisible(
    list(
      data = d,
      models = models,
      results = results
    )
  )
}


HC <- fit_group_models(
  data_all,
  group_value = 0,
  group_name = "HC"
)

OCD <- fit_group_models(
  data_all,
  group_value = 1,
  group_name = "OCD"
)


# ============================================================
# 4. Joint HC-OCD models
# ============================================================

data_joint <- data_all

data_joint$group <- factor(
  data_joint$group,
  levels = c(0, 1),
  labels = c("HC", "OCD")
)

Mgo.joint <- glmer(
  DVgo ~ group * valence * action +
    (1 + valence * action | sID),
  family = binomial,
  data = data_joint,
  control = glmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 10000)
  )
)

Mcorrectgo.joint <- glmer(
  DVcorrectGo ~ group * valence +
    (1 + valence | sID),
  family = binomial,
  data = data_joint[
    data_joint$DVincorrectGo != 1 &
      data_joint$action == 1,
  ],
  control = glmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 10000)
  )
)

Mincorrectgo.joint <- glmer(
  DVincorrectGo ~ group * valence +
    (1 + valence | sID),
  family = binomial,
  data = data_joint[
    data_joint$DVcorrectGo != 1 &
      data_joint$action == 1,
  ],
  control = glmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 10000)
  )
)

Maccuracyofgo.joint <- glmer(
  DVacc ~ group * valence +
    (1 + valence | sID),
  family = binomial,
  data = data_joint[
    data_joint$DVgo == 1 &
      data_joint$action == 1,
  ],
  control = glmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 10000)
  )
)

Mrt.basic.joint <- lmerTest::lmer(
  DVlnrt ~ group * valence * action +
    (1 + valence * action | sID),
  data = data_joint[
    !is.na(data_joint$DVlnrt),
  ],
  control = lmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 10000)
  )
)

Mrt.accuracy.joint <- lmerTest::lmer(
  DVlnrt ~ group * valence * accuracy +
    (1 + valence * accuracy | sID),
  data = data_joint[
    data_joint$action == 1 &
      !is.na(data_joint$DVlnrt),
  ],
  control = lmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 10000)
  )
)

joint_models <- list(
  Mgo_joint = Mgo.joint,
  Mcorrectgo_joint = Mcorrectgo.joint,
  Mincorrectgo_joint = Mincorrectgo.joint,
  Maccuracyofgo_joint = Maccuracyofgo.joint,
  Mrt_basic_joint = Mrt.basic.joint,
  Mrt_accuracy_joint = Mrt.accuracy.joint
)

joint_results <- extract_model_results(joint_models)
joint_singular <- check_models(joint_models)

write.csv(
  joint_results$fixed,
  file.path(
    outdir,
    "EEGpav_mm_group_fixed_effects.csv"
  ),
  row.names = FALSE
)

write.csv(
  joint_results$anova,
  file.path(
    outdir,
    "EEGpav_mm_group_Anova.csv"
  ),
  row.names = FALSE
)

write_xlsx(
  list(
    Fixed_effects = joint_results$fixed,
    TypeII_tests = joint_results$anova,
    Singularity = joint_singular
  ),
  path = file.path(
    outdir,
    "EEGpav_mm_group_results.xlsx"
  )
)

save(
  data_joint,
  Mgo.joint,
  Mcorrectgo.joint,
  Mincorrectgo.joint,
  Maccuracyofgo.joint,
  Mrt.basic.joint,
  Mrt.accuracy.joint,
  file = file.path(
    rootdir,
    "EEGpav_group.RData"
  )
)


# ============================================================
# 5. Joint-model P(Go): Group × Valence × Action
# ============================================================

emm_pgo_joint <- emmeans(
  Mgo.joint,
  ~ group * valence * action,
  at = list(
    valence = c(-1, 1),
    action = c(-1, 1)
  ),
  type = "response"
)

emm_pgo_joint_df <- as.data.frame(emm_pgo_joint) %>%
  mutate(
    Valence = factor(
      valence,
      levels = c(-1, 1),
      labels = c("Aversive", "Appetitive")
    ),
    Action = factor(
      action,
      levels = c(-1, 1),
      labels = c("NoGo", "Go")
    )
  )

write.csv(
  emm_pgo_joint_df,
  file.path(
    figdir,
    "Mgo_joint_groupXvalenceXaction_emmeans.csv"
  ),
  row.names = FALSE
)


emm_pgo_group_contrast <- pairs(
  emmeans(
    Mgo.joint,
    ~ group | valence * action,
    at = list(
      valence = c(-1, 1),
      action = c(-1, 1)
    )
  )
)

write.csv(
  as.data.frame(emm_pgo_group_contrast),
  file.path(
    figdir,
    "Mgo_joint_HC_vs_OCD_by_valence_action.csv"
  ),
  row.names = FALSE
)


emm_pgo_valence_by_group_action <- pairs(
  emmeans(
    Mgo.joint,
    ~ valence | group * action,
    at = list(
      valence = c(-1, 1),
      action = c(-1, 1)
    )
  )
)

write.csv(
  as.data.frame(emm_pgo_valence_by_group_action),
  file.path(
    figdir,
    "Mgo_joint_valence_effect_by_group_action.csv"
  ),
  row.names = FALSE
)


p_pgo <- ggplot(
  emm_pgo_joint_df,
  aes(
    x = Valence,
    y = prob,
    group = group,
    colour = group,
    shape = group
  )
) +
  geom_errorbar(
    aes(
      ymin = asymp.LCL,
      ymax = asymp.UCL
    ),
    position = position_dodge(width = 0.15),
    width = 0.08,
    linewidth = 0.6
  ) +
  geom_line(
    position = position_dodge(width = 0.15),
    linewidth = 0.8
  ) +
  geom_point(
    position = position_dodge(width = 0.15),
    size = 2.6
  ) +
  facet_wrap(~Action) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2)
  ) +
  labs(
    x = "Valence",
    y = "P(Go)",
    colour = NULL,
    shape = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "top"
  )

ggsave(
  file.path(
    figdir,
    "Mgo_joint_groupXvalenceXaction_pGo.png"
  ),
  p_pgo,
  width = 6.5,
  height = 4,
  dpi = 600
)


# ============================================================
# 6. P(correct Go): Group
# ============================================================

emm_correct <- emmeans(
  Mcorrectgo.joint,
  ~ group,
  at = list(
    valence = c(-1, 1)
  ),
  weights = "equal",
  type = "response"
)

emm_correct_df <- as.data.frame(emm_correct)

contrast_correct <- pairs(
  emmeans(
    Mcorrectgo.joint,
    ~ group,
    at = list(
      valence = c(-1, 1)
    ),
    weights = "equal"
  )
)

write.csv(
  emm_correct_df,
  file.path(
    figdir,
    "Mcorrectgo_joint_group_emmeans.csv"
  ),
  row.names = FALSE
)

write.csv(
  as.data.frame(contrast_correct),
  file.path(
    figdir,
    "Mcorrectgo_joint_HC_vs_OCD.csv"
  ),
  row.names = FALSE
)


p_correct <- ggplot(
  emm_correct_df,
  aes(
    x = group,
    y = prob,
    colour = group,
    shape = group
  )
) +
  geom_errorbar(
    aes(
      ymin = asymp.LCL,
      ymax = asymp.UCL
    ),
    width = 0.08,
    linewidth = 0.6
  ) +
  geom_point(size = 2.8) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2)
  ) +
  labs(
    x = "Group",
    y = "P(correct Go)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none"
  )

ggsave(
  file.path(
    figdir,
    "Mcorrectgo_joint_group.png"
  ),
  p_correct,
  width = 4,
  height = 4,
  dpi = 600
)


# ============================================================
# 7. P(correct | Go): Group
# ============================================================

emm_acc <- emmeans(
  Maccuracyofgo.joint,
  ~ group,
  at = list(
    valence = c(-1, 1)
  ),
  weights = "equal",
  type = "response"
)

emm_acc_df <- as.data.frame(emm_acc)

contrast_acc <- pairs(
  emmeans(
    Maccuracyofgo.joint,
    ~ group,
    at = list(
      valence = c(-1, 1)
    ),
    weights = "equal"
  )
)

write.csv(
  emm_acc_df,
  file.path(
    figdir,
    "Maccuracyofgo_joint_group_emmeans.csv"
  ),
  row.names = FALSE
)

write.csv(
  as.data.frame(contrast_acc),
  file.path(
    figdir,
    "Maccuracyofgo_joint_HC_vs_OCD.csv"
  ),
  row.names = FALSE
)


p_acc <- ggplot(
  emm_acc_df,
  aes(
    x = group,
    y = prob,
    colour = group,
    shape = group
  )
) +
  geom_errorbar(
    aes(
      ymin = asymp.LCL,
      ymax = asymp.UCL
    ),
    width = 0.08,
    linewidth = 0.6
  ) +
  geom_point(size = 2.8) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2)
  ) +
  labs(
    x = "Group",
    y = "P(correct | Go)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none"
  )

ggsave(
  file.path(
    figdir,
    "Maccuracyofgo_joint_group.png"
  ),
  p_acc,
  width = 4,
  height = 4,
  dpi = 600
)


# ============================================================
# 8. RT: Group × Action
# ============================================================

emm_rt <- emmeans(
  Mrt.basic.joint,
  ~ group * action,
  at = list(
    valence = c(-1, 1),
    action = c(-1, 1)
  ),
  weights = "equal"
)

emm_rt_df <- as.data.frame(
  confint(
    emm_rt,
    level = 0.95
  )
) %>%
  mutate(
    rt_sec = exp(emmean) - 0.9,
    rt_sec_low = exp(lower.CL) - 0.9,
    rt_sec_high = exp(upper.CL) - 0.9,
    Action = factor(
      action,
      levels = c(-1, 1),
      labels = c("NoGo", "Go")
    )
  )

contrast_rt <- pairs(
  emmeans(
    Mrt.basic.joint,
    ~ group | action,
    at = list(
      valence = c(-1, 1),
      action = c(-1, 1)
    ),
    weights = "equal"
  )
)

write.csv(
  emm_rt_df,
  file.path(
    figdir,
    "Mrt_basic_joint_groupXaction_emmeans.csv"
  ),
  row.names = FALSE
)

write.csv(
  as.data.frame(contrast_rt),
  file.path(
    figdir,
    "Mrt_basic_joint_HC_vs_OCD_by_action.csv"
  ),
  row.names = FALSE
)


p_rt <- ggplot(
  emm_rt_df,
  aes(
    x = Action,
    y = rt_sec,
    group = group,
    colour = group,
    shape = group
  )
) +
  geom_errorbar(
    aes(
      ymin = rt_sec_low,
      ymax = rt_sec_high
    ),
    position = position_dodge(width = 0.15),
    width = 0.08,
    linewidth = 0.6
  ) +
  geom_line(
    position = position_dodge(width = 0.15),
    linewidth = 0.8
  ) +
  geom_point(
    position = position_dodge(width = 0.15),
    size = 2.6
  ) +
  labs(
    x = "Action",
    y = "Model-estimated RT (s)",
    colour = NULL,
    shape = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "top"
  )

ggsave(
  file.path(
    figdir,
    "Mrt_basic_joint_groupXaction.png"
  ),
  p_rt,
  width = 5.2,
  height = 4,
  dpi = 600
)


# ============================================================
# 9. RT: Group × Accuracy
# ============================================================

emm_rt_acc <- emmeans(
  Mrt.accuracy.joint,
  ~ group * accuracy,
  at = list(
    valence = c(-1, 1),
    accuracy = c(-1, 1)
  ),
  weights = "equal"
)

emm_rt_acc_df <- as.data.frame(
  confint(
    emm_rt_acc,
    level = 0.95
  )
) %>%
  mutate(
    rt_sec = exp(emmean) - 0.9,
    rt_sec_low = exp(lower.CL) - 0.9,
    rt_sec_high = exp(upper.CL) - 0.9,
    Accuracy = factor(
      accuracy,
      levels = c(-1, 1),
      labels = c("Error", "Correct")
    )
  )

contrast_rt_acc <- pairs(
  emmeans(
    Mrt.accuracy.joint,
    ~ group | accuracy,
    at = list(
      valence = c(-1, 1),
      accuracy = c(-1, 1)
    ),
    weights = "equal"
  )
)

write.csv(
  emm_rt_acc_df,
  file.path(
    figdir,
    "Mrt_accuracy_joint_groupXaccuracy_emmeans.csv"
  ),
  row.names = FALSE
)

write.csv(
  as.data.frame(contrast_rt_acc),
  file.path(
    figdir,
    "Mrt_accuracy_joint_HC_vs_OCD_by_accuracy.csv"
  ),
  row.names = FALSE
)


p_rt_acc <- ggplot(
  emm_rt_acc_df,
  aes(
    x = Accuracy,
    y = rt_sec,
    group = group,
    colour = group,
    shape = group
  )
) +
  geom_errorbar(
    aes(
      ymin = rt_sec_low,
      ymax = rt_sec_high
    ),
    position = position_dodge(width = 0.15),
    width = 0.08,
    linewidth = 0.6
  ) +
  geom_line(
    position = position_dodge(width = 0.15),
    linewidth = 0.8
  ) +
  geom_point(
    position = position_dodge(width = 0.15),
    size = 2.6
  ) +
  labs(
    x = "Accuracy",
    y = "Model-estimated RT (s)",
    colour = NULL,
    shape = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "top"
  )

ggsave(
  file.path(
    figdir,
    "Mrt_accuracy_joint_groupXaccuracy.png"
  ),
  p_rt_acc,
  width = 5.2,
  height = 4,
  dpi = 600
)


# ============================================================
# 10. Final checks
# ============================================================

cat("\nAnalysis completed.\n")
cat("HC subjects:", length(unique(HC$data$sID)), "\n")
cat("OCD subjects:", length(unique(OCD$data$sID)), "\n")
cat("Main results:", outdir, "\n")
cat("Group figures:", figdir, "\n")
