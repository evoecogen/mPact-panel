####
# mPAct -
# Pseudomonas aeruginosa major clone type panel
#
# Quantitative phenotype correlations
####
# R 4.3.1
# tidyverse 2.0
####



library(tidyverse)
library(ggpubr)
library(broom)
library(svglite)


# Full drug name labeller
drug.full.names <- c(
  'PIP' = 'piperacillin',
  'PIT' = 'piperacillin/tazobactam',
  'CTZ' = 'ceftazidime',
  'CTV' = 'ceftazidime/avibactam',
  'CEP' = 'cefepime',
  'CTT' = 'ceftolozane/tazobactam',
  'CID' = 'cefiderocol',
  'IMI' = 'imipenem',
  'MER' = 'meropenem',
  'AZT' = 'aztreonam',
  'CIP' = 'ciprofloxacin',
  'AMI' = 'amikacin',
  'GEN' = 'gentamicin',
  'TOB' = 'tobramycin',
  'COL' = 'colistin')

# Plot exporting wrapper function
# Same size png and svg of the specified plot. Defaults to golden ratio of height to the input
# width. Path as string with path, but file ending omitted.
print_svgpng <- function(plot, path = "file", width, height = width * 0.618) {
  ggsave(plot, filename = paste0(path, ".svg"), width = width, height = height, units = "in",
         device = svglite, fix_text_size = FALSE)
  ggsave(plot, filename = paste0(path, ".png"), width = width, height = height, units = "in")
}



##### 1) Load data -----
# Strain info
info <- read_csv2(file = "./data/basic_strain_info.csv")
# Growth parameters
growth <- read_csv2(file = "./tables/growth_parameters.csv")
# Disc diffusion
res.dd <- read_csv(file = "./data/disc_diffusion_m9.csv")
# Etest
res.etest <- read_csv(file = "./data/mic_etest_m9.csv")
# Automated sensitivity testing (AST)
res.ast <- read_csv(file = "./data/mic_vitek.csv")



##### 2) Format into singular data frame -----
# Start with basic strain info
cors <- info %>%
  select(strain, phylogroup, exo, plasmids, ICEs, origin_gen, origin_clin, origin_env)
# Add growth parameter means
cors <- growth %>%
  select(strain, growthrate_mean, lagphase_mean, auc_mean, maximumod_mean) %>%
  right_join(., cors, by = "strain")
# Add disc diffusion
cors <- res.dd %>%
  group_by(strain, drug) %>%
  summarize(ddmean = round(mean(diameter, na.rm = T), 2)) %>%
  pivot_wider(names_from = "drug", values_from = "ddmean", names_prefix = "dd.") %>%
  right_join(., cors, by = "strain")
# Add Etest data
cors <- res.etest %>%
  filter(microcol == "n" & time == "24h") %>%
  group_by(strain, drug) %>%
  summarize(micmean = round(mean(mic, na.rm = T), 3)) %>%
  select(strain, drug, micmean) %>%
  pivot_wider(names_from = "drug", values_from = "micmean", names_prefix = "mic.") %>%
  right_join(., cors, by = "strain")
# Add AST data, categorized
cors <- res.ast %>%
  pivot_wider(names_from = "drug", values_from = "mic", names_prefix = "ast.") %>%
  mutate(across(starts_with("ast"), 
                ~ parse_factor(as.character(.), 
                               levels = as.character(sort(unique(res.ast$mic)))))) %>%
  right_join(., cors, by = "strain")

# Tidy the assembled master data frame
cors.master <- cors %>%
  arrange(strain) %>%
  rename(growthrate = growthrate_mean, 
         maximumod = maximumod_mean, 
         lagphase = lagphase_mean, 
         auc = auc_mean) %>%
  select(strain, origin_gen, origin_clin, origin_env, 
         phylogroup, exo, plasmids, ICEs,
         growthrate, lagphase, maximumod, auc,
         dd.AMI, dd.CID, dd.CIP, dd.COL, dd.CTT, dd.CTV, 
         dd.CTZ, dd.GEN, dd.IMI, dd.MER, dd.PIT, dd.TOB,
         mic.AMI, mic.CIP, mic.COL, mic.CTT, mic.CTZ, mic.GEN, mic.MER, mic.PIT, mic.TOB,
         ast.PIP, ast.PIT, ast.CTZ, ast.CEP, ast.AZT, ast.IMI, 
         ast.MER, ast.AMI, ast.GEN, ast.TOB, ast.CIP, ast.COL) %>%
  mutate(across(c(plasmids, ICEs), ~ as.numeric(.x)))



##### 3) Cumulative resistance parameters

# a) Sum of disc diffusion diameters
cumu.res.dd <- cors.master %>%
  select(c(strain, starts_with("dd"))) %>%
  pivot_longer(!strain, names_to = "drug", values_to = "diameter", names_prefix = "dd.") %>%
  group_by(strain) %>%
  summarize(cr.ddd = sum(diameter))

# b) Sum of gradient strip MICs
cumu.res.mic <- cors.master %>%
  select(c(strain, starts_with("mic"))) %>%
  pivot_longer(!strain, names_to = "drug", values_to = "mic", names_prefix = "mic.") %>%
  group_by(strain) %>%
  summarize(cr.mic = sum(mic))

# c) MICs relative to drug breakpoints
# Load clinical breakpoint table
bp <- read_csv(file = "./data/breakpoints.csv")
cumu.rel.mic <- cors.master %>% 
  select(c(strain, starts_with("mic"))) %>%
  pivot_longer(!strain, names_to = "drug", values_to = "mic", names_prefix = "mic.") %>%
  left_join(., bp) %>%
  mutate(fold.res = log2(mic) / log2(micbp)) %>%
  group_by(strain) %>%
  summarize(cr.mic.log2rel = sum(fold.res, na.rm = T))
# This is the sum of log2(mic)/log2(bp), further correcting for the 2x dilution series approach to
# MIC testing.

# d) Compare distribution of the scores
score.comp <- mget(ls(pattern = "cumu.*")) %>%
  reduce(right_join, by = "strain") %>%
  pivot_longer(!strain, values_to = "score", names_to = "method")
scorecomp.pl <- ggplot(score.comp, aes(x = score, color = method, fill = method)) +
  geom_density(alpha = 0.25) +
  scale_fill_manual(values = c("#4E05A6", "#018496", "#F57600"), 
                    aesthetics = c("color", "fill"),
                    name = "Score",
                    labels = c("DDsum",
                               "MICsum",
                               "logMICsum"))
scorecomp.pl
# Save plot; uncomment to overwrite
# print_svgpng(scorecomp.pl, "figures/scorecomp", width = 6)

# Save resistance scores to file
# write.csv2(score.comp, file = "tables/resistance_scores.csv")

# Attach cumulative resistances to main data frame
cors.master <- right_join(cumu.res.dd, cors.master)
cors.master <- right_join(cumu.rel.mic, cors.master)



##### 3) Sorted resistance score by origin -----
score.by.origin <- ggplot(cors.master, aes(x = reorder(strain, cr.ddd), cr.ddd, fill = origin_gen)) +
  geom_col() +
  scale_fill_manual(values = c("#b86fc5", "#5eb96c")) +
  labs(x = "strain", y = "DDsum", fill = "Origin") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
score.by.origin

# Save plot; uncomment to overwrite
# print_svgpng(score.by.origin, "figures/score_origin", width = 6)



##### 4) Resistance vs growth -----
# Subset data frame for easier plotting
cors.growth <- cors.master %>%
  select(strain, growthrate, lagphase, maximumod, auc, 
         cr.ddd, cr.mic.log2rel) %>%
  pivot_longer(cols = c(growthrate, lagphase, maximumod, auc), 
               names_to = "growth.param", 
               values_to = "growth.value") %>%
  pivot_longer(cols = c(cr.ddd, cr.mic.log2rel), 
               names_to = "res.param", 
               values_to = "res.score") %>%
  mutate(growth.param = fct_relevel(
    growth.param, c("growthrate", "lagphase", "maximumod", "auc")))

# Plot DDscore against growth parameters
cor.growth.DDscore <- ggplot(subset(cors.growth, res.param == "cr.ddd"), 
                             aes(res.score, growth.value)) +
  facet_wrap(. ~ growth.param, scales = "free") +
  geom_point() +
  stat_cor() +
  labs(x = "DDscore", y = NULL)
cor.growth.DDscore

# Save plot; uncomment to overwrite
# print_svgpng(cor.growth.DDscore, "figures/cor_growth_DDscore", width = 5, height = 4)

# Plot log-transformed relative MIC score against growth parameters
cor.growth.logMICsum <- ggplot(subset(cors.growth, res.param == "cr.mic.log2rel"), 
                               aes(res.score, growth.value)) +
  facet_wrap(. ~ growth.param, scales = "free") +
  geom_point() +
  stat_cor(aes(label = ..r.label..)) +
  labs(x = "logMICsum", y = NULL)
cor.growth.logMICsum
# No stats due to bias introduced by H15 outlier

# Save plot; uncomment to overwrite
# print_svgpng(cor.growth.logMICsum, "figures/cor_growth_logMICscore", width = 5, height = 4)



##### 5) Trait correlation -----

# a) Resistance vs. basic strain info (phylogroup, origin, mobile counts)
# Create a table to store the test results in and run the first test.
# Resistance scores against phylogroups A/B
test.results <- tidy(wilcox.test(cr.mic.log2rel ~ phylogroup, 
                                 subset(cors.master, phylogroup %in% c("A", "B")))) %>%
  add_column(comparison = "phylogroup A/B - logMICsum",
             stat = "W") %>%
  select(comparison, method, stat, statistic, p.value)

test.results <- bind_rows(test.results, 
                          tidy(wilcox.test(cr.ddd ~ phylogroup,
                                           subset(cors.master, phylogroup %in% c("A", "B")))) %>% 
                            add_column(comparison = "phylogroup ABB - DDsum",
                                       stat = "W") %>% 
                            select(comparison, method, stat, statistic, p.value))
# Corresponding plots
ggplot(cors.master, aes(phylogroup, cr.ddd)) + 
  geom_boxplot()
ggplot(cors.master, aes(phylogroup, cr.mic.log2rel)) + 
  geom_boxplot()


# b) Resistance vs. strain origin
test.results <- bind_rows(test.results,
                          tidy(wilcox.test(cr.mic.log2rel ~ origin_gen, cors.master)) %>% 
                            add_column(comparison = "origin - logMICsum",
                                       stat = "W") %>% 
                            select(comparison, method, stat, statistic, p.value))

test.results <- bind_rows(test.results,
                          tidy(wilcox.test(cr.ddd ~ origin_gen, cors.master)) %>% 
                            add_column(comparison = "origin - DDsum",
                                       stat = "W") %>% 
                            select(comparison, method, stat, statistic, p.value))
# Corresponding plots
ggplot(cors.master, aes(origin_gen, cr.mic.log2rel)) + 
  geom_boxplot()
ggplot(cors.master, aes(origin_gen, cr.ddd)) + 
  geom_boxplot()


# c) Resistance vs. plasmid presence
# Only strains with no or one plasmid are included. There is only one strain with >1 plasmid, so
# a statistical analysis is preoblematic. Strains with no plasmid information are excluded.
test.results <- bind_rows(test.results,
                          tidy(wilcox.test(cr.mic.log2rel ~ plasmids, 
                                           subset(cors.master, plasmids == c("0", "1")))) %>%
                            add_column(comparison = "plasmids - logMICsum",
                                       stat = "W") %>%
                            select(comparison, method, stat, statistic, p.value))

test.results <- bind_rows(test.results,
                          tidy(wilcox.test(cr.ddd ~ plasmids, 
                                           subset(cors.master, plasmids == c("0", "1")))) %>%
                            add_column(comparison = "plasmids - DDsum",
                                       stat = "W") %>%
                            select(comparison, method, stat, statistic, p.value))


# Corresponding plots
ggplot(cors.master, aes(as.factor(plasmids), cr.mic.log2rel)) + 
  geom_boxplot() +
  geom_jitter(height = 0, width = 0.25, alpha = 0.33, color = "blue")
ggplot(cors.master, aes(as.factor(plasmids), cr.ddd)) + 
  geom_boxplot() +
  geom_jitter(height = 0, width = 0.25, alpha = 0.33, color = "blue")


# d) Resistance vs. ICE count
# With values for ICE count ranging from 0 - 6, a regression attempt is possible.
ggplot(subset(cors.master, !is.na(cors.master$ICEs)), aes(ICEs, cr.mic.log2rel)) + 
  geom_point()
test.results <- bind_rows(test.results, 
                          tidy(cor.test(cors.master$cr.mic.log2rel, cors.master$ICEs)) %>%
                            add_column(comparison = "logMICsum ~ ICEs",
                                       stat = "t") %>%
                            select(comparison, method, stat, statistic, p.value))

ggplot(subset(cors.master, !is.na(cors.master$ICEs)), aes(ICEs, cr.ddd)) + 
  geom_point()
test.results <- bind_rows(test.results, 
                          tidy(cor.test(cors.master$cr.ddd, cors.master$ICEs)) %>%
                            add_column(comparison = "DDsum ~ ICEs",
                                       stat = "t") %>%
                            select(comparison, method, stat, statistic, p.value))


# e) Growth vs. origin
test.results <- bind_rows(test.results,
                          tidy(wilcox.test(
                            growthrate ~ origin_gen, 
                            cors.master)) %>% 
                            add_column(comparison = "origin - growthrate",
                                       stat = "W") %>% 
                            select(comparison, method, stat, statistic, p.value))
ggplot(cors.master, aes(origin_gen, growthrate)) +
  geom_boxplot()

test.results <- bind_rows(test.results,
                          tidy(wilcox.test(
                            lagphase ~ origin_gen, 
                            cors.master)) %>% 
                            add_column(comparison = "origin - lagphase",
                                       stat = "W") %>% 
                            select(comparison, method, stat,statistic, p.value))
ggplot(cors.master, aes(origin_gen, lagphase)) +
  geom_boxplot()

test.results <- bind_rows(test.results,
                          tidy(wilcox.test(
                            auc ~ origin_gen, 
                            cors.master)) %>% 
                            add_column(comparison = "origin - auc",
                                       stat = "W") %>% 
                            select(comparison, method, stat, statistic, p.value))
ggplot(cors.master, aes(origin_gen, auc)) +
  geom_boxplot()

test.results <- bind_rows(test.results,
                          tidy(wilcox.test(
                            maximumod ~ origin_gen, 
                            cors.master)) %>% 
                            add_column(comparison = "origin - maximumod",
                                       stat = "W") %>% 
                            select(comparison, method, stat, statistic, p.value))
ggplot(cors.master, aes(origin_gen, maximumod)) +
  geom_boxplot()

# f) Growth vs. phylogroup A/B
test.results <- bind_rows(test.results,
                          tidy(wilcox.test(
                            growthrate ~ phylogroup, 
                            subset(cors.master, phylogroup %in% c("A", "B")))) %>%
                            add_column(comparison = "phylogroup A/B - growthrate",
                                       stat = "W") %>% 
                            select(comparison, method, stat, statistic, p.value))
ggplot(cors.master, aes(phylogroup, growthrate)) +
  geom_boxplot()

test.results <- bind_rows(test.results,
                          tidy(wilcox.test(
                            lagphase ~ phylogroup, 
                            subset(cors.master, phylogroup %in% c("A", "B")))) %>%
                            add_column(comparison = "phylogroup A/B - lagphase",
                                       stat = "W") %>% 
                            select(comparison, method, stat, statistic, p.value))
ggplot(cors.master, aes(phylogroup, lagphase)) +
  geom_boxplot()

test.results <- bind_rows(test.results,
                          tidy(wilcox.test(
                            auc ~ phylogroup, 
                            subset(cors.master, phylogroup %in% c("A", "B")))) %>%
                            add_column(comparison = "phylogroup A/B - auc",
                                       stat = "W") %>% 
                            select(comparison, method, stat, statistic, p.value))
ggplot(cors.master, aes(phylogroup, auc)) +
  geom_boxplot()

test.results <- bind_rows(test.results,
                          tidy(wilcox.test(
                            maximumod ~ phylogroup, 
                            subset(cors.master, phylogroup %in% c("A", "B")))) %>%
                            add_column(comparison = "phylogroup A/B - maximumod",
                                       stat = "W") %>% 
                            select(comparison, method, stat, 
                                   statistic, p.value))
ggplot(cors.master, aes(phylogroup, maximumod)) +
  geom_boxplot()


# g) Growth vs. plasmids
# Because only one strain has 3 plasmids, only strains with one or no plasmids are included.
test.results <- bind_rows(test.results,
                          tidy(wilcox.test(
                            growthrate ~ plasmids, 
                            subset(cors.master, plasmids %in% c("0", "1")))) %>%
                            add_column(comparison = "plasmids - growthrate",
                                       stat = "W") %>% 
                            select(comparison, method, stat, statistic, p.value))
ggplot(cors.master, aes(as.factor(plasmids), growthrate)) +
  geom_boxplot()

test.results <- bind_rows(test.results,
                          tidy(wilcox.test(
                            lagphase ~ plasmids, 
                            subset(cors.master, plasmids %in% c("0", "1")))) %>%
                            add_column(comparison = "plasmids - lagphase",
                                       stat = "W") %>% 
                            select(comparison, method, stat, statistic, p.value))
ggplot(cors.master, aes(as.factor(plasmids), lagphase)) + 
  geom_boxplot()

test.results <- bind_rows(test.results,
                          tidy(wilcox.test(
                            auc ~ plasmids, 
                            subset(cors.master, plasmids %in% c("0", "1")))) %>%
                            add_column(comparison = "plasmids - auc",
                                       stat = "W") %>% 
                            select(comparison, method, stat, statistic, p.value))
ggplot(cors.master, aes(as.factor(plasmids), auc)) +
  geom_boxplot()

test.results <- bind_rows(test.results,
                          tidy(wilcox.test(
                            maximumod ~ plasmids, 
                            subset(cors.master, plasmids %in% c("0", "1")))) %>%
                            add_column(comparison = "plasmids - maximumod",
                                       stat = "W") %>% 
                            select(comparison, method, stat, statistic, p.value))
ggplot(cors.master, aes(as.factor(plasmids), maximumod)) +
  geom_boxplot()


# h) Growth vs. ICE count
test.results <- bind_rows(test.results, 
                          tidy(cor.test(cors.master$growthrate, cors.master$ICEs)) %>%
                            add_column(comparison = "ICEs ~ growthrate",
                                       stat = "t") %>%
                            select(comparison, method, stat, statistic, p.value))
ggplot(cors.master, aes(as.factor(ICEs), growthrate)) +
  geom_boxplot()

test.results <- bind_rows(test.results, 
                          tidy(cor.test(cors.master$lagphase, cors.master$ICEs)) %>%
                            add_column(comparison = "ICEs ~ lagphase",
                                       stat = "t") %>%
                            select(comparison, method, stat, statistic, p.value))
ggplot(cors.master, aes(as.factor(ICEs), lagphase)) +
  geom_boxplot()

test.results <- bind_rows(test.results, 
                          tidy(cor.test(cors.master$auc, cors.master$ICEs)) %>%
                            add_column(comparison = "ICEs ~ auc",
                                       stat = "t") %>%
                            select(comparison, method, stat, statistic, p.value))
ggplot(cors.master, aes(as.factor(ICEs), auc)) +
  geom_boxplot()

test.results <- bind_rows(test.results, 
                          tidy(cor.test(cors.master$maximumod, cors.master$ICEs)) %>%
                            add_column(comparison = "ICEs ~ maximumod",
                                       stat = "t") %>%
                            select(comparison, method, stat, statistic, p.value))
ggplot(cors.master, aes(as.factor(ICEs), maximumod)) +
  geom_boxplot()



###### 5g) Family-wise error rate correction of statistical tests -----
test.results$p.corrected <- p.adjust(test.results$p.value, method = "fdr")
test.out <- test.results %>%
  mutate(across(c(statistic, p.corrected), ~ round(.x, c(2,3)))) %>%
  rename(analysis = comparison,
         `fdr-corrected p` = p.corrected) %>%
  mutate(statistic = as.character(statistic)) %>%
  unite(stat, statistic, sep = " = ", col = "test statistic") %>%
  select(!p.value)

# Export table
# write_excel_csv2(test.out, file = "tables/statistical_tests.csv")



##### 6) Plot sign. correlation between no. of ICEs and lag phase
cor.ice.lag <- ggplot(cors.master, aes(ICEs, lagphase)) +
  geom_point() +
  stat_cor(aes(label = paste(..r.label..))) +
  labs(x = "number of ICEs", y = "lag phase")
cor.ice.lag

# Save plot; uncomment to overwrite
# print_svgpng(cor.ice.lag, "figures/cor_ice_lag", width = 5, height = 4)





