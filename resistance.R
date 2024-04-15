####
# mPAct -
# Pseudomonas aeruginosa major clone type panel
#
# Antimicrobial resistance phenotypes analysis
####
# R 4.3.1
# tidyverse 2.0
####



library(ggpubr)
library(tidyverse)
library(svglite)
library(DescTools)
library(broom)

# Drug order vector for plotting
drug.order <- c("PIP", "PIT",
                "CTZ", "CTV", "CEP", "CTT", "CID",
                "IMI", "MER",
                "AZT",
                "CIP",
                "AMI", "GEN", "TOB",
                "COL")

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
# Disc diffusion on M9 agar
res.dd <- read_csv(file = "./data/disc_diffusion_m9.csv")
# Etests on M9 agar
res.etest <- read_csv(file = "./data/mic_etest_m9.csv")
# Reserve antibiotic Etests
res.reserve <- read_csv(file = "./data/mic_etest_reserve.csv")
# Automated antimicrobial susceptibility testing (AST)
res.ast <- read_csv(file = "./data/mic_vitek.csv")
# Microcolony data
microcol <- read_csv(file = "./data/microcolonies.csv")
# Resistance stability data
stability <- read_csv(file = "./data/resistance_stability.csv")

# Load clinical breakpoint table
bp <- read_csv(file = "./data/breakpoints.csv")

# Load resistance classification info
resint <- read_csv(file = "./data/resistance_interpretation.csv")
resint <- pivot_longer(resint, 
                       cols = 2:ncol(resint), 
                       names_to = "drug", 
                       values_to = "interpretation")



##### 2) Resistance heat maps -----

### 2.1) Automated susceptibility testing (AST) data

# Calculate fold resistance breakpoint ratios
res.ast.rel <- left_join(res.ast, bp) %>%
  mutate(fold.res = mic / micbp)

# Add interpretation data
res.ast.rel <- left_join(res.ast.rel, resint)

# Change antibiotics order to match EUCAST suggested grouping and order
res.ast.rel$drug <- factor(res.ast.rel$drug, levels = drug.order)

# Plot heat map
res.ast.heatmap <- ggplot(res.ast.rel, 
                           aes(drug, strain, fill = log2(fold.res), label = interpretation)) +
  geom_tile(color = "white", lwd = 1) +
  scale_fill_gradient2(low = "#2475c2", 
                       mid = "white", 
                       high = "#d72c08", 
                       limits = c(-10, 9),
                       breaks = c(-8, -4, 0, 4, 8)) +
  geom_text(size = 2) +
  coord_fixed() +
  scale_y_discrete(limits = rev) +
  labs(x = NULL, y = NULL, fill = "MIC, log2 fold\nresistance\nbreakpoint") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")
res.ast.heatmap
# Save plot; uncomment to overwrite
# print_svgpng(res.ast.heatmap, "figures/ast_heatmap", height = 5, width = NA)


### 2.2) Reserve antibiotics Etests

# Calculate medians of technical replicates
res.reserve.med <- res.reserve %>%
  group_by(strain, drug) %>%
  summarize(mic.med = median(mic, na.rm = T)) %>%
  rename(mic = mic.med)

# Calculate fold resistance breakpoint ratios
res.reserve.rel <- left_join(res.reserve.med, bp) %>%
  mutate(fold.res = mic / micbp)

# Add interpretation data
res.reserve.rel <- left_join(res.reserve.rel, resint)

# Change antibiotics order to match EUCAST suggested grouping and order
res.reserve.rel$drug <- factor(res.reserve.rel$drug, levels = drug.order)

# Plot heat map
res.reserve.heatmap <- ggplot(res.reserve.rel, 
                          aes(drug, strain, fill = log2(fold.res), label = interpretation)) +
  geom_tile(color = "white", lwd = 1) +
  scale_fill_gradient2(low = "#2475c2", 
                       mid = "white", 
                       high = "#d72c08", 
                       limits = c(-10, 9),
                       breaks = c(-8, -4, 0, 4, 8)) +
  geom_text(size = 2) +
  coord_fixed() +
  scale_y_discrete(limits = rev) +
  labs(x = NULL, y = NULL, fill = "MIC, log2 fold\nresistance\nbreakpoint") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")
res.reserve.heatmap
# Save plot; uncomment to overwrite
# print_svgpng(res.reserve.heatmap, "figures/reserve_heatmap", height = 5, width = NA)


## 2.3) Etests on M9 agar

# Calculate medians of technical replicates
res.etest.med <- res.etest %>%
  group_by(strain, drug) %>%
  summarize(mic.med = median(mic, na.rm = T)) %>%
  rename(mic = mic.med)

# Calculate fold resistance breakpoint ratios
res.etest.rel <- left_join(res.etest.med, bp) %>%
  mutate(fold.res = mic / micbp)

# Change antibiotics order to match EUCAST suggested grouping and order
res.etest.rel$drug <- factor(res.etest.rel$drug, levels = drug.order)

# Plot heat map
res.etest.heatmap <- ggplot(res.etest.rel, 
                              aes(drug, strain, fill = log2(fold.res))) +
  geom_tile(color = "white", lwd = 1) +
  scale_fill_gradient2(low = "#2475c2", 
                       mid = "white", 
                       high = "#d72c08",
                       limits = c(-10, 9),
                       breaks = c(-8, -4, 0, 4, 8)) +
  coord_fixed() +
  scale_y_discrete(limits = rev) +
  labs(x = NULL, y = NULL, fill = "MIC, log2 fold\nresistance\nbreakpoint") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")
res.etest.heatmap
# Save plot; uncomment to overwrite
# print_svgpng(res.etest.heatmap, "figures/etest_heatmap", height = 5, width = NA)


### 2.4) Disc diffusion

# Calculate medians of technical replicates
res.dd.med <- res.dd %>%
  group_by(strain, drug) %>%
  summarize(dd.med = median(diameter, na.rm = T)) %>%
  rename(dd = dd.med)

# Calculate fold resistance breakpoint ratios
res.dd.rel <- left_join(res.dd.med, bp) %>%
  mutate(fold.res = dd - diameterbp) %>%
  # Remove colistin due to lack of breakpoint in disc diffusion
  filter(drug != "COL")

# Change antibiotics order to match EUCAST suggested grouping and order
res.dd.rel$drug <- factor(res.dd.rel$drug, levels = drug.order)

# Plot heat map
res.dd.heatmap <- ggplot(res.dd.rel, 
                              aes(drug, strain, fill = fold.res)) +
  geom_tile(color = "white", lwd = 1) +
  scale_fill_gradient2(low = "#C7522B", 
                       mid = "white", 
                       high = "#3C5941",
                       limits = c(-30, 30),
                       breaks = c(seq(-20, 20, by = 10))) +
  coord_fixed() +
  scale_y_discrete(limits = rev) +
  labs(x = NULL, y = NULL, fill = "mm distance\nto breakpoint") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")
res.dd.heatmap
# Save plot; uncomment to overwrite
# print_svgpng(res.dd.heatmap, "figures/dd_heatmap", height = 5, width = NA)



##### 3) Resistance change when including microcolonies -----

# Calculate change in MIC caused by microcolonies
microcol.m <- microcol %>%
  group_by(strain, drug, time, microcol) %>%
  summarize(med.mic = median(mic, na.rm = T))
microcol.delta <- microcol.m %>%
  pivot_wider(id_cols = c(strain, drug, time), names_from = microcol, values_from = med.mic) %>%
  mutate(delta.microcol = y / n) %>%
  mutate(across(delta.microcol, ~ round(.x, 3)))
microcol.delta$drug <- factor(microcol.delta$drug, levels = drug.order)

# Full overview plot
delta.microcol.pl.full <- ggplot(microcol.delta, aes(delta.microcol, time)) +
  facet_grid(rows = vars(strain), cols = vars(drug), switch = "y") +
  geom_col(position = position_dodge()) +
  geom_vline(xintercept = 1) +
  scale_y_discrete(limits = rev) +
  labs(x = "log2 deltaMIC from microcolonies", y = NULL)
delta.microcol.pl.full


# Heatmap of MIC change when including microcolonies, 48 h data, continuous coloring
microcol.heatmap <- ggplot(subset(microcol.delta, time == "48h"), 
                           aes(drug, strain, fill = delta.microcol)) +
  geom_tile(color = "white", lwd = 1) +
  scale_fill_continuous(low = "white", high = "#b67b0d") +
  coord_fixed() +
  scale_y_discrete(limits = rev) +
  labs(x = NULL, y = NULL, fill = "fold MIC, log,\nwith microcolonies") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")
microcol.heatmap
# Save plot; uncomment to overwrite
# print_svgpng(microcol.heatmap, "figures/microcol_heatmap", height = 5, width = NA)



##### 4) Resistance measure correlations -----

# Merge disc diffusion and Etest data
# No IMI Etest, CTT missing in VITEK
# PTZ, CTZ, MER, AMI, GEN, TOB, CIP, COL are complete in all 3
res.cors <- res.dd %>%
  group_by(strain, drug) %>%
  summarize(ddmean = round(mean(diameter, na.rm = T), 2)) %>%
  pivot_wider(names_from = "drug", values_from = "ddmean", names_prefix = "dd.")
res.cors <- res.etest %>%
  filter(microcol == "n" & time == "24h") %>%
  group_by(strain, drug) %>%
  summarize(micmean = round(mean(mic, na.rm = T), 3)) %>%
  select(strain, drug, micmean) %>%
  pivot_wider(names_from = "drug", values_from = "micmean", names_prefix = "mic.") %>%
  right_join(., res.cors, by = "strain")
# Add AST data, categorized
res.cors <- res.ast %>%
  pivot_wider(names_from = "drug", values_from = "mic", names_prefix = "ast.") %>%
  mutate(across(starts_with("ast"), 
                ~ parse_factor(as.character(.), 
                               levels = as.character(sort(unique(res.ast$mic)))))) %>%
  right_join(., res.cors, by = "strain")

# Dataset for disc diffusion - Etest correlation
cor.dd.etest <- res.cors %>%
  select(strain, starts_with("dd")) %>%
  pivot_longer(!strain, names_to = "drug", names_prefix = "dd.", values_to = "dd")
cor.dd.etest <- res.cors %>%
  select(strain, starts_with("mic")) %>%
  pivot_longer(!strain, names_to = "drug", names_prefix = "mic.", values_to = "mic") %>%
  right_join(., cor.dd.etest) %>%
  drop_na()
# Change antibiotics order to match EUCAST suggested grouping and order
cor.dd.etest$drug <- factor(cor.dd.etest$drug, levels = drug.order)

# Regression line p values
linreg.pvals <- cor.dd.etest %>%
  group_by(drug) %>%
  do(glance(lm(mic ~ dd, data = .))) %>%
  select(drug, r.squared, statistic, p.value)
# Correct for multiple testing
linreg.pvals$p.corr <- p.adjust(linreg.pvals$p.value, method = "fdr")
linreg.pvals <- mutate(linreg.pvals, sig = case_when(p.corr < 0.001 ~ "***",
                                                     p.corr >= 0.001 & p.corr < 0.01 ~ "**",
                                                     p.corr >= 0.01 & p.corr < 0.05 ~ "*",
                                                     .default = "ns"))

# Plot
cor.dd.mic.plot <- ggplot(cor.dd.etest, aes(dd, log2(mic))) +
  facet_wrap(. ~ drug, nrow = 1, scales = "free", labeller = labeller(drug = drug.full.names)) +
  geom_point() +
  stat_smooth(method = "lm", formula = y ~ x, geom = "smooth", se = FALSE) +
  stat_regline_equation(label.y = 2, aes(label = after_stat(rr.label)), size = 3) +
  theme(panel.grid.major.y = element_blank()) +
  labs(x = "disc diffusion diameter", y = "log2(MIC)")
cor.dd.mic.plot
# Save plot; uncomment to overwrite. Location of R-squared labels needs manual correction.
# print_svgpng(cor.dd.mic.plot, "figures/cor_dd_mic", width = 14, height = 2)



# Dataset for disc diffusion and Etest against categorical VITEK data
cors.ast <- res.cors %>%
  select(strain, starts_with("ast")) %>%
  pivot_longer(!strain, names_to = "drug", names_prefix = "ast.", values_to = "ast") %>%
  right_join(., cor.dd.etest) %>%
  drop_na()
# Change antibiotics order to match EUCAST suggested grouping and order
cors.ast$drug <- factor(cors.ast$drug, levels = drug.order)

# Plot VITEK against disc diffusion
ast.dd.plot <- ggplot(cors.ast, aes(ast, dd)) +
  geom_jitter(width = 0.25, height = 0, alpha = 0.5, shape = 1) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(. ~ drug, nrow = 1, labeller = labeller(drug = drug.full.names)) +
  labs(x = "AST MIC (mg/l)", y = "disc diffusion\ndiameter (mm)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ast.dd.plot
# Save plot; uncomment to overwrite
# print_svgpng(ast.dd.plot, "figures/cor_dd_ast", width = 12, height = 4)

ast.etest.plot <- ggplot(cors.ast, aes(ast, log2(mic))) +
  geom_jitter(width = 0.25, height = 0, alpha = 0.5, shape = 1) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(. ~ drug, nrow = 1, labeller = labeller(drug = drug.full.names)) +
  labs(x = "AST MIC (mg/l)", y = "Etest MIC\n(log2 mg/l)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ast.etest.plot
# Save plot; uncomment to overwrite
# print_svgpng(ast.etest.plot, "figures/cor_etest_ast", width = 12, height = 4)



# Export resistance values to a table
# write_excel_csv2(cors.ast, file = "tables/resistance_measurements.csv")



##### 5) Resistance stability -----

# Get ancestral MICs from the Etest data
stability.ancestral <- res.etest %>%
  filter(microcol == "n", time == "24h") %>%
  filter(drug %in% c("CIP", "GEN", "MER", "PIT", "TOB")) %>%
  group_by(strain, drug) %>%
  summarize(anc_mic = mean(mic, na.rm = T))

# Combine ancestral and evolved MICs into one data frame
stability.all <- stability %>%
  group_by(strain, drug) %>%
  summarize(mean_mic = mean(mic, na.rm = T)) %>%
  right_join(stability.ancestral, .) %>%
  select(strain, drug, mean_mic, anc_mic) %>%
  rename(mic_evol = mean_mic, mic_anc = anc_mic)

# Make a data frame for before and after resistance classification
stability.sir.class.switch <- stability.all %>%
  pivot_longer(cols = starts_with("mic"), names_to = "state", values_to = "mic") %>%
  right_join(select(bp, drug, micbp), ., by = "drug") %>%
  # Special case: MER has an 'I' category at MIC >2 to =<8, check if this applies
  mutate(sir = case_when(drug == "MER" & mic <= 2 ~ "S",
                         drug == "MER" & mic > 2 & mic <= 8 ~ "I",
                         drug == "MER" & mic > 8 ~ "R",
                         # For other drugs, only an 'R' breakpoint exists. For simplicity, 'R' and
                         # 'S' are used for all cases, even when only 'I' exists according to
                         # EUCAST breakpoints
                         drug != "MER" & mic > micbp ~ "R",
                         drug != "MER" & mic <= micbp ~ "S")) %>%
  pivot_wider(id_cols = , names_from = "state", values_from = c(mic, sir)) %>%
  filter(sir_mic_evol != sir_mic_anc) %>%
  select(strain, drug, sir_mic_anc, sir_mic_evol, mic_mic_anc, mic_mic_evol) %>%
  rename(mic_anc = mic_mic_anc, mic_evol = mic_mic_evol, 
         sir_anc = sir_mic_anc, sir_evol = sir_mic_evol)
# Write resistance category changes to a table
# write_excel_csv(stability.sir.class.switch, "tables/resistance_stability_switches.csv")

# Bind ancestral MICs to data frame
stability.calc <- stability.all %>%
  mutate(delta_mic = log2(mic_evol) - log2(mic_anc)) %>%
  mutate(across(drug, ~ factor(.x, levels = drug.order))) %>%
  # Intruduce small values to deltaMIC = 0 instances for log transformation
  mutate(delta_mic_plot = delta_mic + 0.00001)

# Calculate the sum of resistance changes
stability.sums <- stability.calc %>%
  group_by(strain) %>%
  summarize(stab.sum = round(sum(delta_mic), 2))

# Save changes in MIC to a new file
stability.save <- stability.calc %>%
  select(strain, drug, delta_mic) %>%
  mutate(across(delta_mic, ~ round(.x, 2)))
# write_excel_csv(stability.save, "tables/resistance_stability_micchange.csv")

# Plot (barplot)
stability.pl <- ggplot(stability.calc, aes(drug, delta_mic, fill = drug)) +
  facet_wrap(. ~ strain, nrow = 1) +
  geom_col() +
  geom_text(data = stability.sums, 
            mapping = aes(x = 3, y = 1.5, label = stab.sum), 
            inherit.aes = FALSE,
            size = 3) +
  scale_fill_viridis_d(option = "cividis") +
  geom_hline(yintercept = 0) +
  labs(y = "log2(MICevol - MICanc) (mg/l)") +
  theme_minimal() +
  theme(axis.text.x = element_blank(), 
        legend.position = "bottom", 
        panel.grid = element_blank(), 
        axis.title.x = element_blank())
stability.pl

# Save plot; uncomment to overwrite
# print_svgpng(stability.pl, "figures/res_stability", width = 10, height = 4)



##### 6) Time-kill curves -----

# Read and check the data table
tk.raw <- read_csv(file = "./data/timekill_raw.csv")

# Insert minimal value for logarithmic plotting and calculations
tk.raw <- tk.raw %>%
  mutate(cfu = cfu + 1)

# Plots of CTZ in H10 and H13 to illustrate time-kill AUC approach
# MIC of both strains is 2 mg/l, concentration used is 200 mg/l
tk.auc.illus.raw.df <-  filter(tk.raw, 
                               drug == "CTZ", strain %in% c("H10", "H13"), replicate == "rep_2")
tk.auc.illus.raw.pl <- ggplot(tk.auc.illus.raw.df, aes(time, log10(cfu))) +
  facet_wrap(. ~ strain) +
  geom_line() +
  scale_x_continuous(breaks = c(0, 2.5, 5), labels = c("0", "2.5", "5")) +
  labs(x = "time (h)", y = "log10 cfu/ml") +
  ggtitle("Time-kill curve at\n100x MIC ceftazidime")
tk.auc.illus.raw.pl
# Save plot; uncomment to overwrite
# print_svgpng(tk.auc.illus.raw.pl, "figures/tkauc_illus", width = 4, height = 4)

# Calculate the area under the curves per replicate
tk.raw.auc <- tk.raw %>%
  group_by(strain, drug, replicate) %>%
  summarize(auc = AUC(time, log10(cfu), method = "trapezoid"))

# Calculate mean and sd for the replicates
tk.raw.auc.stats <- tk.raw.auc %>%
  ungroup() %>%
  group_by(strain, drug) %>%
  summarize(auc.m = mean(auc, na.rm = T), auc.sd = sd(auc, na.rm = T))

# Save time-kill AUCs as csv
# write_excel_csv(tk.raw.auc.stats, "tables/timekill_auc_stats.csv")

# Split off the reference strain data for the horizontal lines
tk.auc.refs.means <- tk.raw.auc.stats %>%
  filter(strain %in% c("PA14", "PAO1")) %>%
  select(strain, drug, auc.m)

# Filter the remaining strain data for plotting
tk.auc.stats.pl <- tk.raw.auc.stats %>%
  filter( ! strain %in% c("PA14", "PAO1"))
tk.auc.pl <- tk.raw.auc %>%
  filter( ! strain %in% c("PA14", "PAO1"))

# Plot
tk.auc.stats.bar <- ggplot(tk.auc.stats.pl, aes(strain, auc.m)) +
  geom_hline(data = tk.auc.refs.means, 
             aes(yintercept = auc.m, linetype = strain), linewidth = 0.3) +
  geom_col(width = 0.6) +
  geom_point(data = tk.auc.pl, aes(y = auc), 
             shape = 5, size = 1, position = position_dodge2(width = 0.3)) +
  geom_linerange(aes(ymax = auc.m + auc.sd, ymin = auc.m), color = "gray45") +
  facet_wrap(. ~ drug, labeller = labeller(drug = drug.full.names)) +
  scale_linetype_manual(values = c("dotted", "dashed"), name = NULL) +
  labs(y = "AUC") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
tk.auc.stats.bar

# Save plot; uncomment to overwrite
# print_svgpng(tk.auc.stats.bar, "figures/timekill_bar", width = 10)

# Correlate AUC and MIC
tk.raw.mic <- res.dd %>%
  group_by(strain, drug) %>%
  summarize(mean.dd = mean(diameter)) %>%
  right_join(., tk.raw.auc.stats)
tk.auc.raw.dd.cor <- ggplot(tk.raw.mic, aes(mean.dd, auc.m)) +
  facet_wrap(. ~ drug, labeller = labeller(drug = drug.full.names), scales = "free_x") +
  geom_point(size = 1) +
  stat_cor() +
  labs(x = "disc diffusion (mm)", y = "AUC of time-kill curve")
tk.auc.raw.dd.cor

# Save plot; uncomment to overwrite
# print_svgpng(tk.auc.raw.dd.cor, "figures/tk_dd_cors", width = 5, height = 4)




