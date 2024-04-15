####
# mPAct -
# Pseudomonas aeruginosa major clone type panel
#
# Analysis of growth curves (M9 minimal medium, plate reader kinetic)
####
# R 4.3.1
# tidyverse 2.0
####



library(tidyverse)
library(growthrates)
library(DescTools)
library(svglite)


# Plot exporting wrapper function
# Same size png and svg of the specified plot. Defaults to golden ratio of height to the input
# width. Path as string with path, but file ending omitted.
print_svgpng <- function(plot, path = "file", width, height = width * 0.618) {
  ggsave(plot, filename = paste0(path, ".svg"), width = width, height = height, units = "in",
         device = svglite, fix_text_size = FALSE)
  ggsave(plot, filename = paste0(path, ".png"), width = width, height = height, units = "in")
}

# Load data
gc <- read_csv(file = './data/growth_curves.csv')



### Growth curve data table for export ----

gc.exp <- gc %>%
  group_by(strain, time) %>%
  summarize(mean_od = mean(od, na.rm = T), sd_od = sd(od, na.rm = T)) %>%
  mutate(across(ends_with("od"), ~ format(signif(., 3)))) %>%
  rename(od_mean = mean_od, sd = sd_od)
# Export OD per time point table, uncomment to overwrite. 
# write_excel_csv2(gc.exp, file = "./tables/growth_od.csv")



### Analysis ----

# ____Calculate growth parameters per strain and replicate ----

# Model maximum growth rate and lag phase, linear method
gc.linear <- gc
fit.linear <- all_easylinear(od ~ time | strain + replicate, data = gc.linear, h = 10, quota = 0.95)
growth.linear <- results(fit.linear)
hist(growth.linear$r2)
# A cutoff of 0.9 is used for the model fit. Find and remove the replicates below this from the
# base data set.
growth.linear.poor.fit <- filter(growth.linear, r2 < 0.9) %>%
  rownames_to_column(var = "code")
artifacts <- growth.linear.poor.fit$code
print(artifacts) # H15:6, H05:7, H01:16, H01:29, H01:32
gc <- gc %>%
  mutate(code = paste(strain, replicate, sep = ":")) %>%
  filter(! code %in% artifacts)

# Rerun the fit
fit.linear <- all_easylinear(od ~ time | strain + replicate, data = gc, h = 10, quota = 0.95)
growth.linear <- results(fit.linear)
hist(growth.linear$r2)
# All replicates below 0.9 removed.

# Extract maximum OD for each replicate from the raw data
max.od <- gc %>%
  group_by(strain, replicate) %>%
  summarize(max.od = max(od))
ggplot(max.od, aes(strain, max.od)) + 
  geom_boxplot()

# Calculate area under the growth curve for each replicate from the raw data
auc <- gc %>%
  group_by(strain, replicate) %>%
  summarize(auc = AUC(time, od, method = "spline"))
ggplot(auc, aes(strain, auc)) + 
  geom_boxplot()

# Combine modeled and calculated data into one data frame
growth.params <- inner_join(growth.linear, auc) %>%
  select(strain, replicate, mumax, lag, auc) %>%
  inner_join(., max.od) %>%
  rename("growthrate" = mumax,
         "lagphase" = lag,
         "maximumod" = max.od)


# Clean up
rm(fit.linear, gc.linear, growth.linear.poor.fit, artifacts, max.od, auc)



### Collect parameters and export as table----
# Calculate means and standard errors for all parameters per strain
growth.params.stats <- growth.params %>%
  group_by(strain) %>%
  summarize(across(growthrate:maximumod, list(mean = mean, sd = sd))) %>%
  mutate(across(starts_with("growthrate"), ~ round(.x, 4))) %>%
  mutate(across(4:ncol(.), ~ round(.x, 2)))

# Export table as file; uncomment to overwrite
# write_excel_csv2(growth.params.stats, "tables/growth_parameters.csv")



### Growth curve figures----
# Calculate mean and standard error per strain and time point
gc.plot.df <- gc %>%
  group_by(strain, time) %>%
  summarize(mean.od = mean(od),
            sem.od = sd(od, na.rm = T) / sqrt(n())) %>%
  mutate(time.h = time / 60)

# Plot with non-logarithmic y axis
growth.plot <- ggplot(gc.plot.df, aes(time.h, mean.od)) +
  geom_ribbon(aes(ymax = mean.od + sem.od, ymin = mean.od - sem.od), fill = "grey50") +
  geom_line(linewidth = 1) +
  facet_wrap(. ~ strain) +
  labs(title = "Growth curves (M9 minimal medium)", x = "Time (h)", y = "OD")
growth.plot

# Plot with logarithmic y axis
growth.plot.log <- ggplot(gc.plot.df, aes(time.h, log2(mean.od))) +
  geom_ribbon(aes(ymax = log2(mean.od + sem.od), ymin = log2(mean.od - sem.od)), fill = "grey50") +
  geom_line(linewidth = 1) +
  facet_wrap(. ~ strain) +
  labs(title = "Growth curves (M9 minimal medium)", x = "Time (h)", y = "log2 OD")
growth.plot.log



### Main growth data bar plot----
# Split off the reference strain data for the horizontal lines
growth.refs.means <- growth.params.stats %>%
  filter(strain %in% c("PA14", "PAO1")) %>%
  select(strain, ends_with("mean")) %>%
  rename_with(., ~ sub("_mean", "", .x)) %>%
  pivot_longer(cols = 2:5, names_to = "param", values_to = "value") %>%
  # Set the order of the parameters for the plot facets
  mutate(param = factor(param)) %>%
  mutate(param = fct_relevel(param, c("growthrate", "lagphase", "maximumod", "auc")))
  
# Format the remaining strain data into a long format table for plotting
growth.hstrains <- growth.params %>%
  filter( ! strain %in% c("PA14", "PAO1"))
growth.hstrains <- growth.hstrains %>%
    pivot_longer(cols = 3:6, names_to = "param", values_to = "value")
# Set the order of the parameters for the plot facets
growth.hstrains <- growth.hstrains %>%
  mutate(param = factor(param)) %>%
  mutate(param = fct_relevel(param, c("growthrate", "lagphase", "maximumod", "auc")))

# Calculate stats bars and error bars
growth.hstrains.stats <- growth.hstrains %>%
  group_by(strain, param) %>%
  summarize(mean = mean(value, na.rm = T),
            sem = sd(value, na.rm = T) / sqrt(n()))

# Plot
growth.params.barplot <- ggplot(growth.hstrains.stats, aes(strain, mean)) +
  geom_hline(data = growth.refs.means, 
             aes(yintercept = value, linetype = strain), linewidth = 0.3) +
  geom_col(width = 0.6) +
  geom_linerange(aes(ymax = mean + sem, ymin = mean), color = "gray45") +
  geom_point(data = filter(growth.hstrains, !strain %in% c("PA14", "PAO1")), 
             aes(y = value), shape = 5, size = 0.6, position = position_dodge2(width = 0.3)) +
  facet_wrap(vars(param), scales = "free_y") +
  scale_linetype_manual(values = c("dotted", "dashed"), name = NULL) +
  theme(axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))
growth.params.barplot

# Save plot; uncomment to overwrite
# print_svgpng(growth.params.barplot, "figures/growth_parameters_bar", width = 10)

# Clean up
rm(gc.plot.df, growth.hstrains, growth.linear, growth.params.barplot, growth.params.boxplot,
   growth.plot, growth.plot.log)



### cfu count data----

# Load data
cfu <- read_csv2(file = "./data/cfu_stat_phase.csv")

# Function to transform dilution level and count data into cfu values
count_to_cfu <- function(dilution, count, spotvol) {
  cfu <- round((count / spotvol) * 10^(dilution + 3), 0)
  return(cfu)
}

# Calculate cfu counts per ml
cfu.counts <- cfu %>%
  mutate(cfu = mapply(count_to_cfu, dil, count, spotvol = 7)) 

# Calculate means of technical replicates
cfu.bioreps <- cfu.counts %>%
  group_by(strain, biorep) %>%
  summarize(mean.cfu = mean(cfu, na.rm = T))

# Calculate mean cfu/ml per strain and run Kruskal-Wallis-Test
mean.cfu <- cfu.bioreps %>%
  group_by(strain) %>%
  summarize(mcfu = mean(mean.cfu, na.rm = T))
kruskal.test(mcfu ~ strain, data = mean.cfu)
# p = 0.46

# Write a summary table and export
cfu.exp <- cfu.bioreps %>%
  group_by(strain) %>%
  summarize(cfu_mean = round(mean(mean.cfu, na.rm = T), 0), cfu_sd = sd(mean.cfu, na.rm = T))
# write_excel_csv2(cfu.exp, "tables/cfu_stationary.csv")
