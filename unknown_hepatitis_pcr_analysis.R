#### Unknown hepatitis PCR Data analysis ####

library(tidyverse)
library(ggsignif)
library(cowplot)
library(ggpubr)
library(janitor)

# Define thresholds for positive, low-level positive and negative results
not_detected <- 45 # Max Ct value of the PCR
thres <- 38 # Ct value above which we count as low-level positive
set_zero_to <- 50 # Ct value to set samples with no PCR results for plotting purposes

# Function to calculate positive and negative from threshold
pcr_pos <- function(ct, thres, not_detected) {
  
  if (is.na(ct)) {
    return(NA)
    
  }else if (ct > not_detected ){
    return("negative")
    
  } else if (ct >= thres & ct <= not_detected) {
    return("LLP")
    
  } else if (ct < thres) {
    return("positive")
  }
}

pcr_pos <- Vectorize(pcr_pos, vectorize.args = "ct")

#### Import and prepare data ####
data <- read.csv("pcr_data.csv") %>%
  
  # Set negative PCR results to an arbitrary value above the 
  # max number of CT cycles for plotting
  dplyr::mutate(aav2_ct = ifelse(aav2_ct == "ND", 50, aav2_ct),
         adeno_ct = ifelse(adeno_ct == "ND", 50, adeno_ct),
         hhv6_ct = ifelse(hhv6_ct == "ND", 50, hhv6_ct),
         
  # Calculate inverse CT values
         aav2_ct_inv = 1/as.numeric(aav2_ct),
         adeno_ct_inv = 1/as.numeric(adeno_ct), 
         hhv6_ct_inv = 1/as.numeric(hhv6_ct),
  
  # Assign positive/negative/low-level positive to CT values
         aav2_pcr = pcr_pos(ct = aav2_ct, thres = thres, not_detected = not_detected),
         adeno_pcr = pcr_pos(ct = adeno_ct, thres = thres, not_detected = not_detected),
         hhv6_pcr = pcr_pos(ct = hhv6_ct, thres = thres, not_detected = not_detected))


#### Statistical analysis ####

# Liver
all_liver <- data %>%
  dplyr::filter(tissue == "liver") %>%
  dplyr::select(c(group, adeno_ct_inv, aav2_ct_inv, hhv6_ct_inv)) %>% 
  tidyr::pivot_longer(!group, names_to = "virus", values_to = "ct_inv") %>% 
  dplyr::mutate(virus = sub("_ct_inv", "", virus)) %>% 
  dplyr::mutate(virus = factor(virus, levels = c("adeno", "aav2", "hhv6"), 
                        labels = c("AdV", "AAV2", "HHV6")),
         group = factor(group, levels = c("case", "control"), 
                        labels = c("Case", "Control*")))
         
# Wilcox tests, case vs control, for each virus
wilcox.test(ct_inv ~ group, dplyr::filter(all_liver, virus == "AdV"))
wilcox.test(ct_inv ~ group, dplyr::filter(all_liver, virus == "AAV2"))
wilcox.test(ct_inv ~ group, dplyr::filter(all_liver, virus == "HHV6"))

# Blood
control_types_aav2 <- c("case", "Case", "control_liver", "control_stool", "Adenovirus viraemia", "CMV viraemia",
                        "EBV viraemia", "HHV6 viraemia", "Adenovirus and CMV viraemia",
                        "CMV and EBV viraemia", "PIMS-TS syndrome",
                        "SARS-CoV-2 positive", "non_infectious",
                        "Adenovirus with normal ALT",
                        "Adenovirus with normal ALT (blood)", 
                        "Adenovirus with raised ALT", "Critical Illness with raised ALT",
                        "Non-adenovirus raised ALT", "Adenovirus CovidMisC normal ALT", 
                        "Other hepatitis")

control_types_aav2_labels <- c("Case", "Case", "Control*", "Control", "AdV, raised\nALT (GOSH)*", "CMV, raised ALT*",
                               "EBV, raised ALT*", "HHV6, raised ALT*", "AdV and CMV,\nraised ALT*",
                               "CMV and EBV,\nraised ALT*", "MIS-C",
                               "SARS-CoV-2*", "Healthy", 
                               "AdV (elsewhere), normal ALT",
                               "AdV (blood), normal ALT", 
                               "AdV, raised ALT (D/P)", "Critical Illness, raised ALT",
                               "Non-AdV, raised ALT", "AdV CovidMisC, normal ALT", 
                               "Other hepatitis")

control_types <- c("Case", "Control", "control_liver", "control_stool", "Adenovirus viraemia", "CMV viraemia",
                   "EBV viraemia", "HHV6 viraemia", "Adenovirus and CMV viraemia",
                   "CMV and EBV viraemia", "PIMS-TS syndrome", "SARS-CoV-2 positive")

control_types_labels <- c("Case", "Control", "Control*", "Control", "AdV, raised\nALT (GOSH)*", "CMV, raised ALT*",
                          "EBV, raised ALT*", "HHV6, raised ALT*", "AdV and CMV,\nraised ALT*",
                          "CMV and EBV,\nraised ALT*", "MIS-C",
                          "SARS-CoV-2*")

# AAV2 in blood
aav2_blood <- data %>% 
  dplyr::filter(tissue == "blood" & !is.na(aav2_ct)) %>%  
  dplyr::select(c(group, control_type, aav2_ct, aav2_ct_inv, source, age)) %>% 
  dplyr::mutate(control_type = factor(control_type, levels = control_types_aav2,
                               labels = control_types_aav2_labels))

aav2_blood_gosh <- aav2_blood %>% 
  dplyr::filter(!(source %in% c("PERFORM", "DIAMONDS"))) 

aav2_blood_pd <- aav2_blood %>% 
  dplyr::filter(group == "case" | source %in% c("PERFORM", "DIAMONDS"))

# Cases vs all controls
wilcox.test(aav2_ct_inv ~ group, aav2_blood)
wilcox.test(aav2_ct_inv ~ group, aav2_blood_gosh)
wilcox.test(aav2_ct_inv ~ group, aav2_blood_pd)

# Cases vs individual groups of controls
kruskal.test(aav2_ct_inv ~ control_type, aav2_blood)
pairwise.wilcox.test(aav2_blood$aav2_ct_inv, aav2_blood$control_type, 
                     p.adjust.method = "BH")

# Treating all AdV groups in DIAMONDS/PERFORM as single group
pd_adv_control_groups <- c("AdV (elsewhere), normal ALT", 
                           "AdV (blood), normal ALT", 
                           "AdV, raised ALT (D/P)", 
                           "AdV CovidMisC, normal ALT")

aav2_adv_group <- aav2_blood %>% 
  dplyr::mutate(control_type = as.character(control_type)) %>% 
  dplyr::mutate(control_type = ifelse(control_type %in% pd_adv_control_groups, 
                               "adv_group", control_type))

kruskal.test(aav2_ct_inv ~ control_type, aav2_adv_group)
pairwise.wilcox.test(aav2_adv_group$aav2_ct_inv, 
                     aav2_adv_group$control_type, 
                     p.adjust.method = "BH")

# Under 10s only
aav2_blood_u10 <- aav2_blood %>% 
  dplyr::filter(age < 10 | is.na(age))

kruskal.test(aav2_ct_inv ~ control_type, aav2_blood_u10)
pairwise.wilcox.test(aav2_blood_u10$aav2_ct_inv, aav2_blood_u10$control_type, 
                     p.adjust.method = "BH")

aav2_blood_gosh_u10 <- aav2_blood_gosh %>% 
  dplyr::filter(age < 10 | is.na(age))

aav2_blood_pd_u10 <- aav2_blood_pd %>% 
  dplyr::filter(age < 10 | is.na(age))


# AdV in blood

adeno_blood <- data %>% 
  dplyr::filter(tissue == "blood" & !is.na(adeno_ct)) %>%  
  dplyr::select(c(group, control_type, adeno_ct, adeno_ct_inv)) %>% 
  dplyr::mutate(control_type = factor(control_type, levels = control_types, 
                               labels = control_types_labels))

# Cases vs all controls
wilcox.test(adeno_ct_inv ~ group, adeno_blood)

# Cases vs individual control groups
kruskal.test(adeno_ct_inv ~ control_type, adeno_blood)
pairwise.wilcox.test(adeno_blood$adeno_ct_inv, adeno_blood$control_type, 
                     p.adjust.method = "BH")

# HHV6 blood 

hhv6_blood <- data %>% 
  dplyr::filter(tissue == "blood" & !is.na(hhv6_ct)) %>%  
  dplyr::select(c(group, control_type, hhv6_ct, hhv6_ct_inv)) %>% 
  dplyr::mutate(control_type = factor(control_type, levels = control_types, 
                               labels = control_types_labels))

# Case vs all controls
wilcox.test(hhv6_ct_inv ~ group, hhv6_blood)

kruskal.test(hhv6_ct_inv ~ control_type, hhv6_blood)
pairwise.wilcox.test(hhv6_blood$hhv6_ct_inv, hhv6_blood$control_type, 
                     p.adjust.method = "BH")

#### CT Plots ####

# Liver
all_liver_sig <- ggplot2::ggplot(data=all_liver, aes(x=factor(group), y=ct_inv)) +
  geom_point(position=position_jitter(w=0.2, h=0),
             shape = 16, alpha = 0.8, size = 0, color = "white") +
  facet_grid(. ~ virus) +
  theme_light() +
  theme(axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position="none",
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank(),
        strip.text.x = element_blank(),
        plot.margin = margin(2,0,11.5,2, "cm")) +
  ggsignif::geom_signif( comparisons = list(c("Case", "Control*")),
               map_signif_level = TRUE, na.rm = TRUE, tip_length = 0.1) + 
  scale_y_continuous(limits = c(1/51, 1/15))


dat_text <- data.frame(
  virus = c("AdV", "AAV2", "HHV6"),
  label = c("LLP", "", "")) %>% 
  dplyr::mutate(virus = factor(virus, levels = c("AdV", "AAV2", "HHV6"), 
                        labels = c("AdV", "AAV2", "HHV6")))


all_liver_plot <- ggplot2::ggplot(data=all_liver, aes(x=factor(group), y=ct_inv)) +
  geom_boxplot( width = 0.3, outlier.shape = NA) +
  geom_point(position=position_jitter(w=0.2, h=0), shape = 16, alpha = 0.8, size = 2) +
  facet_grid(. ~ virus) +
  theme_light() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 16),
        axis.title.y = element_text(angle = 0, size = 16),
        axis.text.y = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(),
        legend.position = "none",
        strip.text.x = element_text(size = 16),
        plot.caption = element_text(hjust = 0.5)) +
  scale_x_discrete(labels=c("Case", "Control*")) +
  ylab("1/Ct") +
  geom_hline(yintercept = 1/thres, linetype = "dashed") +
  geom_hline(yintercept = 1/not_detected, linetype = "dashed") +
  geom_text(data = dat_text, mapping = aes(x = Inf, y = Inf, label = label, hjust = 8, vjust = 38)) +
  scale_y_continuous(breaks = c(1/50, 1/40, 1/35, 1/30, 1/25, 1/20, 1/15), 
                     labels = c("ND", "1/40", "1/35", "1/30", "1/25", "1/20", "1/15"),
                     limits=c(1/51, 1/14))

# png("Figures/test.png", units="in", width=8, height=6, res=300)
# ggdraw(all_liver_plot) + draw_plot(all_liver_sig)
# dev.off()

# AAV2 in blood - cases and GOSH controls
aav2_blood_gosh_plot <- ggplot2::ggplot(data=aav2_blood_gosh, aes(x=factor(control_type), y=aav2_ct_inv)) +
  geom_boxplot(width = 0.3, outlier.shape = NA, show.legend = FALSE) +
  geom_point(position=position_jitter(w=0.2, h=0), shape = 16, alpha = 0.8, size = 2) +
  theme_light() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 16, angle = 45, hjust = 1),
        axis.title.y = element_text(angle = 0, size=18),
        axis.text.y = element_text(size = 16),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border=element_blank(),
        axis.line = element_line(),
        axis.ticks.y = element_line(color = "black")) +
  ylab("1/Ct") +
  geom_hline(yintercept = 1/thres, linetype = "dashed") +
  annotate("text", x = 0.7, y = 1/thres + 0.002, label = "LLP", size = 5) +
  geom_hline(yintercept = 1/not_detected, linetype = "dashed") +
  ggpubr::stat_compare_means(method = "wilcox.test", label = "p.signif", ref.group = "Case",
                     label.y = 1/17, size = 6,
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.0025, 0.02, 1), 
                                        symbols = c("***", "**", "*", "NS"))) +
  scale_y_continuous(breaks = c(1/50, 1/40, 1/35, 1/30, 1/25, 1/20, 1/15), 
                     labels = c("ND", "1/40", "1/35", "1/30", "1/25", "1/20", "1/15"),
                     limits=c(1/51, 1/14))

# png("Figures/aav2_blood_gosh.png", units="in", width=9, height=5, res=300)
# aav2_blood_gosh_plot
# dev.off()

# AAV2 in blood and PERFORM/DIAMONDS controls
aav2_blood_plot_pd <- ggplot2::ggplot(data=aav2_blood_pd, aes(x=factor(control_type), y=aav2_ct_inv)) +
  geom_boxplot(width = 0.3, outlier.shape = NA, show.legend = FALSE) +
  geom_point(position=position_jitter(w=0.2, h=0),
             shape = 16, alpha = 0.8, size = 2) +
  theme_light() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 16, angle = 45, hjust = 1),
        axis.title.y = element_text(angle = 0, size=18),
        axis.text.y = element_text(size = 16),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border=element_blank(),
        axis.line = element_line(),
        axis.ticks.y = element_line(color = "black")) +
  ylab("1/Ct") +
  geom_hline(yintercept = 1/thres, linetype = "dashed") +
  annotate("text", x = 0.7, y = 1/thres + 0.002, label = "LLP", size = 5) +
  geom_hline(yintercept = 1/not_detected, linetype = "dashed") +
  ggpubr::stat_compare_means(method = "wilcox.test", label = "p.signif", ref.group = "Case",
                     label.y = 1/17, size = 6, 
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 1), 
                                        symbols = c("***", "**", "*", "NS"))) +
  scale_y_continuous(breaks = c(1/50, 1/40, 1/35, 1/30, 1/25, 1/20, 1/15), 
                     labels = c("ND", "1/40", "1/35", "1/30", "1/25", "1/20", "1/15"),
                     limits=c(1/51, 1/14))

# png("Figures/aav2_blood_pd.png", units="in", width=9, height=6, res=300)
aav2_blood_plot_pd
# dev.off()

# AAV2 in blood, cases vs GOSH controls, under 10 only
aav2_blood_gosh_plot_u10 <- ggplot2::ggplot(data=aav2_blood_gosh_u10, aes(x=factor(control_type), y=aav2_ct_inv)) +
  geom_boxplot(width = 0.3, outlier.shape = NA, show.legend = FALSE) +
  geom_point(position=position_jitter(w=0.2, h=0), shape = 16, alpha = 0.8, size = 2) +
  theme_light() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 16, angle = 45, hjust = 1),
        axis.title.y = element_text(angle = 0, size=18),
        axis.text.y = element_text(size = 16),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border=element_blank(),
        axis.ticks.y = element_line(color = "black")) +
  ylab("1/Ct") +
  geom_hline(yintercept = 1/thres, linetype = "dashed") +
  annotate("text", x = 0.7, y = 1/thres + 0.002, label = "LLP", size = 5) +
  geom_hline(yintercept = 1/not_detected, linetype = "dashed") +
  ggpubr::stat_compare_means(method = "wilcox.test", label = "p.signif", ref.group = "Case",
                     label.y = 1/17, size = 6,
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 1), 
                                        symbols = c("***", "**", "*", "NS"))) +
  scale_y_continuous(breaks = c(1/50, 1/40, 1/35, 1/30, 1/25, 1/20, 1/15), 
                     labels = c("ND", "1/40", "1/35", "1/30", "1/25", "1/20", "1/15"),
                     limits=c(1/51, 1/14))

# png("Figures/aav2_blood_gosh_u10.png", units="in", width=9, height=5, res=300)
# aav2_blood_gosh_plot_u10
# dev.off()

aav2_blood_pd_plot_u10 <- ggplot2::ggplot(data=aav2_blood_pd_u10, aes(x=factor(control_type), y=aav2_ct_inv)) +
  geom_boxplot(width = 0.3, outlier.shape = NA, show.legend = FALSE) +
  geom_point(position=position_jitter(w=0.2, h=0),
             shape = 16, alpha = 0.8, size = 2) +
  theme_light() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 16, angle = 45, hjust = 1),
        axis.title.y = element_text(angle = 0, size=18),
        axis.text.y = element_text(size = 16),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border=element_blank(),
        axis.line = element_line(),
        axis.ticks.y = element_line(color = "black")) +
  ylab("1/Ct") +
  geom_hline(yintercept = 1/thres, linetype = "dashed") +
  annotate("text", x = 0.7, y = 1/thres + 0.002, label = "LLP", size = 5) +
  geom_hline(yintercept = 1/not_detected, linetype = "dashed") +
  # labs(caption = "* p < 0.05, *** p < 0.001\nLLP: Low level positive") +
  stat_compare_means(method = "wilcox.test", label = "p.signif", ref.group = "Case",
                     label.y = 1/17, size = 6, 
                     symnum.args = list(cutpoints = c(0, 0.00005, 0.001, 0.01, 1), 
                                        symbols = c("***", "**", "*", "NS"))) +
  scale_y_continuous(breaks = c(1/50, 1/40, 1/35, 1/30, 1/25, 1/20, 1/15), 
                     labels = c("ND", "1/40", "1/35", "1/30", "1/25", "1/20", "1/15"),
                     limits=c(1/51, 1/14))

# png("Figures/aav2_blood_pd_u10.png", units="in", width=8, height=5, res=300)
# aav2_blood_pd_plot_u10
# dev.off()

# AdV in blood
adeno_blood_plot <- ggplot2::ggplot(data=adeno_blood, aes(x=factor(control_type), y=adeno_ct_inv)) +
  geom_boxplot(width = 0.3, outlier.shape = NA, show.legend = FALSE) +
  geom_point(position=position_jitter(w=0.2, h=0), shape = 16, alpha = 0.8, size = 2) +
  theme_light() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 16, angle = 45, hjust = 1),
        axis.title.y = element_text(angle = 0, size=18),
        axis.text.y = element_text(size = 16),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border=element_blank(),
        axis.line = element_line(),
        axis.ticks.x = element_line(color = "black"),
        axis.ticks.y = element_line(color = "black")) +
  ylab("1/Ct") +
  geom_hline(yintercept = 1/thres, linetype = "dashed") +
  annotate("text", x = 8.2, y = 1/thres + 0.002, label = "LLP", size = 5) +
  geom_hline(yintercept = 1/not_detected, linetype = "dashed") +
  scale_y_continuous(breaks = c(1/50, 1/40, 1/35, 1/30, 1/25, 1/20, 1/15), 
                     labels = c("ND", "1/40", "1/35", "1/30", "1/25", "1/20", "1/15"),
                     limits=c(1/51, 1/14)) +
  stat_compare_means(method = "wilcox.test", label = "p.signif", ref.group = "Case",
                     label.y = 1/17, size = 6,
                     # thresholds adjusted manually since the package doesn't incorporate multiple hypothesis correction
                     symnum.args = list(cutpoints = c(0, 0.00001, 0.003, 0.02, 1), 
                                        symbols = c("***", "**", "*", "NS")))


# png("Figures/adeno_blood_all.png", units="in", width=9, height=5, res=300)
# adeno_blood_plot
# dev.off()

hhv6_blood_plot <- ggplot2::ggplot(data=hhv6_blood, aes(x=factor(control_type), y=hhv6_ct_inv)) +
  geom_boxplot(width = 0.3, outlier.shape = NA, show.legend = FALSE) +
  geom_point(position=position_jitter(w=0.2, h=0),
             shape = 16, alpha = 0.8, size = 2) +
  theme_light() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 16, angle = 45, hjust = 1),
        axis.title.y = element_text(angle = 0, size=18),
        axis.text.y = element_text(size = 16),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border=element_blank(),
        axis.line = element_line(),
        axis.ticks.y = element_line(color = "black")) +
  ylab("1/Ct") +
  geom_hline(yintercept = 1/thres, linetype = "dashed") +
  annotate("text", x = 0.65, y = 1/thres + 0.002, label = "LLP", size = 5) +
  geom_hline(yintercept = 1/not_detected, linetype = "dashed") +
  scale_y_continuous(breaks = c(1/50, 1/40, 1/35, 1/30, 1/25, 1/20, 1/15), 
                     labels = c("ND", "1/40", "1/35", "1/30", "1/25", "1/20", "1/15"),
                     limits=c(1/51, 1/14)) +
  stat_compare_means(method = "wilcox.test", label = "p.signif", ref.group = "Case",
                     label.y = 1/17, size = 6,
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.002, 0.01, 1),
                                        symbols = c("***", "**", "*", "NS")))

# png("Figures/hhv6_blood_all.png", units="in", width=8, height=5, res=300)
# hhv6_blood_plot
# dev.off()

#### Stacked bar chart for comparison of positive and negative results ####

# Prepare data
data_blood <- dplyr::filter(data, tissue == "blood")

bar_blood_adeno <- janitor::tabyl(data_blood, control_type, adeno_pcr) %>% 
  tidyr::pivot_longer(!control_type, names_to = "result", values_to = "count") %>% 
  dplyr::mutate(comparison = "AdV Blood") %>% 
  dplyr::rename(group = control_type)

bar_blood_aav2 <- janitor::tabyl(data_blood, control_type, aav2_pcr) %>% 
  tidyr::pivot_longer(!control_type, names_to = "result", values_to = "count") %>% 
  dplyr::mutate(comparison = "AAV2 Blood") %>% 
  dplyr::rename(group = control_type)

bar_blood_hhv6 <- janitor::tabyl(data_blood, control_type, hhv6_pcr) %>% 
  tidyr::pivot_longer(!control_type, names_to = "result", values_to = "count") %>% 
  dplyr::mutate(comparison = "HHV6 Blood") %>% 
  dplyr::rename(group = control_type)

data_liver <- data %>% 
  dplyr::filter(tissue == "liver")

bar_liver_adeno <- janitor::tabyl(data_liver, group, adeno_pcr) %>% 
  tidyr::pivot_longer(!group, names_to = "result", values_to = "count") %>% 
  dplyr::mutate(comparison = "AdV Liver") %>% 
  dplyr::mutate(group = ifelse(group == "control", "control_liver", group))

bar_liver_aav2 <- janitor::tabyl(data_liver, group, aav2_pcr) %>% 
  tidyr::pivot_longer(!group, names_to = "result", values_to = "count") %>% 
  dplyr::mutate(comparison = "AAV2 Liver") %>% 
  dplyr::mutate(group = ifelse(group == "control", "control_liver", group))

bar_liver_hhv6 <- janitor::tabyl(data_liver, group, hhv6_pcr) %>% 
  tidyr::pivot_longer(!group, names_to = "result", values_to = "count") %>% 
  dplyr::mutate(comparison = "HHV6 Liver") %>%
  dplyr::mutate(group = ifelse(group == "control", "control_liver", group))

bar_aav2 <- dplyr::full_join(bar_blood_aav2, bar_liver_aav2) %>% 
  dplyr::filter(!is.na(result) & result != "NA_" & group %in% control_types_aav2 & count != 0) %>% 
  dplr::mutate(comparison = factor(comparison, levels = c("AAV2 Liver", "AAV2 Blood", "AAV2 Stool"),
                             labels = c("AAV2\nLiver", "AAV2\nBlood", "AAV2\nStool")),
               result = factor(result, levels = c("negative", "LLP", "positive"), 
                         labels = c("Negative", "Low-level positive", "Positive")),
               group = factor(group, levels = control_types_aav2, labels = control_types_aav2_labels))

bar_adeno <- dplyr::full_join(bar_blood_adeno, bar_liver_adeno) %>% 
  dplyr::filter(!is.na(result) & result != "NA_" & group %in% control_types & count != 0) %>% 
  dplyr::mutate(comparison = factor(comparison, levels = c("AdV Liver", "AdV Blood", "AdV Stool"),
                             labels = c("AdV\nLiver", "AdV\nBlood", "AdV\nStool")), 
                result = factor(result, levels = c("negative", "LLP", "positive"), 
                         labels = c("Negative", "Low-level positive", "Positive")), 
                group = factor(group, levels = control_types, labels = control_types_labels))

bar_hhv6 <- dplyr::full_join(bar_blood_hhv6, bar_liver_hhv6) %>% 
  dplyr::filter(!is.na(result) & result != "NA_" & group %in% control_types & count != 0) %>% 
  dplyr::mutate(comparison = factor(comparison, levels = c("HHV6 Liver", "HHV6 Blood", "HHV6 Stool"),
                             labels = c("HHV6\nLiver", "HHV6\nBlood", "HHV6\nStool")),
                result = factor(result, levels = c("negative", "LLP", "positive"), 
                         labels = c("Negative", "Low-level positive", "Positive")),
                group = factor(group, levels = control_types, labels = control_types_labels))

# AAV2 plot
pcr_results_aav2 <- ggplot2::ggplot(data = bar_aav2, aes(x = group, count, fill = result)) +
  geom_col(position = "stack") +
  facet_grid(cols = vars(comparison), scales = "free_x", space = "free_x") +
  theme_light() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.y.right = element_text(angle = 0),
        legend.position = "top") +
  scale_y_continuous(expand = c(0,0), limits= c(0,25)) +
  ylab("Number of Samples") +
  scale_fill_manual(name = "PCR result", 
                    = c("firebrick", "lightblue", "dodgerblue4"), 
                    guide=guide_legend(reverse=T))

# png("Figures/pcr_results_aav2.png", units="in", width=7, height=4, res=300)
# pcr_results_aav2
# dev.off()

# AdV plot
pcr_results_adeno <- ggplot2::ggplot(data = bar_adeno, aes(x = group, count, fill = result)) +
  geom_col(position = "stack") +
  facet_grid(cols = vars(comparison), scales = "free", space = "free") +
  theme_light() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(angle = 0),
        legend.position = "top") +
  scale_y_continuous(expand = c(0,0), limits= c(0,25)) +
  ylab("Number of Samples") +
  scale_fill_manual(name = "PCR result", 
                    values = c("firebrick", "lightblue", "dodgerblue4"), 
                    guide=guide_legend(reverse=T))

# png("Figures/pcr_results_adeno.png", units="in", width=6, height=4, res=300)
# pcr_results_adeno
# dev.off()

# HHV6 plot
pcr_results_hhv6 <- ggplot2::ggplot(data = bar_hhv6, aes(x = group, count, fill = result)) +
  geom_col(position = "stack") +
  facet_grid(cols = vars(comparison), scales = "free_x", space = "free_x") +
  theme_light() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.y.right = element_text(angle = 0),
        legend.position = "top") +
  scale_y_continuous(expand = c(0,0), limits= c(0,25)) +
  ylab("Number of Samples") +
  scale_fill_manual(name = "PCR result", 
                    values = c("firebrick", "dodgerblue4", "lightblue"), 
                    guide=guide_legend(reverse=T))

# png("Figures/pcr_results_hhv6.png", units="in", width=6, height=4, res=300)
# pcr_results_hhv6
# dev.off()





  
