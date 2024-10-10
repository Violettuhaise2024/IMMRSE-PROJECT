###################################################
###          Install required packages          ###
###################################################
library(readxl)
library(janitor)
library(tidyverse)
library(ggplot2)
library(jpeg)

###################################################
###                 Import Data                 ###
###################################################
rm(list=ls())
detach(package:plotly, unload=TRUE)
setwd("D:/Computer/IMMRSE/Data-from-pipeline/Drug_resistance/Round3/")

data <- read.csv("databases/sample_genotypes_R3.csv")
#data_meta <- bind_rows(data_meta_HTS, data_meta_LTS)
#metadata <- read_csv("./")

###################################################
###                     Plots                   ###
###################################################
##### TABULATE Ns #####
# N_collected <- as.data.frame(table(metadata$site))
# names(N_collected) <- c("site", "N_collected")
# 
# N_genotyped <- as.data.frame(table(data_meta$site))
# names(N_genotyped) <- c("site", "N_genotyped")
# 
# N_na <- rowSums(is.na(data_meta[2:84]))
# N_na2 <- as.data.frame(cbind(data_meta$barcode, data_meta$site, N_na))
# N_filt <- subset(N_na2, N_na < 40)
# names(N_filt) <- c("barcode","site", "N_na")
# N_40 <- as.data.frame(table(N_filt$site))
# names(N_40) <- c("site", "N_40percentSuccess")
# 
# tabs1 <- merge(N_collected, N_genotyped, by = "site", all.x = TRUE)
# tabs <- merge(tabs1, N_40, by = "site", all.x = TRUE)
# 
# write.csv(tabs, "./231229/databases/Sample_sizes_combined_231229.csv", row.names = FALSE)

##### TABULATE Ns #####
meta_long <- pivot_longer(data, cols = arps10_V127M:dhps_S436C) 
#names(meta_long) <- c("SampleName", "site", "monthyear", "agecategories", "sex", "malariatestdone", "malariatestresults", "fpcollected", "SNP", "Genotype")

names(meta_long) <- c("SampleName", "site", "SNP", "Genotype")

tab_meta <- as.data.frame(xtabs(~SNP + Genotype + site, meta_long))
tab_meta_wide <- spread(tab_meta, Genotype, Freq)
names(tab_meta_wide) <- c("SNP", "Site", "WT_N", "Mix_N", "Mut_N")
tab_meta_wide$N <- tab_meta_wide$WT_N + tab_meta_wide$Mix_N + tab_meta_wide$Mut_N
tab_meta_wide$WT_f <- round(tab_meta_wide$WT_N/tab_meta_wide$N * 100, 1)
tab_meta_wide$Mix_f <- round(tab_meta_wide$Mix_N/tab_meta_wide$N * 100, 1)
tab_meta_wide$Mut_f <- round(tab_meta_wide$Mut_N/tab_meta_wide$N * 100, 1)
tab_meta_wide <- subset(tab_meta_wide, select = c("SNP", "Site", "N", "WT_N", "Mix_N", "Mut_N", "WT_f", "Mix_f", "Mut_f"))
tab_meta_wide <- tab_meta_wide %>% mutate(District = case_when(Site=="Aboke" ~ "Kole",
                                                               Site=="Aduku" ~ "Kwania",
                                                               Site=="Alebtong" ~ "Alebtong",
                                                               Site=="Amolatar" ~ "Amolatar",
                                                               Site=="Atiak" ~ "Amuru",
                                                               Site=="Bikurungu" ~ "Rukungiri",
                                                               Site=="Budondo" ~ "Jinja",
                                                               Site=="Busitema" ~ "Busia",
                                                               Site=="Chahafi" ~ "Kisoro",
                                                               Site=="Kamwezi" ~ "Rukiga",
                                                               Site=="Karambi" ~ "Kasese",
                                                               Site=="Kasambya" ~ "Mubende",
                                                               Site=="Kigandalo" ~ "Mayuge",
                                                               Site=="Kigorobya" ~ "Hoima",
                                                               Site=="Kihihi" ~ "Kanungu",
                                                               Site=="Kiyunga" ~ "Mukono",
                                                               Site=="Lalogi" ~ "Omoro",
                                                               Site=="Lobule" ~ "Koboko",
                                                               Site=="Lokolia" ~ "Kaabong",
                                                               Site=="Maziba" ~ "Kabale",
                                                               Site=="Metu" ~ "Moyo",
                                                               Site=="Muko" ~ "Rubanda",
                                                               Site=="Nadunget" ~ "Moroto",
                                                               Site=="Nagongera" ~ "Tororo",
                                                               Site=="Namokora" ~ "Kitgum",
                                                               Site=="Nawaikoke" ~ "Kaliro",
                                                               Site=="Orum" ~ "Otuke",
                                                               Site=="Padibe" ~ "Lamwo",
                                                               Site=="Patongo" ~ "Agago",
                                                               Site=="Opia" ~ "Arua",
                                                               Site=="Awach" ~ "Gulu",
                                                               Site=="Ayipe" ~ "Koboko",
                                                               Site=="Diima" ~ "Kiryandongo",
                                                               Site=="Kyatiri" ~ "Masindi",
                                                               TRUE ~ NA))

tab_meta_long <- pivot_longer(tab_meta_wide, cols = WT_f:Mut_f, names_to = "genotype", values_to = "prev")

tab_meta_long <- tab_meta_long %>% 
  mutate(Region = case_when(
    District %in% c("Kole", "Kwania", "Alebtong", "Amolatar", "Amuru", "Omoro", "Kitgum", "Otuke", "Lamwo", "Agago", "Gulu") ~ "North",
    District %in% c("Jinja", "Busia", "Mayuge", "Tororo", "Kaliro", "Mukono") ~ "East",
    District %in% c("Rukungiri", "Kasese", "Hoima", "Kanungu", "Mubende", "Kiryandongo", "Masindi") ~ "West",
    District %in% c("Koboko", "Moyo", "Arua") ~ "North West",
    District %in% c("Kaabong", "Moroto") ~ "North East",
    District %in% c("Mubende", "Mukono") ~ "Central",
    District %in% c("Kisoro", "Rukiga", "Kabale", "Rubanda", "Kanungu") ~ "South West",
    TRUE ~ NA_character_
  ))



write.csv(tab_meta_wide, "./databases/IMMRSE-U_tabulations_R3.csv", row.names = FALSE)

tab_meta_long$Site <- factor(tab_meta_long$Site, levels = unique(tab_meta_long$Site[order(tab_meta_long$Region)]))

############## All loci
#plot_tri <- 
  ggplot(tab_meta_long, aes(x= Site, y = prev, fill = factor(genotype, levels=c("WT_f", "Mix_f", "Mut_f")))) +
  geom_col() +
  facet_wrap(facets = tab_meta_long$SNP) +
  scale_fill_manual(name="Genotype", values = c("gray", "#F8766D", "#00BFC4"), labels = c("% WT", "% Mix", "% Mut")) +
  ylab("Prevalence (%)") +
  xlab("Site") +
  scale_x_discrete(guide = guide_axis(angle = 90))
#plot_tri


############## All K13
k13 <- subset(tab_meta_long, str_detect(SNP, "k13_"))
#plot_all_k13 <-
  ggplot(k13, aes(x = Site, y = prev, fill = factor(genotype, levels=c("WT_f", "Mix_f", "Mut_f")))) +
  geom_col() +
  facet_wrap(facets = k13$SNP) +
  scale_fill_manual(name="Genotype", values = c("gray", "#F8766D", "#00BFC4"), labels = c("% WT", "% Mix", "% Mut")) +
  ylab("Prevalence (%)") +
  xlab("Site") +
  scale_x_discrete(guide = guide_axis(angle = 90))
#plot_all_k13

############## All Polymorphic K13
k13_variable <- subset(tab_meta_long, SNP=="k13_R515K" | SNP=="k13_C469Y" | SNP=="k13_C469F" | SNP=="k13_R561H" | SNP=="k13_P441A" | SNP=="k13_P441L" | SNP=="k13_A578S" | SNP=="k13_C580Y" | SNP=="k13_A675V" | SNP== "k13_F442L")

plot_k13_variable <-
  ggplot(k13_variable, aes(x = Site, y = prev, fill = factor(genotype, levels=c("WT_f", "Mix_f", "Mut_f")))) +
  geom_col() +
  facet_wrap(facets = k13_variable$SNP) +
  scale_fill_manual(name="Genotype", values = c("gray", "#F8766D", "#00BFC4"), labels = c("% WT", "% Mix", "% Mut")) +
  ylab("Prevalence (%)") +
  xlab("Site") +
  scale_x_discrete(guide = guide_axis(angle = 90))
plot_k13_variable
ggsave(plot = plot_k13_variable, "./figures/k13_polymorphic_loci.svg", width = 15, height = 7.5)

############## K13 mutations of interest
k13_ofinterest <- subset(tab_meta_long, (SNP=="k13_C469Y" | SNP=="k13_C469F" | SNP=="k13_P441L" | SNP=="k13_P441A" | SNP=="k13_R561H" | SNP=="k13_R622I" | SNP=="k13_A675V" | SNP=="k13_C580Y"))
  
plot_k13_ofinterest <-
  ggplot(k13_ofinterest, aes(x = Site, y = prev, fill = factor(genotype, levels=c("WT_f", "Mix_f", "Mut_f")))) +
    geom_col() +
    facet_wrap(facets = factor(k13_ofinterest$SNP, levels=c("k13_A675V", 
                                                            "k13_C469Y", 
                                                            "k13_C469F", 
                                                            "k13_R561H", 
                                                            "k13_P441L",
                                                            "k13_P441A",
                                                            "k13_C580Y",
                                                            "k13_R622I")), 
                                                  labeller = as_labeller(c(k13_C469Y = "C469Y",
                                                                           k13_A675V = "A675V",
                                                                           k13_C469F = "C469F", 
                                                                           k13_675 = "A675V",
                                                                           k13_R561H = "R561H", 
                                                                           k13_P441L = "P441L",
                                                                           k13_P441A = "P441A",
                                                                           k13_C580Y = "C580Y",
                                                                           k13_R622I = "R622I"))) +
    scale_fill_manual(name="Genotype", values = c("gray", "#F8766D", "#00BFC4"), labels = c("% WT", "% Mix", "% Mut")) +
    ylab("Prevalence (%)") +
    xlab("Site") +
    scale_x_discrete(guide = guide_axis(angle = 90))
plot_k13_ofinterest
ggsave(plot = plot_k13_ofinterest, "./figures/k13_loi.svg", width = 15, height = 7.5)


############## K13 mutations of interest
k13_ofinterest1 <- subset(tab_meta_long, (SNP=="k13_C469Y" | SNP=="k13_C469F" | SNP=="k13_P441L" | SNP=="k13_P441A" | SNP=="k13_R561H" | SNP=="k13_C580Y" | SNP=="k13_A675V"))

plot_k13_ofinterest1 <-
  ggplot(k13_ofinterest1, aes(x = Site, y = prev, fill = factor(genotype, levels=c("WT_f", "Mix_f", "Mut_f")))) +
  geom_col() +
  facet_wrap(facets = factor(k13_ofinterest1$SNP, levels=c("k13_A675V", 
                                                          "k13_C469Y", 
                                                          "k13_C469F", 
                                                          "k13_R561H", 
                                                          "k13_P441L",
                                                          "k13_P441A",
                                                          "k13_C580Y")), 
             labeller = as_labeller(c(k13_C469Y = "C469Y", 
                                      k13_C469F = "C469F", 
                                      k13_A675V = "A675V",
                                      k13_R561H = "R561H", 
                                      k13_P441L = "P441L",
                                      k13_P441A = "P441A",
                                      k13_C580Y = "C580Y"))) +
  scale_fill_manual(name="Genotype", values = c("gray", "#F8766D", "#00BFC4"), labels = c("% WT", "% Mix", "% Mut")) +
  ylab("Prevalence (%)") +
  xlab("Site") +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  coord_cartesian(ylim = c(0, 50))
plot_k13_ofinterest1
ggsave(plot = plot_k13_ofinterest1, "./figures/k13_loi_axis50.svg", height = 15, width = 7.5)


############## K13 mutations of interest 50% axis
k13_ofinterest2 <- subset(tab_meta_long, (SNP=="k13_C469Y" | SNP=="k13_R622I" | SNP=="k13_P441L" | SNP=="k13_R561H" | SNP=="k13_C580Y" | SNP=="k13_A675V") & (genotype=="Mix_f" | genotype=="Mut_f"))

jpeg("./figures/2024_R3_pfK13_Interest_prev.jpeg", width = 25, height = 15, units = "cm", res = 300, pointsize = 8)
  ggplot(k13_ofinterest2, aes(x = Site, y = prev, fill = factor(genotype, levels=c("Mix_f", "Mut_f")))) +
  geom_col() +
  facet_wrap(facets = factor(k13_ofinterest2$SNP, levels=c("k13_A675V", 
                                                          "k13_C469Y", 
                                                          "k13_R622I", 
                                                          "k13_R561H", 
                                                          "k13_P441L",
                                                          "k13_C580Y")), 
             labeller = as_labeller(c(k13_C469Y = "C469Y", 
                                      k13_R622I = "R622I", 
                                      k13_A675V = "A675V",
                                      k13_R561H = "R561H", 
                                      k13_P441L = "P441L",
                                      k13_C580Y = "C580Y"))) +
  scale_fill_manual(name="Genotype", values = c("#F8766D", "#00BFC4"), labels = c("% Mix", "% Mut")) +
  ylab("Prevalence (%)") +
  xlab("Site") +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  coord_cartesian(ylim = c(0, 50))
dev.off()

# plot_k13_ofinterest_v2
# ggsave(plot = plot_k13_ofinterest_v2, "./figures/k13_loi_axis50_no441A.svg", height = 15, width = 7.5)


#Katairo #####K13 Second graph
k13_secondgraph <- subset(tab_meta_long, (SNP=="k13_R515K" | SNP=="k13_C580Y" | SNP=="k13_R539T" | SNP=="k13_A578S" | SNP== "k13_F442L") )

jpeg("./figures/2024_R3_pfK13_2_prev.jpeg", width = 20, height = 15, units = "cm", res = 300, pointsize = 8)

ggplot(k13_secondgraph, aes(x = Site, y = prev, fill = factor(genotype, levels=c("WT_f", "Mix_f", "Mut_f")))) +
  geom_col() +
  facet_wrap(~SNP, ncol = 2, 
             labeller = labeller(SNP = c("k13_R515K" = "R515K",
                                         "k13_C580Y" = "C580Y",
                                         "k13_R539T" = "R539T",
                                         "k13_A578S" = "A578S",
                                         "k13_F442L" = "F442L"))) +
  scale_fill_manual(name="Genotype", values = c("gray", "#F8766D", "#00BFC4"), labels = c("% WT", "% Mix", "% Mut")) +
  theme_bw() +
  ylab("Prevalence (%)") +
  xlab("Site") +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  coord_cartesian(ylim = c(0, 10)) +
  theme(legend.position = "bottom",
        legend.direction = "horizontal") +
  geom_hline(yintercept = 5, linetype = "dotted", color = "black", size = 0.5)

dev.off()


#Katairo ###K13 Second graph

k13_firstgraph <- subset(tab_meta_long, (SNP=="k13_C469Y" | SNP=="k13_C469F" | SNP=="k13_P441L" | SNP=="k13_R561H" | SNP=="k13_A675V" | SNP=="k13_P441A") )

jpeg("./figures/2024_R3_pfK13_1_prev.jpeg", width = 20, height = 15, units = "cm", res = 300, pointsize = 8)

ggplot(k13_firstgraph, aes(x = Site, y = prev, fill = factor(genotype, levels=c("WT_f", "Mix_f", "Mut_f")))) +
  geom_col() +
  facet_wrap(~SNP, ncol = 2, 
             labeller = labeller(SNP = c("k13_A675V" = "A675V",
                                         "k13_C469F" = "C469F",
                                         "k13_C469Y" = "C469Y",
                                         "k13_P441A" = "P441A",
                                         "k13_P441L" = "P441L",
                                         "k13_R561H" = "R561H" ))) +
  scale_fill_manual(name="Genotype", values = c("gray", "#F8766D", "#00BFC4"), labels = c("% WT", "% Mix", "% Mut")) +
  theme_bw() +
  ylab("Prevalence (%)") +
  xlab("Site") +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  coord_cartesian(ylim = c(0, 50)) +
  theme(legend.position = "bottom",
        legend.direction = "horizontal") +
  geom_hline(yintercept = c(20, 40), linetype = "dotted", color = "black", size = 0.5)

dev.off()

############## K13 mutations of interest 50% axis
# k13_ofinterest3 <- subset(tab_meta_long, (SNP=="k13_C469Y" | SNP=="k13_C469F" | SNP=="k13_P441L" | SNP=="k13_R561H" | SNP=="k13_675") & (genotype=="Mix_f" | genotype=="Mut_f"))
# 
# plot_k13_ofinterest_v3 <-
#   ggplot(k13_ofinterest3, aes(x = Site, y = prev, fill = factor(genotype, levels=c("Mix_f", "Mut_f")))) +
#   geom_col() +
#   facet_wrap(facets = factor(k13_ofinterest3$SNP, levels=c("k13_675", 
#                                                            "k13_C469Y", 
#                                                            "k13_C469F", 
#                                                            "k13_R561H", 
#                                                            "k13_P441L")), 
#              labeller = as_labeller(c(k13_C469Y = "C469Y", 
#                                       k13_C469F = "C469F", 
#                                       k13_675 = "A675V",
#                                       k13_R561H = "R561H", 
#                                       k13_P441L = "P441L"))) +
#   scale_fill_manual(name="Genotype", values = c("#F8766D", "#00BFC4"), labels = c("% Mix", "% Mut")) +
#   ylab("Prevalence (%)") +
#   xlab("Site") +
#   scale_x_discrete(guide = guide_axis(angle = 90)) +
#   coord_cartesian(ylim = c(0, 50))
# plot_k13_ofinterest_v3
# ggsave(plot = plot_k13_ofinterest_v3, "./231229/figures/k13_loi_axis50_no580.tiff", height = 15, width = 7.5)
# 
# 
# 
# ############## K13 mutations of interest 50% axis
# k13_ofinterest4 <- subset(tab_meta_long, (SNP=="k13_C469Y" | SNP=="k13_C469F" | SNP=="k13_P441L" | SNP=="k13_R561H" | SNP=="k13_675") )
# 
# jpeg("2023_pfK13_prev.jpeg", width = 20, height = 15, units = "cm", res = 300, pointsize = 8)
# 
# ggplot(k13, aes(x = Site, y = prev, fill = factor(genotype, levels=c("WT_f", "Mix_f", "Mut_f")))) +
#   geom_col() +
#   facet_wrap(~SNP, ncol = 2, 
#              labeller = labeller(SNP = c("k13_675" = "A675V",
#                                          "k13_C469F" = "C469F",
#                                          "k13_C469Y" = "C469Y",
#                                          "k13_P441A" = "P441A",
#                                          "k13_P441L" = "P441L",
#                                          "k13_R561H" = "R561H" ))) +
#   scale_fill_manual(name="Genotype", values = c("gray", "#F8766D", "#00BFC4"), labels = c("% WT", "% Mix", "% Mut")) +
#   theme_bw()+
#   ylab("Prevalence (%)") +
#   xlab("Site") +
#   scale_x_discrete(guide = guide_axis(angle = 90)) +
#   theme(legend.position = "bottom",
#         legend.direction = "horizontal")
# dev.off()
# 
# 
# jpeg("./231229/figures/k13_ofinterest.jpeg", width = 20, height = 15, units = "cm", res = 300, pointsize = 8)
# 
#   ggplot(k13_ofinterest4, aes(x = Site, y = prev, fill = factor(genotype, levels=c("WT_f", "Mix_f", "Mut_f"))))  +
#   geom_col() +
#   facet_wrap(facets = factor(k13_ofinterest4$SNP, levels=c("k13_675", 
#                                                            "k13_C469Y", 
#                                                            "k13_C469F", 
#                                                            "k13_R561H", 
#                                                            "k13_P441L")),
#              ncol = 2,
#              labeller = as_labeller(c(k13_C469Y = "C469Y", 
#                                       k13_C469F = "C469F", 
#                                       k13_675 = "A675V",
#                                       k13_R561H = "R561H", 
#                                       k13_P441L = "P441L"))) +
#   scale_fill_manual(name="Genotype", values = c("gray", "#F8766D", "#00BFC4"), labels = c("% WT", "% Mix", "% Mut")) +
#   ylab("Prevalence (%)") +
#   xlab("Site") +
#   scale_x_discrete(guide = guide_axis(angle = 90)) +
#   coord_cartesian(ylim = c(0, 50)) +
#   theme(legend.position = "bottom",
#           legend.direction = "horizontal")
# dev.off()
#     
# plot_k13_ofinterest_v4
# ggsave(plot = plot_k13_ofinterest_v4, "./231229/figures/k13_loi_axis50_no580_W.svg", height = 7.5, width = 15)


############## DHFR mutations
dhfr <- subset(tab_meta_long, SNP=="dhfr_N51I" | SNP=="dhfr_C59R" | SNP=="dhfr_S108N" | SNP=="dhfr_I164L")

jpeg("./figures/2024_R3_dhfr_prev.jpeg", width = 20, height = 15, units = "cm", res = 300, pointsize = 8)

ggplot(dhfr, aes(x = Site, y = prev, fill = factor(genotype, levels=c("WT_f", "Mix_f", "Mut_f")))) +
  geom_col() +
  facet_wrap(~SNP, ncol = 2, 
             labeller = labeller(SNP = c("dhfr_N51I" = "N51I",
                                         "dhfr_C59R" = "C59R",
                                         "dhfr_S108N" = "S108N",
                                         "dhfr_I164L" = "I164L"))) +
  scale_fill_manual(name="Genotype", values = c("gray", "#F8766D", "#00BFC4"), labels = c("% WT", "% Mix", "% Mut")) +
  theme_bw() +
  ylab("Prevalence (%)") +
  xlab("Site") +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  coord_cartesian(ylim = c(0, 100)) +
  theme(legend.position = "bottom",
        legend.direction = "horizontal") +
  geom_hline(yintercept = c(50, 75), linetype = "dotted", color = "black", size = 0.5)

dev.off()


# plot_dhfr <-
#   ggplot(dhfr, aes(x = Site, y = prev, fill = factor(genotype, levels=c("WT_f", "Mix_f", "Mut_f")))) +
#   geom_col() +
#   facet_wrap(facets = factor(dhfr$SNP, levels=c("dhfr_51", 
#                                                 "dhfr_C59R", 
#                                                 "dhfr_S108N", 
#                                                 "dhfr_164")), 
#              labeller = as_labeller(c(dhfr_51 = "N51I", 
#                                       dhfr_C59R = "C59R", 
#                                       dhfr_S108N = "S108N",
#                                       dhfr_164 = "I164L"))) +
#   scale_fill_manual(name="Genotype", values = c("gray", "#F8766D", "#00BFC4"), labels = c("% WT", "% Mix", "% Mut")) +
#   ylab("Prevalence (%)") +
#   xlab("Site") +
#   scale_x_discrete(guide = guide_axis(angle = 90))
# plot_dhfr
# ggsave(plot = plot_dhfr, "./231229/figures/dhfr_loi.svg", height = 7.5, width = 15)

############## DHPS mutations
dhps <- subset(tab_meta_long, SNP=="dhps_S436A" | SNP=="dhps_S436C" | SNP=="dhps_S436H" | SNP=="dhps_G437A" |SNP=="dhps_K540E" | SNP=="dhps_A581G" | SNP=="dhps_A613S")

jpeg("./figures/2024_R3_dhps_prev.jpeg", width = 20, height = 15, units = "cm", res = 300, pointsize = 8)

ggplot(dhps, aes(x = Site, y = prev, fill = factor(genotype, levels=c("WT_f", "Mix_f", "Mut_f")))) +
  geom_col() +
  facet_wrap(~SNP, ncol = 2, 
             labeller = labeller(SNP = c("dhps_S436A" = "S436A",
                                         "dhps_S436C" = "S436C",
                                         "dhps_S436H" = "S436H",
                                         "dhps_G437A" = "G437A", 
                                         "dhps_K540E" = "K540E", 
                                         "dhps_A581G" = "A581G",
                                         "dhps_A613S" = "A613S"))) +
  scale_fill_manual(name="Genotype", values = c("gray", "#F8766D", "#00BFC4"), labels = c("% WT", "% Mix", "% Mut")) +
  theme_bw() +
  ylab("Prevalence (%)") +
  xlab("Site") +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  coord_cartesian(ylim = c(0, 100)) +
  theme(legend.position = "bottom",
        legend.direction = "horizontal") +
  geom_hline(yintercept = c(50, 75), linetype = "dotted", color = "black", size = 0.5)

dev.off()


plot_dhps <-
  ggplot(dhps, aes(x = Site, y = prev, fill = factor(genotype, levels=c("WT_f", "Mix_f", "Mut_f")))) +
  geom_col()+
  facet_wrap(facets = factor(dhps$SNP, levels=c("dhps_S436A",
                                                "dhps_S436C",
                                                "dhps_S436H",
                                                "dhps_437",
                                                "dhps_540",
                                                "dhps_581",
                                                "dhps_A613S")),
             labeller = as_labeller(c(dhps_S436A = "S436A",
                                      dhps_S436C = "S436C",
                                      dhps_S436H = "S436H",
                                      dhps_437 = "G437A",
                                      dhps_540 = "K540E",
                                      dhps_581 = "A581G",
                                      dhps_A613S = "A613S"))) +
  scale_fill_manual(name="Genotype", values = c("gray", "#F8766D", "#00BFC4"), labels = c("% WT", "% Mix", "% Mut")) +
  ylab("Prevalence (%)") +
  xlab("Site") +
  scale_x_discrete(guide = guide_axis(angle = 90))
plot_dhps
# ggsave(plot = plot_dhps, "./231229/figures/dhps_loi.svg", height = 7.5, width = 15)

############## Transporter mutations
transporters <- subset(tab_meta_long, (str_detect(SNP, "crt") | str_detect(SNP, "mdr1"))  & SNP!="crt_72" & SNP!="crt_73" & SNP!="crt_K76R" & SNP!="mdr1_N86F" & SNP!="mdr1_N1042Y")
plot_trans <-
  ggplot(transporters, aes(x = Site, y = prev, fill = factor(genotype, levels=c("WT_f", "Mix_f", "Mut_f")))) +
  geom_col()+
  facet_wrap(facets = factor(transporters$SNP, levels=c("crt_E75N", 
                                                        "crt_K76T",
                                                        "mdr1_N86Y",
                                                        "mdr1_Y184F",
                                                        "mdr1_S1034C",
                                                        "mdr1_N1042D",
                                                        "mdr1_D1246Y")), 
             labeller = as_labeller(c(crt_E75N= "E75N", 
                                      crt_K76T = "K76T",
                                      mdr1_N86Y = "N86Y",
                                      mdr1_Y184F = "Y184F",
                                      mdr1_S1034C = "S1034C",
                                      mdr1_N1042D = "N1042D",
                                      mdr1_D1246Y = "D1246Y"))) +
  scale_fill_manual(name="Genotype", values = c("gray", "#F8766D", "#00BFC4"), labels = c("% WT", "% Mix", "% Mut")) +
  ylab("Prevalence (%)") +
  xlab("Site") +
  scale_x_discrete(guide = guide_axis(angle = 90))
plot_trans
#ggsave(plot = plot_trans, "./figures/transporters_loi.svg", height = 15, width = 7.5)

#############################
transporters_2 <- subset(tab_meta_long,  SNP=="crt_K76T" | SNP=="mdr1_N86Y" | SNP=="mdr1_Y184F" | SNP=="mdr1_S1034C" | SNP=="mdr1_N1042D" | SNP=="mdr1_D1246Y")
plot_trans2 <-
  ggplot(transporters_2, aes(x = Site, y = prev, fill = factor(genotype, levels=c("WT_f", "Mix_f", "Mut_f")))) +
  geom_col()+
  facet_wrap(facets = factor(transporters_2$SNP, levels=c("crt_K76T",
                                                        "mdr1_N86Y",
                                                        "mdr1_Y184F",
                                                        "mdr1_1034",
                                                        "mdr1_N1042D",
                                                        "mdr1_D1246Y")), 
             labeller = as_labeller(c(crt_K76T = "K76T",
                                      mdr1_N86Y = "N86Y",
                                      mdr1_Y184F = "Y184F",
                                      mdr1_S1034C = "S1034C",
                                      mdr1_N1042D = "N1042D",
                                      mdr1_D1246Y = "D1246Y"))) +
  scale_fill_manual(name="Genotype", values = c("gray", "#F8766D", "#00BFC4"), labels = c("% WT", "% Mix", "% Mut")) +
  ylab("Prevalence (%)") +
  xlab("Site") +
  scale_x_discrete(guide = guide_axis(angle = 90))
plot_trans2
ggsave(plot = plot_trans2, "./figures/transporters_loi_v2.svg", height = 7.5, width = 15)

########################### Katairo###############
transporters_3 <- subset(tab_meta_long,  SNP=="crt_K76T" | SNP=="mdr1_N86Y" | SNP=="mdr1_N1042D" | SNP=="mdr1_D1246Y")

jpeg("./figures/2024_R3_transporter_prev.jpeg", width = 20, height = 15, units = "cm", res = 300, pointsize = 8)

ggplot(transporters_3, aes(x = Site, y = prev, fill = factor(genotype, levels=c("WT_f", "Mix_f", "Mut_f")))) +
  geom_col() +
  facet_wrap(~SNP, ncol = 2, 
             labeller = labeller(SNP = c("crt_K76T" = "K76T",
                                         "mdr1_N86Y" = "N86Y",
                                         "mdr1_N1042D" = "N1042D",
                                         "mdr1_D1246Y" = "D1246Y"))) +
  scale_fill_manual(name="Genotype", values = c("gray", "#F8766D", "#00BFC4"), labels = c("% WT", "% Mix", "% Mut")) +
  theme_bw() +
  ylab("Prevalence (%)") +
  xlab("Site") +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  coord_cartesian(ylim = c(0, 50)) +
  theme(legend.position = "bottom",
        legend.direction = "horizontal") +
  geom_hline(yintercept = c(10, 25), linetype = "dotted", color = "black", size = 0.5)

dev.off()



##############Antifolates###################
antifolates <- subset(tab_meta_long, SNP=="dhps_G437A" | SNP=="dhps_A581G" | SNP=="dhfr_I164L" | SNP=="dhps_S436A")

jpeg("./figures/2024_R3_antifolates_prev.jpeg", width = 20, height = 15, units = "cm", res = 300, pointsize = 8)

ggplot(antifolates, aes(x = Site, y = prev, fill = factor(genotype, levels=c("WT_f", "Mix_f", "Mut_f")))) +
  geom_col() +
  facet_wrap(~SNP, ncol = 2, 
             labeller = labeller(SNP = c("dhps_S436A" = "dhps S436A", 
                                         "dhfr_I164L" = "dhfr I164L",
                                         "dhps_S436H" = "S436H",
                                         "dhps_G437A" = "dhps G437A", 
                                         "dhps_K540E" = "K540E", 
                                         "dhps_A581G" = "A581G",
                                         "dhps_A613S" = "A613S"))) +
  scale_fill_manual(name="Genotype", values = c("gray", "#F8766D", "#00BFC4"), labels = c("% WT", "% Mix", "% Mut")) +
  theme_bw() +
  ylab("Prevalence (%)") +
  xlab("Site") +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  coord_cartesian(ylim = c(0, 75)) +
  theme(legend.position = "bottom",
        legend.direction = "horizontal") +
  geom_hline(yintercept = c(25, 50, 75), linetype = "dotted", color = "black", size = 0.5)

dev.off()
