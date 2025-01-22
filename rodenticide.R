library(tidyverse)
library(dplyr)
library(ggthemes)
library(conflicted)

conflicts_prefer(stats::fisher.test)
conflicts_prefer(stats::chisq.test)

setwd("/Users/trent/Desktop/CSU/rodenticide paper/")
data <- read.csv("liversamples_clean.csv")
lookup <- read.csv("lookup.csv")
generations <- read.csv("generations.csv")

generations <- generations %>%
  mutate(number_livers = as.numeric(number_livers))

generations$percentage <- generations$number_livers/28

exposure <- read.csv("exposure.csv")

#figure 2
ggplot(data)+
  geom_bar(aes(x=total_detected, fill = factor(species)), position = position_dodge2(width = 0.7, preserve = "single"))+
  geom_text(aes(x = total_detected, label = ..count.., group = species), stat = 'count', position = position_dodge2(width = 0.7, preserve = "single"),
            vjust = -0.3, size = 3) +
  labs(x = "Total anticoagulant rodenticides detected per sample", y = "Liver samples", title="", fill="Species")+
  scale_fill_manual(values = c("gray60", "tan"))+
  scale_x_continuous(breaks = unique(data$total_detected), labels = unique(data$total_detected))+
  scale_y_continuous(breaks=seq(0,12,2), expand = c(0, 0), limits = c(0, 8))+
  theme_classic()+
  theme(legend.position = "right", 
        legend.background = element_rect(color = "black", size = 0.5),
        legend.text = element_text(size = 7),
        plot.margin = margin(1, 1, 1, 1, "cm"),
        plot.title = element_text(hjust = 0.5, vjust = 5, size = 15),
        strip.text = element_text(face = "italic", size = 10),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 9))

#ggsave("figure2.pdf", plot = last_plot(), width = 8, height = 6, path = "/Users/trent/Desktop/CSU/rodenticide paper/final_plots")

#figure 3
exposure_long <- exposure %>%
  pivot_longer(X1L2022:X28L2022, names_to = "sample_id", values_to = "ppm")%>%
  mutate(sample_id = str_remove_all(sample_id, "X"))%>%
  left_join(data %>% select(sample_id, species), by = "sample_id") %>%
  mutate(generation = as.factor(generation))

plotdata <- exposure_long %>%
  mutate(ppm = ifelse(is.na(ppm), 0, ppm)) %>%
  group_by(species, generation, rodenticide) %>%
  summarize(n = sum(ppm > 0) + sum(ppm == 0 & is.na(ppm))) %>%
  ungroup()


ggplot(plotdata, aes(x = factor(rodenticide), y = n, fill = species)) +
  geom_col() +
  #geom_text(aes(label = n),
  #vjust = -0.3, color = "black", size = 3) +
  scale_fill_manual(values = c("gray60", "tan"))+
  scale_y_continuous(breaks = seq(0, 10, 2), expand = c(0, 0), limits = c(0, 10))+
  labs(
    x = "",
    y = "Number of livers exposed to anticoagulant rodenticides",
    fill = "Species",
    title = "") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 15),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.title.y = element_text(size = 10),
        legend.position = "right",
        legend.box.background = element_rect(color = "black", size = 1))+
  facet_wrap(~generation, scales = "free_x", strip.position = "bottom")

#ggsave("figure3.pdf", plot = last_plot(), width = 8, height = 6, path = "/Users/trent/Desktop/CSU/rodenticide paper/final_plots")

#exposure test 1
exposure_test <- data %>%
  select(sample_id, 11:18)%>%
  mutate(across(.cols = 2:9, 
                .fns = ~ifelse(. == "not detected", 0, 1)))%>%
  rowwise() %>%
  mutate(first_gen = sum(c(chlor_ppm, dic_ppm, dipha_ppm, warf_ppm), na.rm = TRUE),
         first_gen = ifelse(first_gen > 0, 1, 0)) %>%
  ungroup() %>%
  rowwise() %>%
  mutate(second_gen = sum(c(brod_ppm, brom_ppm, difen_ppm, difet_ppm), na.rm = TRUE),
         second_gen = ifelse(second_gen > 0, 1, 0)) %>%
  ungroup() %>%
  select(first_gen, second_gen) %>%
  group_by(first_gen, second_gen) %>%
  summarise(Count = n(), .groups = 'drop') %>%
  pivot_wider(names_from = second_gen, values_from = Count, values_fill = list(Count = 0)) %>%
  select(-first_gen) %>%
  as.matrix()
dimnames(exposure_test) <- list(first_gen = c("No", "Yes"), second_gen = c("No", "Yes"))

exposure_test
fisher.test(exposure_test)
chisq.test(exposure_test)

#exposure test 2
#redundant code, for error checking
exposure_test2 <- data %>%
  select(sample_id, 11:20)%>%
  mutate(across(.cols = 2:9, 
                .fns = ~ifelse(. == "not detected", 0, 1)))%>%
  rowwise() %>%
  mutate(first_gen = sum(c(chlor_ppm, dic_ppm, dipha_ppm, warf_ppm), na.rm = TRUE),
         first_gen = ifelse(first_gen > 0, 1, 0)) %>%
  ungroup() %>%
  rowwise() %>%
  mutate(second_gen = sum(c(brod_ppm, brom_ppm, difen_ppm, difet_ppm), na.rm = TRUE),
         second_gen = ifelse(second_gen > 0, 1, 0)) %>%
  ungroup()%>%
  select(sample_id, first_gen, second_gen, total_detected)

exposure_test2_table <- table(exposure_test2$first_gen, exposure_test2$second_gen)
exposure_test2_table

fisher.test(exposure_test2_table)
chisq.test(exposure_test2_table, correct = FALSE)

