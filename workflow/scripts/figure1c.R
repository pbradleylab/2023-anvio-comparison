#Meant to be run after figure1.R!
library(png)
library(cowplot)

# Condense Duplicate Values to make plotting easier since we really don't need
# to plot values that are exactly the same.
strays_anst_unq <- strays_anst %>% distinct(V5, .keep_all = TRUE)
not_standard_mare_unq <- unique(subset(mare, as.numeric(V6) > 0.01)) %>% distinct(V6, .keep_all = TRUE)
strays_mare_unq <- strays_mare %>% distinct(V6, .keep_all = TRUE)
anst_unq <- anst %>% distinct(V5, .keep_all = TRUE)
anra_unq <- anra %>% distinct(V5, .keep_all = TRUE)
mare_unq <- mare %>% distinct(V6, .keep_all = TRUE)

# Generate evals dataframes to compare strays
colnames(anst_unq)[colnames(anst_unq) == "V5"] <- "Anvio (Adaptive+Stray)"
colnames(anra_unq)[colnames(anra_unq) == "V5"] <- "Anvio (Adaptive)"
colnames(mare_unq)[colnames(mare_unq) == "V6"] <- "MicrobeAnnotator"
colnames(strays_anst_unq)[colnames(strays_anst_unq) == "V5"] <- "Anvio (Adaptive+Stray)"
colnames(strays_mare_unq)[colnames(strays_mare_unq) == "V6"] <- "MicrobeAnnotator"
colnames(not_standard_mare_unq)[colnames(not_standard_mare_unq) == "V6"] <- "MicrobeAnnotator"

mare_unq$MicrobeAnnotator <- as.numeric(mare_unq$MicrobeAnnotator)
strays_mare_unq$MicrobeAnnotator <- as.numeric(strays_mare_unq$MicrobeAnnotator)
not_standard_mare_unq$MicrobeAnnotator <- as.numeric(not_standard_mare_unq$MicrobeAnnotator)

evals <- bind_rows(dplyr::select(anst_unq, `Anvio (Adaptive+Stray)`, `Anvio (Adaptive+Stray)`),
                   dplyr::select(anra_unq, `Anvio (Adaptive)`, `Anvio (Adaptive)`),
                   dplyr::select(mare_unq, MicrobeAnnotator, MicrobeAnnotator),
                   .id = "Source") %>% 
  pivot_longer(cols = c(`Anvio (Adaptive+Stray)`, `Anvio (Adaptive)`, MicrobeAnnotator)) %>%
  drop_na()

evals_stray <- bind_rows(dplyr::select(strays_anst_unq, `Anvio (Adaptive+Stray)`, `Anvio (Adaptive+Stray)`),
                         dplyr::select(not_standard_mare_unq, MicrobeAnnotator, MicrobeAnnotator), 
                         .id = "Source") %>% 
  pivot_longer(cols = c(`Anvio (Adaptive+Stray)`, MicrobeAnnotator)) %>%
  drop_na()

evals_subset <- evals %>% 
  filter(name %in% c("Anvio (Adaptive+Stray)", "Anvio (Adaptive)"))
evals_stray_subset <- evals_stray %>% 
  filter(name %in% c("Anvio (Adaptive+Stray)", "Anvio (Adaptive)"))

#evals_stray_subset <- rbind(evals_stray_subset, c(0, "Anvio (Adaptive)", 0.0))

# Combine the datasets and add a column to distinguish between them
evals_combined <- rbind(
  mutate(evals_subset, set = "Subset"),
  mutate(evals_stray_subset, set = "Stray")
)

# Define custom colors
custom_colors <- c("#2a9d8f", "#2a9d8f", "#6DB9E4", "black")

# Plotting  
p1 <- ggplot() +
  geom_jitter(data = evals, aes(x = name, y = value, color = name),
              position = position_jitter(width = 0.2), alpha = 1) +
  geom_jitter(data = evals_stray, aes(x = name, y = value, color = "Stray"),
              position = position_jitter(width = 0.2), alpha = 1) +
  scale_color_manual(values = custom_colors, 
                     labels = c("Anvio", "Anvio", "MicrobeAnnotator", "Stray")) + 
  labs(title = "Stray Recovery Methods E-Values",
       x = "Category",
       y = "E-Value",
       color = "Method")

# Create the plot
p2 <- ggplot(data = evals_combined, aes(x = name, y = value, color = set)) +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.6) +
  scale_color_manual(values = c("Subset" = custom_colors[1], "Stray" = "black")) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    plot.background = element_rect(color = "black", fill = NA, size = 1, linetype = "dotted"), # Dotted border around entire plot
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  ) +
  labs(x = "Category",
       y = "E-Value") 

c <- p1 + 
  annotation_custom(ggplotGrob(p2), xmin = .5, xmax = 2.5, 
                    ymin = 2.5, ymax = 10) +
  geom_rect(aes(xmin = .5, xmax = 2.5, ymin = 2.5, ymax = 10), 
            color='black', linetype='dashed', alpha=0) +
  geom_path(aes(x,y,group=grp), 
            data=data.frame(x = c(1,.5,2,2.5), y=c(0.3,2.5,0.3,2.5),grp=c(1,1,2,2)),
            linetype='dashed')

first_row = plot_grid(a, b, labels = c('A','B'), nrow=1)
second_row = plot_grid(c, labels = c('C'))
gg_all = plot_grid(first_row, second_row, labels=c('', ''), ncol=1)



# Make multiple plots in same pane
#combined_plot <- plot_grid(a, c, b, ncol = 2, rel_heights = c(2, 1), align = "v")

