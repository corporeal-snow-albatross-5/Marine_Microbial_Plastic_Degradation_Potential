# ============================================================
# Environmental distribution heatmap
# FULL DATASET from Table S1
# ============================================================

library(tidyverse)
library(viridis)

# ------------------------------------------------------------
# 1. Environmental counts from the full dataset
# ------------------------------------------------------------

env_counts <- tribble(
  ~Protein, ~Environment, ~Count,
  
  # PETase (UUT36764.1) - n = 83
  "PETase (UUT36764.1)", "Plastisphere / marine biofilm", 15,
  "PETase (UUT36764.1)", "Surface seawater", 15,
  "PETase (UUT36764.1)", "Marine host-associated", 14,
  "PETase (UUT36764.1)", "Marine sediment", 13,
  "PETase (UUT36764.1)", "Deep-sea / chemosynthetic sediment", 9,
  "PETase (UUT36764.1)", "Deep seawater", 8,
  "PETase (UUT36764.1)", "Intertidal / coastal seawater", 6,
  "PETase (UUT36764.1)", "Intertidal sediment / beach sand", 2,
  "PETase (UUT36764.1)", "Polar sea ice / cryosphere", 1,
  
  # PETase (UUT36763.1) - n = 48
  "PETase (UUT36763.1)", "Plastisphere / marine biofilm", 12,
  "PETase (UUT36763.1)", "Deep-sea / chemosynthetic sediment", 9,
  "PETase (UUT36763.1)", "Deep seawater", 8,
  "PETase (UUT36763.1)", "Marine host-associated", 8,
  "PETase (UUT36763.1)", "Surface seawater", 4,
  "PETase (UUT36763.1)", "Intertidal / coastal seawater", 2,
  "PETase (UUT36763.1)", "Intertidal sediment / beach sand", 2,
  "PETase (UUT36763.1)", "Marine sediment", 2,
  "PETase (UUT36763.1)", "Polar sea ice / cryosphere", 1,
  
  # MCL PHA depolymerase (AAQ72538.1) - n = 11
  "MCL PHA depolymerase (AAQ72538.1)",
  "Intertidal sediment / beach sand", 3,
  
  "MCL PHA depolymerase (AAQ72538.1)",
  "Marine host-associated", 3,
  
  "MCL PHA depolymerase (AAQ72538.1)",
  "Marine sediment", 2,
  
  "MCL PHA depolymerase (AAQ72538.1)",
  "Surface seawater", 2,
  
  "MCL PHA depolymerase (AAQ72538.1)",
  "Intertidal / coastal seawater", 1,
  
  # PHB depolymerase (AAB40611.1) - n = 167
  "PHB depolymerase (AAB40611.1)", "Surface seawater", 58,
  "PHB depolymerase (AAB40611.1)", "Marine host-associated", 33,
  "PHB depolymerase (AAB40611.1)", "Plastisphere / marine biofilm", 24,
  "PHB depolymerase (AAB40611.1)", "Marine sediment", 11,
  "PHB depolymerase (AAB40611.1)", "Deep seawater", 10,
  "PHB depolymerase (AAB40611.1)", "Deep-sea / chemosynthetic sediment", 10,
  "PHB depolymerase (AAB40611.1)", "Intertidal / coastal seawater", 8,
  "PHB depolymerase (AAB40611.1)", "Intertidal sediment / beach sand", 6,
  "PHB depolymerase (AAB40611.1)", "Mangrove sediment", 5,
  "PHB depolymerase (AAB40611.1)", "Polar sea ice / cryosphere", 2,
  
  # P3HB depolymerase (LC127088.1) - n = 188
  "P3HB depolymerase (LC127088.1)", "Surface seawater", 78,
  "P3HB depolymerase (LC127088.1)", "Marine host-associated", 44,
  "P3HB depolymerase (LC127088.1)", "Plastisphere / marine biofilm", 18,
  "P3HB depolymerase (LC127088.1)", "Intertidal sediment / beach sand", 13,
  "P3HB depolymerase (LC127088.1)", "Deep seawater", 11,
  "P3HB depolymerase (LC127088.1)", "Deep-sea / chemosynthetic sediment", 6,
  "P3HB depolymerase (LC127088.1)", "Mangrove sediment", 6,
  "P3HB depolymerase (LC127088.1)", "Intertidal / coastal seawater", 5,
  "P3HB depolymerase (LC127088.1)", "Marine sediment", 5,
  "P3HB depolymerase (LC127088.1)", "Polar sea ice / cryosphere", 2,
  
  # Laccase (UVG67878.1) - n = 140
  "Laccase (UVG67878.1)", "Marine host-associated", 70,
  "Laccase (UVG67878.1)", "Surface seawater", 42,
  "Laccase (UVG67878.1)", "Deep seawater", 6,
  "Laccase (UVG67878.1)", "Marine sediment", 6,
  "Laccase (UVG67878.1)", "Polar sea ice / cryosphere", 6,
  "Laccase (UVG67878.1)", "Deep-sea / chemosynthetic sediment", 4,
  "Laccase (UVG67878.1)", "Intertidal / coastal seawater", 3,
  "Laccase (UVG67878.1)", "Intertidal sediment / beach sand", 2,
  "Laccase (UVG67878.1)", "Mangrove sediment", 1,
  
  # GPX (WP_026826575.1) - n = 45
  "GPX (WP_026826575.1)", "Marine host-associated", 11,
  "GPX (WP_026826575.1)", "Marine sediment", 8,
  "GPX (WP_026826575.1)", "Intertidal sediment / beach sand", 6,
  "GPX (WP_026826575.1)", "Surface seawater", 5,
  "GPX (WP_026826575.1)", "Mangrove sediment", 4,
  "GPX (WP_026826575.1)", "Deep seawater", 3,
  "GPX (WP_026826575.1)", "Deep-sea / chemosynthetic sediment", 3,
  "GPX (WP_026826575.1)", "Intertidal / coastal seawater", 3,
  "GPX (WP_026826575.1)", "Plastisphere / marine biofilm", 2,
  
  # Oxidoreductase (UOM43036.1) - n = 14
  "Oxidoreductase (UOM43036.1)", "Marine sediment", 5,
  "Oxidoreductase (UOM43036.1)", "Marine host-associated", 3,
  "Oxidoreductase (UOM43036.1)", "Plastisphere / marine biofilm", 3,
  "Oxidoreductase (UOM43036.1)", "Deep seawater", 2,
  "Oxidoreductase (UOM43036.1)", "Mangrove sediment", 1,
  
  # GPX (WP_003945816.1) - n = 22
  "GPX (WP_003945816.1)", "Marine host-associated", 13,
  "GPX (WP_003945816.1)", "Marine sediment", 4,
  "GPX (WP_003945816.1)", "Intertidal / coastal seawater", 2,
  "GPX (WP_003945816.1)", "Intertidal sediment / beach sand", 1,
  "GPX (WP_003945816.1)", "Mangrove sediment", 1,
  "GPX (WP_003945816.1)", "Surface seawater", 1
)

# ------------------------------------------------------------
# 2. Define plotting order
# ------------------------------------------------------------

environment_order <- c(
  "Plastisphere / marine biofilm",
  "Marine host-associated",
  "Mangrove sediment",
  "Polar sea ice / cryosphere",
  "Deep-sea / chemosynthetic sediment",
  "Intertidal sediment / beach sand",
  "Intertidal / coastal seawater",
  "Deep seawater",
  "Marine sediment",
  "Surface seawater"
)

protein_order <- c(
  "PETase (UUT36764.1)",
  "PETase (UUT36763.1)",
  "MCL PHA depolymerase (AAQ72538.1)",
  "PHB depolymerase (AAB40611.1)",
  "P3HB depolymerase (LC127088.1)",
  "Laccase (UVG67878.1)",
  "GPX (WP_026826575.1)",
  "Oxidoreductase (UOM43036.1)",
  "GPX (WP_003945816.1)"
)

# ------------------------------------------------------------
# 3. Calculate percentage within each gene dataset
# ------------------------------------------------------------

heatmap_data <- env_counts %>%
  complete(
    Protein = protein_order,
    Environment = environment_order,
    fill = list(Count = 0)
  ) %>%
  group_by(Protein) %>%
  mutate(
    Total = sum(Count),
    Percent = (Count / Total) * 100
  ) %>%
  ungroup() %>%
  mutate(
    # Zero-count cells become NA so they plot as gray
    Percent_plot = ifelse(Count == 0, NA, Percent),
    
    Protein = factor(
      Protein,
      levels = rev(protein_order)
    ),
    
    Environment = factor(
      Environment,
      levels = environment_order
    )
  )

# ------------------------------------------------------------
# 4. Check totals
# ------------------------------------------------------------

heatmap_data %>%
  distinct(Protein, Total) %>%
  arrange(Protein)

# ------------------------------------------------------------
# 5. Create heatmap
# ------------------------------------------------------------

p <- ggplot(
  heatmap_data,
  aes(
    x = Environment,
    y = Protein,
    fill = Percent_plot
  )
) +
  
  geom_tile(
    color = "white",
    linewidth = 0.7
  ) +
  
  # Percentage and n within each populated cell
  geom_text(
    aes(
      label = ifelse(
        Count == 0,
        "",
        paste0(
          round(Percent, 1),
          "%\n(n=", Count, ")"
        )
      )
    ),
    size = 3,
    lineheight = 0.9
  ) +
  
  # Plasma palette used previously
  # Empty cells = light gray
  scale_fill_viridis_c(
    option = "C",
    direction = -1,
    na.value = "grey85",
    limits = c(
      0,
      max(heatmap_data$Percent, na.rm = TRUE)
    ),
    name = "% of sequences"
  ) +
  
  labs(
    x = "Isolation environment",
    y = NULL,
    title = paste(
      "Environmental distribution of",
      "candidate plastic-degradation-associated proteins"
    )
  ) +
  
  theme_classic(base_size = 12) +
  
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      size = 10,
      color = "black"
    ),
    
    axis.text.y = element_text(
      size = 10,
      color = "black"
    ),
    
    axis.ticks = element_blank(),
    
    panel.grid = element_blank(),
    
    plot.title = element_text(
      size = 14,
      face = "bold",
      hjust = 0
    ),
    
    legend.title = element_text(
      face = "bold"
    ),
    
    legend.text = element_text(
      size = 10
    )
  )

p

# ------------------------------------------------------------
# 6. Save publication-quality versions
# ------------------------------------------------------------
setwd("/Users/sabrinaelkassas/Downloads")
ggsave(
  filename = "Environmental_distribution_heatmap_full_dataset.png",
  plot = p,
  width = 11,
  height = 7,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = "Environmental_distribution_heatmap_full_dataset.svg",
  plot = p,
  width = 11,
  height = 7,
  units = "in",
  bg = "white"
)
