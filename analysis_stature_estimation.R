# ==============================================================================
# ==============================================================================
#
# THESIS: BIAS ANALYSIS IN SKELETAL STATURE ESTIMATION
# Comparative Evaluation of Osteometric Methods — Goldman Skeletal Collection
#
# Author:   [Joaquin Moreno Ibarra]
# Date:     2025
# Contact:  [jmoreno2017@udec.cl]
#
# Description:
#   This script applies eleven stature estimation methods (Pearson 1899 through
#   Menéndez et al. 2018) to long bone measurements from the Goldman dataset.
#   It performs descriptive statistics, parametric assumption testing, Bland-
#   Altman agreement analysis, a frequentist linear mixed model (LMM), and a
#   Bayesian generalized linear mixed model (BGLMM).
#
# Input:  goldman_dataset.xlsx
# Output: CSV tables, PNG figures
#
# ==============================================================================


# ==============================================================================
# 0. DEPENDENCIES AND DATA LOADING
# ==============================================================================

library(conflicted)
conflicted::conflicts_prefer(dplyr::filter)

library(dplyr)
library(tidyr)
library(readr)
library(tidyverse)
library(readxl)

# --- Load dataset and clean outliers (codes > 600 mm treated as missing) ------
datos_goldman <- readxl::read_excel("goldman_dataset.xlsx") %>%
  mutate(across(
    c(p_LMF, p_LMT, p_LMH, p_LMR),
    ~ if_else(. > 600, NA_real_, .)
  ))


# ==============================================================================
# 1. STATURE ESTIMATION FUNCTIONS (output in cm)
# ==============================================================================

# ------------------------------------------------------------------------------
# Pearson (1899)
# Reference population: European (Terry Collection)
# ------------------------------------------------------------------------------
est_pearson <- function(hueso, sexo_num, longitud_mm) {
  long_cm <- longitud_mm / 10
  if (hueso == "p_LMF") {
    dplyr::if_else(sexo_num == 0, 81.306 + 1.880 * long_cm, 72.844 + 1.945 * long_cm)
  } else if (hueso == "p_LMT") {
    dplyr::if_else(sexo_num == 0, 78.664 + 2.376 * long_cm, 74.774 + 2.352 * long_cm)
  } else if (hueso == "p_LMH") {
    dplyr::if_else(sexo_num == 0, 70.641 + 2.894 * long_cm, 71.475 + 2.754 * long_cm)
  } else if (hueso == "p_LMR") {
    dplyr::if_else(sexo_num == 0, 85.925 + 3.271 * long_cm, 81.224 + 3.343 * long_cm)
  } else { rep(NA_real_, length(longitud_mm)) }
}

# ------------------------------------------------------------------------------
# Telkka (1950)
# Reference population: Finnish (European)
# ------------------------------------------------------------------------------
est_telkka <- function(hueso, sexo_num, longitud_mm) {
  long_cm <- longitud_mm / 10
  if (hueso == "p_LMF") {
    dplyr::if_else(sexo_num == 0, 169.4 + 2.1 * (long_cm - 45.5), 156.8 + 1.8 * (long_cm - 41.8))
  } else if (hueso == "p_LMT") {
    dplyr::if_else(sexo_num == 0, 169.4 + 2.1 * (long_cm - 36.2), 156.8 + 1.9 * (long_cm - 33.1))
  } else if (hueso == "p_LMR") {
    dplyr::if_else(sexo_num == 0, 169.4 + 3.4 * (long_cm - 22.7), 156.8 + 3.1 * (long_cm - 20.8))
  } else if (hueso == "p_LMH") {
    dplyr::if_else(sexo_num == 0, 169.4 + 2.8 * (long_cm - 32.9), 156.8 + 2.7 * (long_cm - 30.7))
  } else { rep(NA_real_, length(longitud_mm)) }
}

# ------------------------------------------------------------------------------
# Trotter & Gleser (1952) — European reference population
# ------------------------------------------------------------------------------
est_trotter_eur <- function(hueso, sexo_num, longitud_mm) {
  if (hueso == "p_LMF") {
    dplyr::if_else(
      sexo_num == 0,
      ((2.38 * longitud_mm) + 614.1) / 10,
      ((2.19 * longitud_mm) + 541.0) / 10
    )
  } else { rep(NA_real_, length(longitud_mm)) }
}

# ------------------------------------------------------------------------------
# Trotter & Gleser (1958) — African-American reference population
# ------------------------------------------------------------------------------
est_trotter_afro <- function(hueso, sexo_num, longitud_mm) {
  if (hueso == "p_LMF") {
    dplyr::if_else(
      sexo_num == 0,
      ((2.11 * longitud_mm) + 703.0) / 10,
      ((2.49 * longitud_mm) + 597.6) / 10
    )
  } else { rep(NA_real_, length(longitud_mm)) }
}

# ------------------------------------------------------------------------------
# Genovés (1967)
# Reference population: Mesoamerican (Mexican indigenous)
# ------------------------------------------------------------------------------
est_genoves <- function(hueso, sexo_num, longitud_mm) {
  long_cm <- longitud_mm / 10
  if (hueso == "p_LMF") {
    dplyr::if_else(sexo_num == 0, (2.26 * long_cm) + 66.379, (2.59 * long_cm) + 49.742)
  } else if (hueso == "p_LMT") {
    dplyr::if_else(sexo_num == 0, (1.96 * long_cm) + 93.752, (2.72 * long_cm) + 63.781)
  } else { rep(NA_real_, length(longitud_mm)) }
}

# ------------------------------------------------------------------------------
# Del Angel & Cisneros (2004)
# Reference population: Contemporary Mexican
# ------------------------------------------------------------------------------
est_del_angel <- function(hueso, sexo_num, longitud_mm) {
  long_cm <- longitud_mm / 10
  if (hueso == "p_LMF") {
    dplyr::if_else(sexo_num == 0, 63.89 + 2.262 * long_cm, 47.25 + 2.588 * long_cm)
  } else if (hueso == "p_LMT") {
    dplyr::if_else(sexo_num == 0, 91.26 + 1.958 * long_cm, 94.09 + 1.919 * long_cm)
  } else if (hueso == "p_LMR") {
    dplyr::if_else(sexo_num == 0, 98.22 + 2.668 * long_cm, 66.88 + 3.926 * long_cm)
  } else if (hueso == "p_LMH") {
    dplyr::if_else(sexo_num == 0, 83.52 + 2.505 * long_cm, 32.35 + 4.160 * long_cm)
  } else { rep(NA_real_, length(longitud_mm)) }
}

# ------------------------------------------------------------------------------
# Béguelin (2011)
# Reference population: Southern South American (Patagonia)
# Sex-combined formula
# ------------------------------------------------------------------------------
est_beguelin <- function(hueso, sexo_num, longitud_mm) {
  long_cm <- longitud_mm / 10
  if      (hueso == "p_LMF") { 75.48 + (2.06 * long_cm)
  } else if (hueso == "p_LMT") { 71.60 + (2.54 * long_cm)
  } else if (hueso == "p_LMH") { 99.74 + (2.19 * long_cm)
  } else if (hueso == "p_LMR") { 103.11 + (2.61 * long_cm)
  } else { rep(NA_real_, length(longitud_mm)) }
}

# ------------------------------------------------------------------------------
# Belmonte (2011)
# Reference population: Spanish (European)
# Applies only to females (sexo_num == 1), tibia only
# ------------------------------------------------------------------------------
est_belmonte <- function(hueso, sexo_num, longitud_mm) {
  long_cm <- longitud_mm / 10
  if (hueso == "p_LMT") {
    dplyr::if_else(sexo_num == 1, 75.798 + (2.435 * long_cm), NA_real_)
  } else { rep(NA_real_, length(longitud_mm)) }
}

# ------------------------------------------------------------------------------
# Ross (2011)
# Reference population: South African (mixed ancestry)
# ------------------------------------------------------------------------------
est_ross <- function(hueso, sexo_num, longitud_mm) {
  if (hueso == "p_LMH") {
    dplyr::if_else(
      sexo_num == 0,
      (820.36 + (2.53 * longitud_mm)) / 10,
      (989.28 + (1.91 * longitud_mm)) / 10
    )
  } else if (hueso == "p_LMF") {
    dplyr::if_else(
      sexo_num == 0,
      (510.32 + (2.07 * longitud_mm)) / 10,
      (813.85 + (1.76 * longitud_mm)) / 10
    )
  } else if (hueso == "p_LMT") {
    dplyr::if_else(
      sexo_num == 0,
      (356.48 + (2.26 * longitud_mm)) / 10,
      (1026.97 + (1.41 * longitud_mm)) / 10
    )
  } else { rep(NA_real_, length(longitud_mm)) }
}

# ------------------------------------------------------------------------------
# Abarca (2013)
# Reference population: Chilean
# ------------------------------------------------------------------------------
est_abarca <- function(hueso, sexo_num, longitud_mm) {
  long_cm <- longitud_mm / 10
  if (hueso == "p_LMF") {
    dplyr::if_else(sexo_num == 0, 47.2 + 2.62 * long_cm, 74.25 + 2.10 * long_cm)
  } else { rep(NA_real_, length(longitud_mm)) }
}

# ------------------------------------------------------------------------------
# Menéndez et al. (2018)
# Reference population: South American (Argentina)
# ------------------------------------------------------------------------------
est_menendez <- function(hueso, sexo_num, longitud_mm) {
  long_cm <- longitud_mm / 10
  if (hueso == "p_LMH") {
    dplyr::if_else(sexo_num == 0, 51.6404 + (3.5756 * long_cm), 29.2104 + (4.2294 * long_cm))
  } else if (hueso == "p_LMF") {
    dplyr::if_else(sexo_num == 0, 58.5371 + (2.4211 * long_cm), 25.7316 + (3.1379 * long_cm))
  } else if (hueso == "p_LMT") {
    dplyr::if_else(sexo_num == 0, 62.1694 + (2.7730 * long_cm), 51.5941 + (3.0067 * long_cm))
  } else { rep(NA_real_, length(longitud_mm)) }
}


# ==============================================================================
# 2. APPLICATION OF FORMULAS TO THE DATASET
# ==============================================================================

analisis_estatura_completo <- datos_goldman %>%
  mutate(
    ID_Individuo = row_number(),

    # 1. Pearson (1899)
    Est_Pearson_Femur  = est_pearson("p_LMF", sexo_biológico, p_LMF),
    Est_Pearson_Tibia  = est_pearson("p_LMT", sexo_biológico, p_LMT),
    Est_Pearson_Humero = est_pearson("p_LMH", sexo_biológico, p_LMH),
    Est_Pearson_Radio  = est_pearson("p_LMR", sexo_biológico, p_LMR),

    # 2. Telkka (1950)
    Est_Telkka_Femur  = est_telkka("p_LMF", sexo_biológico, p_LMF),
    Est_Telkka_Tibia  = est_telkka("p_LMT", sexo_biológico, p_LMT),
    Est_Telkka_Humero = est_telkka("p_LMH", sexo_biológico, p_LMH),
    Est_Telkka_Radio  = est_telkka("p_LMR", sexo_biológico, p_LMR),

    # 3. Trotter & Gleser — European (1952)
    Est_TrotterEur_Femur  = est_trotter_eur("p_LMF", sexo_biológico, p_LMF),

    # 4. Trotter & Gleser — African-American (1958)
    Est_TrotterAfro_Femur = est_trotter_afro("p_LMF", sexo_biológico, p_LMF),

    # 5. Genovés (1967)
    Est_Genoves_Femur = est_genoves("p_LMF", sexo_biológico, p_LMF),
    Est_Genoves_Tibia = est_genoves("p_LMT", sexo_biológico, p_LMT),

    # 6. Del Angel & Cisneros (2004)
    Est_DelAngel_Femur  = est_del_angel("p_LMF", sexo_biológico, p_LMF),
    Est_DelAngel_Tibia  = est_del_angel("p_LMT", sexo_biológico, p_LMT),
    Est_DelAngel_Humero = est_del_angel("p_LMH", sexo_biológico, p_LMH),
    Est_DelAngel_Radio  = est_del_angel("p_LMR", sexo_biológico, p_LMR),

    # 7. Béguelin (2011)
    Est_Beguelin_Femur  = est_beguelin("p_LMF", sexo_biológico, p_LMF),
    Est_Beguelin_Tibia  = est_beguelin("p_LMT", sexo_biológico, p_LMT),
    Est_Beguelin_Humero = est_beguelin("p_LMH", sexo_biológico, p_LMH),
    Est_Beguelin_Radio  = est_beguelin("p_LMR", sexo_biológico, p_LMR),

    # 8. Belmonte (2011) — females, tibia only
    Est_Belmonte_Tibia = est_belmonte("p_LMT", sexo_biológico, p_LMT),

    # 9. Ross (2011)
    Est_Ross_Femur  = est_ross("p_LMF", sexo_biológico, p_LMF),
    Est_Ross_Tibia  = est_ross("p_LMT", sexo_biológico, p_LMT),
    Est_Ross_Humero = est_ross("p_LMH", sexo_biológico, p_LMH),

    # 10. Abarca (2013)
    Est_Abarca_Femur = est_abarca("p_LMF", sexo_biológico, p_LMF),

    # 11. Menéndez et al. (2018)
    Est_Menendez_Femur  = est_menendez("p_LMF", sexo_biológico, p_LMF),
    Est_Menendez_Tibia  = est_menendez("p_LMT", sexo_biológico, p_LMT),
    Est_Menendez_Humero = est_menendez("p_LMH", sexo_biológico, p_LMH)
  )


# ==============================================================================
# 3. DESCRIPTIVE STATISTICS
# ==============================================================================

# Pivot to long format (only estimation columns)
datos_largos <- analisis_estatura_completo %>%
  select(starts_with("Est_")) %>%
  pivot_longer(
    cols      = everything(),
    names_to  = "Metodo_Hueso",
    values_to = "Estatura_Estimada"
  ) %>%
  drop_na(Estatura_Estimada)

# Compute descriptive summary table
tabla_descriptiva <- datos_largos %>%
  group_by(Metodo_Hueso) %>%
  summarise(
    N                  = n(),
    Media              = round(mean(Estatura_Estimada), 2),
    Desviacion_Estandar = round(sd(Estatura_Estimada), 2),
    Minimo             = round(min(Estatura_Estimada), 2),
    Maximo             = round(max(Estatura_Estimada), 2)
  ) %>%
  arrange(Media)

print(tabla_descriptiva, n = 35)

write.csv2(tabla_descriptiva, "Estadistica_Descriptiva_Definitiva.csv", row.names = FALSE)


# ==============================================================================
# 4. PARAMETRIC ASSUMPTION TESTING
# ==============================================================================

# --- 4.1 Kolmogorov-Smirnov normality test ------------------------------------
prueba_ks <- ks.test(
  datos_largos$Estatura_Estimada, "pnorm",
  mean = mean(datos_largos$Estatura_Estimada),
  sd   = sd(datos_largos$Estatura_Estimada)
)

message("--- Kolmogorov-Smirnov Normality Test ---")
print(prueba_ks)

# --- 4.2 Levene's test for homogeneity of variances ---------------------------
library(car)

prueba_levene <- leveneTest(
  Estatura_Estimada ~ as.factor(Metodo_Hueso),
  data = datos_largos
)

message("--- Levene's Test for Equality of Variances ---")
print(prueba_levene)


# ==============================================================================
# 5. INFERENTIAL ANALYSIS — WELCH ANOVA AND GAMES-HOWELL POST-HOC
# ==============================================================================

# --- 5.1 Welch's one-way ANOVA (robust to unequal variances) ------------------
anova_welch <- oneway.test(
  Estatura_Estimada ~ Metodo_Hueso,
  data      = datos_largos,
  var.equal = FALSE
)

message("--- Welch's One-Way ANOVA ---")
print(anova_welch)

# --- 5.2 Games-Howell post-hoc test -------------------------------------------
library(rstatix)

posthoc_gh <- games_howell_test(datos_largos, Estatura_Estimada ~ Metodo_Hueso)

# Filter significantly different pairs (adjusted p < 0.05)
diferencias_significativas <- posthoc_gh %>%
  filter(p.adj < 0.05)

message(paste(
  "There are", nrow(diferencias_significativas),
  "method pairs that differ significantly (p.adj < 0.05)."
))

# Display the 10 most extreme differences
head(diferencias_significativas %>% arrange(p.adj), 10)


# ==============================================================================
# 6. BLAND-ALTMAN AGREEMENT ANALYSIS (ALL PAIRWISE COMBINATIONS)
# ==============================================================================

library(ggplot2)
library(stringr)

# --- 6.1 Generate all pairwise combinations -----------------------------------
columnas_est <- grep("^Est_", names(analisis_estatura_completo), value = TRUE)
pares        <- combn(columnas_est, 2)

lista_resultados <- list()

for (i in 1:ncol(pares)) {
  metodo1 <- pares[1, i]
  metodo2 <- pares[2, i]

  df_par <- analisis_estatura_completo %>%
    select(all_of(metodo1), all_of(metodo2)) %>%
    drop_na()

  # Only compute when at least 30 shared observations are available
  if (nrow(df_par) > 30) {
    m1         <- df_par[[metodo1]]
    m2         <- df_par[[metodo2]]
    diferencia <- m1 - m2
    sesgo      <- mean(diferencia)
    sd_diff    <- sd(diferencia)
    lim_inf    <- sesgo - (1.96 * sd_diff)
    lim_sup    <- sesgo + (1.96 * sd_diff)
    amplitud   <- lim_sup - lim_inf

    lista_resultados[[i]] <- data.frame(
      Metodo_A       = metodo1,
      Metodo_B       = metodo2,
      N_Comun        = nrow(df_par),
      Sesgo_Medio    = round(sesgo,    2),
      SD_Diferencia  = round(sd_diff,  2),
      Lim_Inferior_95 = round(lim_inf, 2),
      Lim_Superior_95 = round(lim_sup, 2),
      Amplitud_Error = round(amplitud, 2)
    )
  }
}

# --- 6.2 Compile and export results -------------------------------------------
tabla_bland_altman_total <- bind_rows(lista_resultados) %>%
  arrange(Amplitud_Error)

message("--- TOP 10: Most Concordant Method Pairs ---")
print(head(tabla_bland_altman_total, 10))

message("--- TOP 10: Most Discordant Method Pairs ---")
print(tail(tabla_bland_altman_total, 10))

write.csv2(tabla_bland_altman_total, "Matriz_BlandAltman_Total.csv", row.names = FALSE)

# --- 6.3 Figure 5: Forest plot — best vs worst pairs -------------------------
top_10_mejores          <- head(tabla_bland_altman_total, 10)
top_10_mejores$Categoria <- "A. Top 10: Highest Agreement"

top_10_peores           <- tail(tabla_bland_altman_total, 10)
top_10_peores$Categoria  <- "B. Top 10: Lowest Agreement"

datos_ba_extremos <- bind_rows(top_10_mejores, top_10_peores) %>%
  mutate(
    Metodo_A_Limpio = str_replace(Metodo_A, "Est_", ""),
    Metodo_B_Limpio = str_replace(Metodo_B, "Est_", ""),
    Comparacion     = paste(Metodo_A_Limpio, "vs", Metodo_B_Limpio)
  )

datos_ba_extremos$Comparacion <- factor(
  datos_ba_extremos$Comparacion,
  levels = datos_ba_extremos$Comparacion[order(datos_ba_extremos$Amplitud_Error)]
)

g5 <- ggplot(datos_ba_extremos, aes(x = Sesgo_Medio, y = Comparacion, color = Categoria)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
  geom_errorbarh(aes(xmin = Lim_Inferior_95, xmax = Lim_Superior_95),
                 height = 0.3, linewidth = 1.2) +
  geom_point(size = 3.5) +
  facet_wrap(~Categoria, scales = "free_y", ncol = 1) +
  scale_color_manual(values = c(
    "A. Top 10: Highest Agreement" = "steelblue",
    "B. Top 10: Lowest Agreement"  = "darkred"
  )) +
  theme_bw() +
  labs(
    title    = "Forensic Error Range: Best vs Worst Scenarios",
    subtitle = "95% Limits of Agreement (Bland-Altman) within-individual. Range in cm.",
    x        = "Mean Bias and Limits of Agreement (cm)",
    y        = "Comparison (Method A vs Method B)"
  ) +
  theme(
    legend.position  = "none",
    plot.title       = element_text(face = "bold", size = 15),
    plot.subtitle    = element_text(size = 12, color = "gray30"),
    axis.text.y      = element_text(size = 10),
    strip.text       = element_text(face = "bold", size = 12, color = "white"),
    strip.background = element_rect(fill = "gray20")
  )

print(g5)
ggsave("Grafico_5_BlandAltman_Extremos.png", plot = g5, width = 12, height = 9, dpi = 300)
message("Saved: Grafico_5_BlandAltman_Extremos.png")


# ==============================================================================
# 7. LINEAR MIXED MODEL (LMM — FREQUENTIST)
# ==============================================================================

library(lme4)
library(sjPlot)

# --- 7.1 Fit mixed model ------------------------------------------------------
# Fixed effects: method, bone, country of origin
# Random effect: random intercept by individual (repeated measurements)
modelo_mixto <- lmer(
  Estatura_Estimada ~ Metodo + Hueso + país + (1 | ID_Individuo),
  data = datos_modelo
)

summary(modelo_mixto)

# --- 7.2 Extract coefficient names for plotting -------------------------------
nombres_coef <- names(fixef(modelo_mixto))
print(nombres_coef)

coef_metodos <- nombres_coef[str_detect(nombres_coef, "^Metodo")]
coef_paises  <- nombres_coef[str_detect(nombres_coef, "^país")]

# --- 7.3 Figure 6a: Effect of estimation methods ------------------------------
g6a <- plot_model(
  modelo_mixto,
  type         = "est",
  terms        = coef_metodos,
  show.values  = TRUE,
  value.offset = 0.3,
  colors       = "#2166AC",
  title        = "Effect of Estimation Methods on Stature (LMM)",
  axis.title   = "Difference in cm relative to reference method",
  vline.color  = "gray30",
  dot.size     = 2.5,
  line.size    = 0.8
) +
  scale_y_discrete(labels = function(x) str_replace(x, "^Metodo", "")) +
  theme_bw(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", size = 15, hjust = 0, color = "#2166AC"),
    plot.subtitle    = element_text(size = 11, color = "gray40"),
    axis.text.y      = element_text(size = 10),
    axis.text.x      = element_text(size = 10),
    panel.grid.minor = element_blank()
  ) +
  labs(
    subtitle = "LMM coefficients with 95% confidence intervals.",
    caption  = "Dashed line = no difference relative to reference method."
  )

print(g6a)
ggsave("Grafico_6a_LMM_Metodos.png", plot = g6a, width = 11, height = 8, dpi = 300)
message("Saved: Grafico_6a_LMM_Metodos.png")

# --- 7.4 Figure 6b: Effect of country of origin -------------------------------
g6b <- plot_model(
  modelo_mixto,
  type         = "est",
  terms        = coef_paises,
  show.values  = TRUE,
  value.offset = 0.3,
  colors       = "#1A9641",
  title        = "Effect of Country of Origin on Stature Estimation (LMM)",
  axis.title   = "Difference in cm relative to reference country",
  vline.color  = "gray30",
  dot.size     = 2.5,
  line.size    = 0.8
) +
  scale_y_discrete(labels = function(x) str_replace(x, "^país", "")) +
  theme_bw(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", size = 15, hjust = 0, color = "#1A9641"),
    plot.subtitle    = element_text(size = 11, color = "gray40"),
    axis.text.y      = element_text(size = 10),
    axis.text.x      = element_text(size = 10),
    panel.grid.minor = element_blank()
  ) +
  labs(
    subtitle = "LMM coefficients with 95% confidence intervals.",
    caption  = "Dashed line = no difference relative to reference country."
  )

print(g6b)
ggsave("Grafico_6b_LMM_Paises.png", plot = g6b, width = 11, height = 8, dpi = 300)
message("Saved: Grafico_6b_LMM_Paises.png")


# ==============================================================================
# 8. BAYESIAN GENERALIZED LINEAR MIXED MODEL (BGLMM) WITH brms
# ==============================================================================

library(brms)
library(tidybayes)
library(scales)

# --- 8.1 Fit Bayesian model ---------------------------------------------------
# Note: This step may take several minutes depending on hardware.
# Four MCMC chains, 2000 iterations each (1000 warmup).
message("Fitting Bayesian model. This may take several minutes...")

modelo_bayesiano <- brm(
  formula = Estatura_Estimada ~ Metodo + Hueso + país + (1 | ID_Individuo),
  data    = datos_modelo,
  family  = gaussian(),
  chains  = 4,
  cores   = 4,       # Adjust to number of available CPU cores
  iter    = 2000,
  warmup  = 1000,
  seed    = 12345
)

message("--- Bayesian Model Summary (BGLMM) ---")
summary(modelo_bayesiano)

# --- 8.2 Extract fixed effects posterior summary ------------------------------
posterior_summary <- as.data.frame(fixef(modelo_bayesiano)) %>%
  tibble::rownames_to_column("Termino") %>%
  rename(
    Estimacion = Estimate,
    Error      = Est.Error,
    Lim_Inf    = `Q2.5`,
    Lim_Sup    = `Q97.5`
  ) %>%
  mutate(
    Tipo = case_when(
      str_detect(Termino, "^Metodo") ~ "Metodo",
      str_detect(Termino, "^Hueso")  ~ "Hueso",
      str_detect(Termino, "^país")   ~ "País",
      TRUE                           ~ "Intercepto"
    ),
    Etiqueta = Termino %>%
      str_replace("^Metodo", "") %>%
      str_replace("^Hueso",  "") %>%
      str_replace("^país",   "") %>%
      str_replace("Intercept", "Intercept")
  ) %>%
  filter(Tipo != "Intercepto")

# --- 8.3 Auxiliary function: Bayesian forest plot ----------------------------
construir_g7 <- function(datos, color_principal, titulo, titulo_eje) {

  datos <- datos %>%
    mutate(Etiqueta = reorder(Etiqueta, Estimacion))

  ggplot(datos, aes(x = Estimacion, y = Etiqueta)) +

    # Practical equivalence band: ±2 cm
    annotate("rect",
             xmin = -2, xmax = 2,
             ymin = -Inf, ymax = Inf,
             fill = "gray90", alpha = 0.5) +

    geom_vline(xintercept = 0, linetype = "dashed",
               color = "gray40", linewidth = 0.8) +

    # 95% credibility interval
    geom_errorbarh(
      aes(xmin = Lim_Inf, xmax = Lim_Sup),
      height    = 0.35,
      linewidth = 1.1,
      color     = color_principal,
      alpha     = 0.75
    ) +

    # Posterior mean
    geom_point(size = 4, color = color_principal) +

    # Numeric label
    geom_text(
      aes(label = sprintf("%.1f", Estimacion)),
      hjust = -0.35,
      size  = 3.2,
      color = "gray20"
    ) +

    scale_x_continuous(
      name   = titulo_eje,
      breaks = pretty_breaks(n = 6),
      expand = expansion(mult = c(0.05, 0.12))
    ) +

    labs(
      title    = titulo,
      subtitle = "Posterior mean with 95% credibility interval.\nGrey band = practical equivalence zone (±2 cm).",
      y        = NULL,
      caption  = "Dashed line = no difference relative to reference category."
    ) +

    theme_bw(base_size = 12) +
    theme(
      plot.title         = element_text(face = "bold", size = 15,
                                        hjust = 0, color = color_principal),
      plot.subtitle      = element_text(size = 11, color = "gray40", hjust = 0),
      plot.caption       = element_text(size = 9,  color = "gray55", hjust = 0),
      axis.text.y        = element_text(size = 10, color = "gray10"),
      axis.text.x        = element_text(size = 10),
      panel.grid.minor   = element_blank(),
      panel.grid.major.y = element_blank(),
      plot.margin        = margin(15, 20, 10, 15)
    )
}

# --- 8.4 Figure 7a: Bayesian effects — estimation methods --------------------
g7a <- construir_g7(
  datos           = posterior_summary %>% filter(Tipo == "Metodo"),
  color_principal = "#2166AC",
  titulo          = "Bayesian Effects — Estimation Methods",
  titulo_eje      = "Posterior difference in cm relative to reference method"
)
print(g7a)
ggsave("Grafico_7a_Bayes_Metodos.png", plot = g7a, width = 11, height = 8, dpi = 300)
message("Saved: Grafico_7a_Bayes_Metodos.png")

# --- 8.5 Figure 7b: Bayesian effects — country of origin ---------------------
g7b <- construir_g7(
  datos           = posterior_summary %>% filter(Tipo == "País"),
  color_principal = "#1A9641",
  titulo          = "Bayesian Effects — Country of Origin",
  titulo_eje      = "Posterior difference in cm relative to reference country"
)
print(g7b)
ggsave("Grafico_7b_Bayes_Paises.png", plot = g7b, width = 11, height = 8, dpi = 300)
message("Saved: Grafico_7b_Bayes_Paises.png")

# --- 8.6 Figure 7c: Bayesian effects — skeletal element ----------------------
g7c <- construir_g7(
  datos           = posterior_summary %>% filter(Tipo == "Hueso"),
  color_principal = "#D73027",
  titulo          = "Bayesian Effects — Skeletal Element",
  titulo_eje      = "Posterior difference in cm relative to reference bone"
)
print(g7c)
ggsave("Grafico_7c_Bayes_Huesos.png", plot = g7c, width = 11, height = 8, dpi = 300)
message("Saved: Grafico_7c_Bayes_Huesos.png")

# ==============================================================================
# END OF SCRIPT
# ==============================================================================
