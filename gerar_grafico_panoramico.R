# ==============================================================================
# Script R Oficial: Gráfico 1 (Event Plot Panorâmico de Largura Total) - 70ª RBras
# ==============================================================================
# Características:
# 1. Preenche 100% da largura útil da coluna da direita (1276 bp)
# 2. Altura vertical contida (proporção panorâmica ~ 3.5 : 1) para dar espaço à Seção 4
# 3. Sem títulos internos de eixos nem legenda flutuante (conforme manual da RBras)
# 4. Paleta oficial da RBras: Azul #073A70, Cinza ardósia #52606D, Grade #E9EAEB
# ==============================================================================

suppressPackageStartupMessages({
  library(survival)
  library(ggplot2)
  library(dplyr)
})

# 1. Dados CGD
data(cgd, package = "survival")

# Tratamento dos dados e ordenação por tempo de acompanhamento (v_i)
cgd_clean <- cgd %>%
  mutate(treat_label = factor(treat, levels = c("placebo", "rIFN-g"),
                              labels = c("Placebo", "rIFN-γ")))

resumo_indiv <- cgd_clean %>%
  group_by(id) %>%
  summarise(
    n_eventos      = sum(status),
    tempo_followup = max(tstop),
    treat          = first(treat_label),
    .groups        = "drop"
  ) %>%
  arrange(tempo_followup) %>%
  mutate(ordem_grafico = row_number())

cgd_eventos <- cgd_clean %>%
  filter(status == 1) %>%
  left_join(resumo_indiv %>% select(id, ordem_grafico), by = "id")

# 2. Construção do gráfico panorâmico
fig1_wide <- ggplot() +
  # Linhas individuais de seguimento dos 128 pacientes
  geom_segment(data = resumo_indiv,
               aes(x = 0, xend = tempo_followup, y = ordem_grafico, yend = ordem_grafico),
               color = "#C7CACC", linewidth = 0.45, alpha = 0.85) +
  # Infecções recorrentes: círculos azuis institucionais (#073A70)
  geom_point(data = cgd_eventos,
             aes(x = tstop, y = ordem_grafico),
             shape = 16, color = "#073A70", size = 2.0, alpha = 0.9) +
  # Fim do seguimento / Censura (v_i): triângulos cinza ardósia (#52606D)
  geom_point(data = resumo_indiv,
             aes(x = tempo_followup, y = ordem_grafico),
             shape = 17, color = "#52606D", size = 1.7) +
  # Escalas
  scale_x_continuous(breaks = seq(0, 400, by = 50), limits = c(0, 415), expand = c(0.008, 0)) +
  scale_y_continuous(breaks = seq(0, 128, by = 32), limits = c(0, 132), expand = c(0.01, 0)) +
  # Rótulos diretos discretos na área livre inferior direita (sem caixa de legenda)
  annotate("point", x = 280, y = 20, shape = 16, color = "#073A70", size = 2.8) +
  annotate("text",  x = 295, y = 20, label = "Infecção Recorrente", hjust = 0,
           color = "#073A70", fontface = "bold", size = 4.2) +
  annotate("point", x = 280, y = 9,  shape = 17, color = "#52606D", size = 2.8) +
  annotate("text",  x = 295, y = 9,  label = "Censura (v[i])", hjust = 0,
           parse = TRUE, color = "#52606D", fontface = "bold", size = 4.2) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title         = element_blank(),
    plot.subtitle      = element_blank(),
    axis.title         = element_blank(),
    axis.text          = element_text(color = "#323F4B", size = 11, face = "bold"),
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(color = "#E9EAEB", linewidth = 0.45),
    panel.border       = element_rect(color = "#C7CACC", fill = NA, linewidth = 0.7),
    legend.position    = "none",
    plot.margin        = margin(t = 4, r = 8, b = 4, l = 6)
  )

# Salvar nas dimensões panorâmicas: 36 cm x 10 cm (proporção 3.6:1) a 300 DPI
ggsave("fig1_event_history_panoramic.png", fig1_wide,
       width = 36.0, height = 10.0, units = "cm", dpi = 300)
ggsave("fig1_event_history_panoramic.pdf", fig1_wide,
       width = 36.0, height = 10.0, units = "cm")
cat("Gráfico panorâmico gerado com sucesso!\n")
