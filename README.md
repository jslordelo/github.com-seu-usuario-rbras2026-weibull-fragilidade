# Modelo de Regressao Weibull Modificado com Fragilidade Gama Compartilhada para Processos de Contagem em Eventos Recorrentes

Repositorio complementar oficial do trabalho apresentado na 70a Reuniao Anual da Regiao Brasileira da Sociedade Internacional de Biometria (RBras), Salvador/BA (2026).

Este repositorio disponibiliza o artigo completo em PDF, o codigo-fonte em LaTeX, o script R autossuficiente e reprodutivel, bem como as figuras em alta resolucao (300 DPI) geradas para a analise de eventos recorrentes da Doenca Granulomatosa Cronica (CGD).

---

## Autores e Filiacao

* **Jairo Santo Lordelo** (Autor principal / Apresentador)  
  *Programa de Pos-Graduacao em Estatistica e Ciencia de Dados (PGECD) -- IME/UFBA*
* **Prof. Dra. Giovana Oliveira Silva** (Orientadora)  
  *Departamento de Estatistica -- Instituto de Matematica e Estatistica (IME/UFBA)*

---

## Resumo do Projeto

Na modelagem de eventos recorrentes em saude, os modelos tradicionais baseados na distribuicao Weibull classica impoem taxas de falha estritamente monotonas e assumem independencia condicional entre as observacoes. 

Para superar essas restricoes, este projeto propoe uma extensao do modelo de regressao para processos de contagem sob Processo de Poisson Nao-Homogeneo (PPNH) fundamentado na distribuicao Weibull Modificada (MW) (Lai et al., 2003; Conceicao et al., 2025), incorporando um efeito aleatorio multiplicativo de fragilidade Gama compartilhada.

### Principais Destaques Metodologicos e Clinicos
1. **Flexibilidade Basal em Banheira (Bathtub):** A distribuicao MW com parametros basais de forma gamma = 0,663 < 1 e aceleracao tau = 0,0035 > 0 (p < 0,01) captura formalmente taxas de falha com decrescimo inicial e elevacao tardia, com ponto de minimo estimado em 43,2 dias.
2. **Marginalizacao Analitica Exata:** A integracao do efeito aleatorio latente e obtida analiticamente em forma fechada via Transformada de Laplace, garantindo convergencia e estabilidade numerica ao estimador quasi-Newton (BFGS).
3. **Casos Particulares Aninhados:** O modelo proposto engloba como casos particulares o modelo MW Base (sem fragilidade), a Weibull com fragilidade Gama (taxa monotona) e a Weibull classica independente.
4. **Correcao de Sobredispersao Amostral:** O diagnostico empirico de Poisson nos dados da CGD revelou uma razao de dispersao amostral Var[N_i(t)]/E[N_i(t)] = 2,027 > 1. A estimativa da variancia da fragilidade (phi = 0,504; p = 0,027 via TRV na fronteira parametrica) absorve essa sobredispersao e ajusta os erros-padrao dos coeficientes.
5. **Eficacia Terapeutica:** O tratamento com interferon gama (rIFN-gamma) demonstrou forte efeito protetor, reduzindo em 63,4% a taxa instantanea de infeccoes graves recorrentes frente ao placebo (HR = 0,366; p = 0,0008).

---

## Estrutura do Repositorio

```text
├── README.md                           # Visao geral e instrucoes de execucao
├── analise_cgd_weibull_modificada.R    # Script R autossuficiente e reprodutivel
├── artigo_weibull_modificada_cgd_v9.pdf# Artigo completo diagramado (14 paginas)
├── artigo_weibull_modificada_cgd_v9.tex# Codigo-fonte LaTeX do artigo
├── figures/                            # Figuras em alta resolucao (300 DPI)
│   ├── fig1_event_history_wide.png     # Grafico 1: Event plot panoramico (0-400 dias, n=128)
│   ├── fig2_mcf.png                    # Grafico 2: Funcao de Media Acumulada (Nelson-Aalen)
│   └── fig3_poisson_diagnosis.png      # Grafico 3: Diagnostico empirico de sobredispersao
└── poster/                             # Arquivos do poster oficial (70a RBras)
    └── poster_rbras_completo.tex       # Codigo LaTeX com rbras-poster.cls
```

---

## Sintese dos Resultados Empiricos (Dados CGD)

| Parametro / Criterio | Weibull Base | Weibull + Gama | MW Base (Conceicao et al.) | MW + Gama (Proposto) | Razao de Risco (HR) | p-valor |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| Intercepto (beta_0) | -6,9266 | -6,7837 | -4,5624 | **-4,6627** (0,9274) | -- | <0,001 |
| Idade no ingresso (beta_2) | -0,0371 | -0,0367 | -0,0377 | **-0,0369** (0,0163) | 0,964 | 0,0236 |
| Tratamento: rIFN-gamma (beta_3) | -1,0556 | -0,9892 | -1,0701 | **-1,0032** (0,2996) | **0,366** | **0,0008** |
| Corticosteroides: Sim (beta_4) | 1,3711 | 1,5530 | 1,4005 | **1,5550** (0,7603) | 4,735 | 0,0408 |
| Heranca: Autossomica (beta_5) | 0,6471 | 0,5942 | 0,6412 | **0,5785** (0,3298) | 1,783 | 0,0793 |
| Forma basal (gamma) | 1,2544 | 1,2232 | 0,6284 | **0,6633** (0,1995) | -- | 0,0009 |
| Aceleracao basal (tau) | -- | -- | 0,0039 | **0,0035** (0,0012) | -- | 0,0035 |
| Variancia da Fragilidade (phi) | -- | 0,6064 | -- | **0,5044** (0,3230) | -- | 0,0588 |
| **Log-Verossimilhanca (MLL)** | -528,97 | -525,44 | -523,88 | **-521,42** | TRV: p = 0,027 | |
| **Criterio de Akaike (AIC)** | 1071,94 | 1066,89 | 1063,75 | **1060,84** | Delta AIC = -2,91 | |
| **Criterio Bayesiano (BIC)** | 1091,91 | 1089,70 | 1086,59 | **1086,51** | Delta BIC = -0,08 | |

---

## Instrucoes de Execucao

### Pre-requisitos
O script requer o **R (versao >= 4.0.0)** e os seguintes pacotes do CRAN:
```r
install.packages(c("survival", "ggplot2", "dplyr", "scales"))
```

### Execucao do Script
1. Clone ou baixe os arquivos deste repositorio.
2. Abra o arquivo `analise_cgd_weibull_modificada.R` no RStudio ou terminal.
3. Execute o script:
```r
source("analise_cgd_weibull_modificada.R")
```

### O que o script realiza automaticamente:
* Carrega o conjunto de dados `cgd` nativo do pacote `survival`.
* Calcula o diagnostico de sobredispersao (Var/E = 2,027).
* Gera e exporta em PDF e PNG (300 DPI) as tres figuras do estudo.
* Ajusta os quatro modelos parametricos via otimizacao numerica BFGS (`optim`).
* Inverte a matriz Hessiana observada para obter os erros-padrao assintoticos.
* Executa o Teste da Razao de Verossimilhanca (TRV) na fronteira parametrica (mistura 1/2 qui-quadrado).
* Imprime no console a tabela comparativa formatada com estimativas, erros-padrao, HRs e criterios AIC/BIC.

---

## Citacao

Se utilizar este modelo, codigos ou materiais, por favor cite:

```bibtex
@inproceedings{lordelo2026rbras,
  author    = {Lordelo, Jairo Santo and Silva, Giovana Oliveira},
  title     = {Modelo de Regressao Weibull Modificado com Fragilidade Gama para Processos de Contagem em Eventos Recorrentes},
  booktitle = {Anais da 70a Reuniao Anual da Regiao Brasileira da Sociedade Internacional de Biometria (RBras)},
  year      = {2026},
  address   = {Salvador, BA, Brasil}
}
```

---

## Licenca

Este projeto e distribuido sob a licenca **MIT**.
