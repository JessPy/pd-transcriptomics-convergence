# >>>>>>>>>>>>>>>>>>>>>>>
# Arquivo: main.R
# Descrição: Orquestrador do pipeline para DE de RNA-Seq + Análise de Enriquecimento
#
# Uso:
#   Rscript main.R
#   OU source("main.R") de uma sessão R
#
# Arquivos de entrada esperados:
#   data/GSE90469/counts/LRRK2_mRNA_counts_with_symbols.csv
#     Colunas: gene_id, gene_symbol, <amostra_1>, <amostra_2>, ...
#   data/GSE90469/counts/metadata.csv
#     Colunas: sample, condition (Control | Mutant)
#
# Estrutura de saída:
#   results/   ← Tabelas TSV (resultados DE, DEGs, contagens normalizadas, enriquecimento)
#   plots/     ← Figuras PNG (CQ, volcano, MA, GO, GSEA)
# <<<<<<<<<<<<<<<<<<<<<<<

# >>>>>>>>>>>>>>>>>>>>>>>
# 0. Inicialização
# <<<<<<<<<<<<<<<<<<<<<<<
cat(rep("=", 60), "\n", sep = "")
cat("  Pipeline de DE de RNA-Seq + Enriquecimento\n")
cat(format(Sys.time(), "  Iniciado em: %Y-%m-%d %H:%M:%S\n"))
cat(rep("=", 60), "\n", sep = "")

# Ajustando caminhos de source relativos à raiz do projeto
source("bin/deseq/src/utils/utils.R")
source("bin/deseq/src/r_scripts/deseq2_pipeline.R")
source("bin/deseq/src/r_scripts/visualization.R")
source("bin/deseq/src/r_scripts/gsea_pipeline.R")

# >>>>>>>>>>>>>>>>>>>>>>>
# Resolver conflitos de namespace introduzidos por pacotes do Bioconductor
# (AnnotationDbi::select e stats::filter mascaram equivalentes do dplyr)
# <<<<<<<<<<<<<<<<<<<<<<<
select <- dplyr::select
filter <- dplyr::filter
rename <- dplyr::rename

# Definindo diretórios de saída dentro da estrutura data/results
out_tables <- "data/results/tables"
out_plots  <- "data/results/plots"
setup_dirs(c(out_tables, out_plots))

# >>>>>>>>>>>>>>>>>>>>>>>
# 1. Configuração da análise
# <<<<<<<<<<<<<<<<<<<<<<<
# Para adicionar mais coortes, anexe entradas a esta lista.
# Cada entrada aciona uma execução completa de DE + enriquecimento.
cohorts <- list(
  GBA = list(
    counts_file     = "data/processed/GBA/GBA_mRNA_counts_with_symbols.csv",
    meta_file       = "data/processed/GBA/GBA_metadata.csv",
    comparisons     = list(
      "IVS_het" = c("GBA1 +/+ (Controles Saudáveis / Wild Type)", "GBA1 IVS/+ (Mutante Heterozigoto)"),
      "IVS_hom" = c("GBA1 +/+ (Controles Saudáveis / Wild Type)", "GBA1 IVS/IVS (Mutante Homozigoto Nulo)")
    ),
    combined_design = ~ condition
  ),
  LRRK2 = list(
    counts_file = "data/processed/LRRK2/LRRK2_mRNA_counts_with_symbols.csv",
    meta_file   = "data/processed/LRRK2/LRRK2_metadata.csv",
    comparisons = c("Isogenic", "Non-Isogenic")
  ),
  SNCA = list(
    counts_file = "data/processed/SNCA/SNCA_mRNA_counts_with_symbols.csv",
    meta_file   = "data/processed/SNCA/SNCA_metadata.csv",
    comparisons = NULL
  )
)

# >>>>>>>>>>>>>>>>>>>>>>>
# Limites globais — aplicados consistentemente em todos os coortes
# <<<<<<<<<<<<<<<<<<<<<<<
thresholds <- list(
  lfc_threshold  = 1,      # |log2FC| ≥ 1 para chamar um DEG
  padj_threshold = 0.05,   # FDR < 0.05
  min_count      = 10,     # pré-filtragem: contagem bruta mínima por gene
  padj_ora       = 0.05,   # corte de FDR para GO ORA
  padj_gsea      = 0.05,   # corte de FDR para GSEA
  min_gs_size    = 15,     # tamanho mínimo do conjunto de genes para GSEA
  max_gs_size    = 500     # tamanho máximo do conjunto de genes para GSEA
)

# >>>>>>>>>>>>>>>>>>>>>>>
# 2. Executar pipeline para cada coorte
# <<<<<<<<<<<<<<<<<<<<<<<

# Helper para rodar o fluxo completo para um alvo (independente se subset ou combined)
run_full_pipeline <- function(target_name, counts_file, meta_file, 
                              subset_col = NULL, subset_val = NULL, 
                              design_formula = ~ condition, thresholds) {
  # Passo 1 — Expressão Diferencial (DESeq2)
  de_results <- run_deseq2_analysis(
    target_name    = target_name,
    counts_file    = counts_file,
    meta_file      = meta_file,
    design_formula = design_formula,
    subset_col     = subset_col,
    subset_val     = subset_val,
    lfc_threshold  = thresholds$lfc_threshold,
    padj_threshold = thresholds$padj_threshold,
    min_count      = thresholds$min_count
  )

  # Passo 2 — Visualização (Volcano + MA)
  generate_volcano_plot(
    res_shrunken   = de_results$res_shrunken,
    target_name    = target_name,
    lfc_threshold  = thresholds$lfc_threshold,
    padj_threshold = thresholds$padj_threshold
  )

  generate_ma_plot(
    res_shrunken   = de_results$res_shrunken,
    target_name    = target_name,
    padj_threshold = thresholds$padj_threshold
  )

  # Passo 3 — Enriquecimento (GO ORA + GSEA)
  run_enrichment_analysis(
    res_shrunken = de_results$res_shrunken,
    res_raw      = de_results$res_raw,
    degs         = de_results$degs,
    target_name  = target_name,
    padj_ora     = thresholds$padj_ora,
    padj_gsea    = thresholds$padj_gsea,
    min_gs_size  = thresholds$min_gs_size,
    max_gs_size  = thresholds$max_gs_size
  )
}

pipeline_status <- lapply(names(cohorts), function(cohort_id) {
  cfg <- cohorts[[cohort_id]]

  cat("\n", rep("-", 60), "\n", sep = "")
  log_info("Processando coorte: ", cohort_id)
  cat(rep("-", 60), "\n", sep = "")

  tryCatch({
    # A. Análises Individuais (Subsets)
    if (!is.null(cfg$comparisons)) {
      # Suporta tanto vetores (legado) quanto listas nomeadas (para multi-genótipos)
      comp_names <- if (is.list(cfg$comparisons)) names(cfg$comparisons) else cfg$comparisons
      
      for (comp_idx in seq_along(comp_names)) {
        c_name <- comp_names[comp_idx]
        c_val  <- if (is.list(cfg$comparisons)) cfg$comparisons[[c_name]] else c_name
        
        log_info("Rodando análise para comparação: ", c_name)
        run_full_pipeline(
          target_name = paste0(cohort_id, "_", c_name),
          counts_file = cfg$counts_file,
          meta_file   = cfg$meta_file,
          subset_col  = "comparison",
          subset_val  = c_val,
          thresholds  = thresholds
        )
      }
    }

    # B. Análise Combinada (Pool)
    log_info("Rodando análise combinada para o coorte: ", cohort_id)
    
    # Se combined_design for fornecido, usa-o. Caso contrário, se houver comparações
    # que NÃO são uma lista (como no caso LRRK2), tenta controlar por comparação.
    # Para listas (GBA), cai no default ~ condition para evitar colinearidade.
    comb_design <- if (!is.null(cfg$combined_design)) {
      cfg$combined_design
    } else if (!is.null(cfg$comparisons) && !is.list(cfg$comparisons)) {
      ~ comparison + condition
    } else {
      ~ condition
    }

    run_full_pipeline(
      target_name    = paste0(cohort_id, "_Combined"),
      counts_file    = cfg$counts_file,
      meta_file      = cfg$meta_file,
      design_formula = comb_design,
      thresholds     = thresholds
    )

    log_info("Coorte concluído: ", cohort_id)
    list(cohort = cohort_id, status = "SUCCESS", error = NA)

  }, error = function(e) {
    log_error("Coorte FALHOU: ", cohort_id)
    log_error("Razão: ", conditionMessage(e))
    list(cohort = cohort_id, status = "FAILED", error = conditionMessage(e))
  })
})

# >>>>>>>>>>>>>>>>>>>>>>>
# 3. Relatório de resumo
# <<<<<<<<<<<<<<<<<<<<<<<
cat("\n", rep("=", 60), "\n", sep = "")
cat("  Resumo do Pipeline\n")
cat(rep("=", 60), "\n", sep = "")

for (s in pipeline_status) {
  status_icon <- if (s$status == "SUCCESS") "✓" else "✗"
  cat(sprintf("  %s  %-20s %s\n", status_icon, s$cohort, s$status))
  if (!is.na(s$error)) cat("       Erro:", s$error, "\n")
}

save_session_info(out_dir = "results")
cat(format(Sys.time(), "\n  Finalizado em: %Y-%m-%d %H:%M:%S\n"))
cat(rep("=", 60), "\n", sep = "")
