
## Análise das células Balf em indivíduos com COVID-19

<!-- badges: start -->

<!-- badges: end -->


Este repositório contém os scripts utilizados no fluxo de análise
realizado. A organização é a seguinte:

1.  Identificação e controle de qualidade das células sequenciadas em
    cada amostra (indivíduo).
      - `00__descr.R`
    
      - `01__get_qc_plots.R`
    
      - `02__read_apply_qc.R`
2.  Normalização das células dentro de cada paciente.
      - `03__normalization.R`
3.  Seleção das *features*.
      - `04__feature_selection.R`
4.  Investigação e correção do efeito de lote (paciente/indivíduo).
      - `05__batch_correction.R`
5.  Análise de agrupamento para identificar subpopulações
    (tipos/estados) de células.
      - `06__cluster.R`
6.  Análise de anotação dos agrupamentos para identificar o “tipo” das
    células com maior predominância em cada grupo.
      - `07__cluster_annotation.R`
7.  Identificação dos genes marcadores de cada subpopulação de células.
      - `08__cluster_markers.R`
8.  Análise DE/DS para comparação global e dentro e cada subpopulação de
    células entre as condições dos indívudos (HC, S e M).
      - `09__DE_mixed_model.R`
    
      - `09__DE_pseduo_bulk_edgeR.R`
    
      - `10__DE_comparison.R`

As matrizes de contagens dos 12 indivíduos estão disponíveis no GEO pelo
código
[GSM4339769](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4339769)
