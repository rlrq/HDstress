# HD stress in *Brassica rapa*

Script files for network analysis:
- network_analysis.py
  - Main script for analysis of network dispersion and master regulatory nodes
- cx2_to_tsv.py
  - Converts cx2 format to tsv file
    - To be used by plot_for_manuscript.R for plotting network layouts
  - Cytoscape's network layout seems cleaner than the layouts available through ggraph, so we use it instead
  - Open *.edges in Cytoscape and exporting as cx2 format to obtain cx2 files
- plot_for_manuscript.R
  - Generate plots for figures 5, 6, S14, and S15

Script files for GO enrichment analysis of BLAST-mappable and -unmappable genes
- orthofam_analysis.sh
  - Filter orthofam data for Brassica and Arabidopsis genes only
- orthofam_analysis.R
  - Main script for assigning GO terms to Brassica genes using different methods and analysing enrichment with different universes
- plot_enrichment.R
  - Generate plots for figures S5, S6, S7

Script files for GO hierarchy analysis:
- extract_go.R
  - Subset go.obo to extract relevant GO terms and all of their parental nodes (direct and indirect)
- plot_go.R
  - Generate plots for figures 3B, 4C, S3, S11
