library(tidyverse)
library(ontologyIndex)

extrafont::font_import(pattern = "Arial", prompt = FALSE)
theme_set(theme_bw(base_family = "Arial"))
theme_set(theme_void(base_family = "Arial"))

mkpath <- function(...){paste(..., sep = '/')}

machine_name <- Sys.info()[["nodename"]]
if (machine_name == "chaelab-rlrq"){
    dir_root <- "/home/rachelle/OneDrive/NUS/CEY"
} else if (machine_name == "rlrq-home") {
    dir_root <- "/mnt/d/OneDrive_doysd/OneDrive - Default Directory/NUS/CEY"
}

dir_proj <- mkpath(dir_root, "/zzOtherzz/XiaoMei/HDstress")


file_ext <- function(fname){
  ext <- str_extract(fname, "(?<=/.)[^.]+$")
  if (is.na(ext)){return('')}
  return(ext)
}
save_plot <- function(fname, p, w = 9.91, h = 5.44, units = "in", fmt = NA, 
                      dir = mkpath(dir_proj, "plots", "enrichment"), ...){
  dir.create(dir, showWarnings = FALSE)
  if(is.na(fmt)){
    if (file_ext(fname) != ''){
      print(h)
      ggsave(paste0(dir, "/", fname), 
             p, w = w, h = h, units = units) ##9.91x5.44in if legend.pos at side
      return()
    } else {
      fmt = c("pdf", "png")
    }
  }
  for (ext in fmt){
    ggsave(paste0(dir, "/", fname, ".", ext), 
           p, w = w, h = h, units = units) ##9.91x5.44in if legend.pos at side
  }
}


## read files
f_go <- mkpath(dir_root, "/data/GO/go.obo")
go <- get_OBO(f_go)

f_enrichment <- mkpath(dir_proj, "results/enrichment", "enrichment.txt")
df.enrichment <- read.table(f_enrichment, header = TRUE, sep = '\t') %>%
    dplyr::mutate(group = ifelse(group == "ECR", "E_CR", group)) %>%
    dplyr::mutate(group = factor(group, levels = c("E_CR", "C", "HD", "CHD", "noblast_merged", "merged")),
                  universe = ifelse(universe == "orthofam_all", "orthofam_Chiifu3.5",
                                    ifelse(universe == "blast_all", "blast_Chiifu3.5", universe))) %>%
    dplyr::filter(! GO.id %in% c("GO:0009987")) ## drop 'cellular process'

f_gosumm_blast <- mkpath(dir_proj, "results/enrichment", "go_summary.brassica.blast.tsv")
f_gosumm_of <- mkpath(dir_proj, "results/enrichment", "go_summary.brassica.orthofam.tsv")
df.gosumm.blast <- read.table(f_gosumm_blast, header = TRUE, sep = '\t')
df.gosumm.of <- read.table(f_gosumm_of, header = TRUE, sep = '\t')

f_of_counts_species <- mkpath(dir_proj, "results/enrichment", "orthofam.counts.ath-bra.tsv")
f_of_counts_bra <- mkpath(dir_proj, "results/enrichment", "orthofam.counts.brassica.tsv")
df.of.counts.species <- read.table(f_of_counts_species, header = TRUE, sep = '\t')
df.of.counts.bra <- read.table(f_of_counts_bra, header = TRUE, sep = '\t')

f_of_ath <- mkpath(dir_proj, "data/orthofam", "orthofam.ath.db3.16.0.tsv")
df.of.ath <- read.table(f_of_ath, header = TRUE, sep = '\t')

## some qc related stuff

## for GO terms that are in both blast & orthofam annotations, plot N.bra blast vs. N.bra.norm & N.bra & uniq N.ath orthofam
df.of.bra.uniqath <- df.of.counts.bra %>%
    dplyr::filter(annotated.orthofam) %>%
    dplyr::select(orthofam, group) %>%
    dplyr::distinct() %>%
    dplyr::left_join(df.of.ath, by = "orthofam", relationship = "many-to-many") %>%
    dplyr::select(orthofam, group, GOALL, N.ath.go, N.ath, fraction.ath.go) %>%
    dplyr::distinct()
df.toplot <- df.gosumm.of %>%
    dplyr::group_by(group, GOALL) %>%
    dplyr::summarise(N.uniqGO.multByBra = sum(N.bra),
                     N.normByAth.multByBra = sum(N.bra.norm)) %>% dplyr::ungroup() %>%
    dplyr::left_join(df.of.bra.uniqath %>%
                     dplyr::group_by(group, GOALL) %>% dplyr::summarise(N.multByAth = sum(N.ath.go)) %>%
                     dplyr::ungroup(), by = c("group", "GOALL")) %>%
    tidyr::gather("N.of.method", "N.of", c(N.uniqGO.multByBra, N.normByAth.multByBra, N.multByAth)) %>%
    dplyr::distinct() %>%
    dplyr::left_join(df.gosumm.blast %>%
                     dplyr::group_by(group, GOALL) %>% dplyr::summarise(N.blast = sum(N.bra)) %>% dplyr::ungroup(),
                     by = c("group", "GOALL")) %>%
    dplyr::mutate(group = ifelse(group == "ECR", "CR",
                          ifelse(group == "C", "G_control",
                          ifelse(group == "HD", "G_HD",
                          ifelse(group == "CHD", "G_cHD", group))))) %>%
    dplyr::mutate(group = factor(group, levels = c("CR", "G_control", "G_HD", "G_cHD")))
p.go.count.overview <- ggplot(df.toplot %>% dplyr::filter(!is.na(group))) +
    geom_point(aes(x = N.blast, y = N.of, colour = N.of.method), pch = 16, alpha = 0.3) +
    ## geom_smooth(aes(x = N.blast, y = N.of, colour = N.of.method), method = "auto") +
    guides(colour = guide_legend(override.aes = list(alpha = 1), title = "Normalisation method")) +
    facet_wrap(~group, scales = "free") +
    theme_bw() +
    xlab("frequency of GO term in blast-able brassica genes") +
    ylab("frequency of GO term in brassica genes mapped via orthofam")
save_plot("go_count.overview",
          p.go.count.overview +
          ggtitle("GO counting methods via orthofam against blast",
                  subtitle = paste0("# brassica genes annotated by blast: 26,690;\n",
                                    "# brassica genes annotated via orthofam: 31,028\n",
                                    "# brassica genes annotated only by blast: 2,416\n",
                                    "# brassica genes annotated only by orthofam: 6,754")) +
          theme(legend.position = "bottom"),
          h = 7, w = 7, fmt = c("png", "pdf"))
save_plot("FigureSX_Fig7supp_go_count_overview",
          p.go.count.overview +
          ggtitle("GO counting methods via orthofam against blast",
                  subtitle = paste0("# brassica genes annotated by blast: 26,690;\n",
                                    "# brassica genes annotated via orthofam: 31,028\n",
                                    "# brassica genes annotated only by blast: 2,416\n",
                                    "# brassica genes annotated only by orthofam: 6,754")) +
          theme(legend.position = "bottom"),
          h = 7, w = 7, fmt = c("png", "pdf"),
          dir = mkpath(dir_proj, "plots/figures"))

## do some parsing for enrichment analysis
df.enrichment.groupings <- df.enrichment %>%
    dplyr::group_by(universe, has.blast, treatmentDEG, group) %>%
    dplyr::summarise(N.GO = n()) %>%
    dplyr::ungroup()

## compare enrichment by group of blastable genes vs. all orthofam (as quality control for using ath in orthofam to get annotations)
max.pvalue <- 0.05
df.enrichment.best10.pvalue <- df.enrichment %>%
    dplyr::filter(pvalue < max.pvalue) %>%
    dplyr::group_by(universe, has.blast, treatmentDEG, group) %>%
    dplyr::slice_min(pvalue, n = 10, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(GO.name = sapply(GO.id, function(goid){go$name[[goid]]}))
df.enrichment.best10.Nnorm <- df.enrichment %>%
    dplyr::filter(pvalue < max.pvalue) %>%
    dplyr::group_by(universe, has.blast, treatmentDEG, group) %>%
    dplyr::slice_max(N.norm, n = 10, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(GO.name = sapply(GO.id, function(goid){go$name[[goid]]}))

df.toplot <- df.enrichment.best10.pvalue %>%
    dplyr::filter(universe %in% c("blast_all", "orthofam_Chiifu3.5", "orthofam_Z1") & treatmentDEG)
p.enrich.pvalue.control <- ggplot(df.toplot %>%
                                  dplyr::mutate(GeneRatio = sapply(GeneRatio, function(x){eval(parse(text = x))}),
                                                NnormRatio = sapply(1:nrow(.), function(i){
                                                    .[[i,"N.norm"]]/as.numeric(str_split(.[i,"GeneRatio"], '/')[[1]][[2]])
                                                }))) +
    geom_point(aes(x = NnormRatio, y = tidytext::reorder_within(GO.name, N.norm, list(universe, group)),
                   colour = p.adjust, size = N.norm), pch = 16) +
    scale_colour_gradient(high = "blue", low = "red") +
    tidytext::scale_y_reordered(labels = function(x){str_wrap(tidytext::reorder_func(x), width = 40)}) +
    ggh4x::facet_nested(universe~group, scales = "free", independent = "all") +
    theme_bw() + ylab("GO term")
df.toplot <- df.enrichment.best10.Nnorm %>%
    dplyr::filter(universe %in% c("blast_all", "orthofam_Chiifu3.5", "orthofam_Z1") & treatmentDEG)
p.enrich.N.control <- ggplot(df.toplot %>%
                             dplyr::mutate(GeneRatio = sapply(GeneRatio, function(x){eval(parse(text = x))}),
                                           NnormRatio = sapply(1:nrow(.), function(i){
                                               .[[i,"N.norm"]]/as.numeric(str_split(.[i,"GeneRatio"], '/')[[1]][[2]])[[1]]
                                           }))) +
    geom_point(aes(x = NnormRatio, y = tidytext::reorder_within(GO.name, N.norm, list(universe, group)),
                   colour = p.adjust, size = N.norm), pch = 16) +
    scale_colour_gradient(high = "blue", low = "red") +
    tidytext::scale_y_reordered(labels = function(x){str_wrap(tidytext::reorder_func(x), width = 40)}) +
    ggh4x::facet_nested(universe~group, scales = "free", independent = "all") +
    theme_bw() + ylab("GO term")

save_plot("enrich.pvalue.control",
          p.enrich.pvalue.control +
          ggtitle("treatment DEG against X universe",
                  subtitle = paste0("GO:0009987 (cellular process) dropped; p < ", max.pvalue)),
          w = 20, h = 10, fmt = c("png", "pdf"))
save_plot("enrich.Nnorm.control",
          p.enrich.N.control +
          ggtitle("treatment DEG against X universe",
                  subtitle = paste0("GO:0009987 (cellular process) dropped; p < ", max.pvalue)),
          w = 20, h = 10, fmt = c("png", "pdf"))

axis.lab.width <- 50
df.toplot <- df.enrichment.best10.pvalue %>%
    dplyr::filter(universe %in% c("blast_Chiifu3.5", "orthofam_Chiifu3.5", "orthofam_Z1")) %>%
    dplyr::mutate(treatmentDEG = ifelse(is.na(treatmentDEG), "all", ifelse(treatmentDEG, "treatmentDEG", "non-treatmentDEG"))) %>%
    dplyr::mutate(group = as.character(group)) %>%
    dplyr::mutate(group = sapply(group, function(x){ifelse(is.na(x), "intermediate",
                                                    ifelse(x=="E_CR", "X_CR",
                                                    ifelse(x=="C", "X_control",
                                                    ifelse(x=="HD", "X_HD",
                                                    ifelse(x=="CHD", "X_cHD", x)))))})) %>%
    dplyr::mutate(group = factor(group, levels = c("noblast_merged", "X_CR", "X_control", "X_HD", "X_cHD", "intermediate")),
                  mapped.by = sapply(universe, function(x){
                      ifelse(str_detect(x, "orthofam"), "mapped by Orthofam", "mapped by BLAST")}))
p.enrich.pvalue.full <- ggplot(
    df.toplot %>%
    dplyr::mutate(GeneRatio = sapply(GeneRatio, function(x){eval(parse(text = x))}),
                  NnormRatio = sapply(1:nrow(.), function(i){
                      .[[i,"N.norm"]]/as.numeric(str_split(.[i,"GeneRatio"], '/')[[1]][[2]])
                  }))) +
    geom_point(aes(x = NnormRatio, y = tidytext::reorder_within(GO.name, N.norm,
                                                                list(universe, group, treatmentDEG)),
                   colour = p.adjust, size = N.norm), pch = 16) +
    scale_colour_gradient(high = "blue", low = "red") +
    tidytext::scale_y_reordered(labels = function(x){
        str_wrap(tidytext::reorder_func(x), width = axis.lab.width)}) +
    ## facet_wrap(~universe*has.blast*group, scales = "free", nrow = 3) +
    ggh4x::facet_nested(group~treatmentDEG*mapped.by*universe, independent = "all", scales = "free") +
    theme_bw() + ylab("GO term") +
    theme(axis.text.x = element_text(angle = 90))
df.toplot <- df.enrichment.best10.Nnorm %>%
    dplyr::filter(universe %in% c("blast_Chiifu3.5", "orthofam_Chiifu3.5", "orthofam_Z1")) %>%
    dplyr::mutate(treatmentDEG = ifelse(is.na(treatmentDEG), "all",
                                 ifelse(treatmentDEG, "treatmentDEG", "non-treatmentDEG"))) %>%
    dplyr::mutate(group = as.character(group)) %>%
    dplyr::mutate(group = sapply(group, function(x){ifelse(is.na(x), "intermediate",
                                                    ifelse(x=="E_CR", "X_CR",
                                                    ifelse(x=="C", "X_control",
                                                    ifelse(x=="HD", "X_HD",
                                                    ifelse(x=="CHD", "X_cHD", x)))))})) %>%
    dplyr::mutate(group = factor(group, levels = c("noblast_merged", "X_CR", "X_control",
                                                   "X_HD", "X_cHD", "intermediate")),
                  mapped.by = sapply(universe, function(x){
                      ifelse(str_detect(x, "orthofam"), "mapped by Orthofam", "mapped by BLAST")}))
p.enrich.N.full <- ggplot(
    df.toplot %>%
    dplyr::mutate(GeneRatio = sapply(GeneRatio, function(x){eval(parse(text = x))}),
                  NnormRatio = sapply(1:nrow(.), function(i){
                      .[[i,"N.norm"]]/as.numeric(str_split(.[i,"GeneRatio"], '/')[[1]][[2]])
                  }))) +
    geom_point(aes(x = NnormRatio,
                   y = tidytext::reorder_within(GO.name, N.norm, list(universe, group, treatmentDEG)),
                   colour = p.adjust, size = N.norm), pch = 16) +
    scale_colour_gradient(high = "blue", low = "red") +
    tidytext::scale_y_reordered(labels = function(x){
        str_wrap(tidytext::reorder_func(x), width = axis.lab.width)}) +
    ## facet_wrap(~universe*has.blast*group, scales = "free", nrow = 3) +
    ggh4x::facet_nested(group~treatmentDEG*mapped.by*universe, independent = "all", scales = "free") +
    theme_bw() + ylab("GO term") +
    theme(axis.text.x = element_text(angle = 90))

save_plot("enrich.pvalue.full",
          p.enrich.pvalue.full +
          theme(axis.title.y = element_blank()) +
          ggtitle("GO enrichment with different mapping methods & universes",
                  subtitle = paste0("facets top to bottom: treatment DEG, map method, universe\n",
                                    "GO:0009987 (cellular process) dropped; p < ", max.pvalue)),
          w = 25, h = 15, fmt = c("png", "pdf"))
save_plot("enrich.Nnorm.full",
          p.enrich.N.full +
          theme(axis.title.y = element_blank()) +
          ggtitle("GO enrichment with different mapping methods & universes",
                  subtitle = paste0("facets top to bottom: treatment DEG, map method, universe\n",
                                    "GO:0009987 (cellular process) dropped; p < ", max.pvalue)),
          w = 25, h = 15, fmt = c("png", "pdf"))

save_plot("FigureSX_Fig7supp_map_method_overview_best10pvalue",
          p.enrich.pvalue.full +
          theme(axis.title.y = element_blank()) +
          ggtitle("GO enrichment with different mapping methods & universes",
                  subtitle = paste0("facets top to bottom: treatment DEG, map method, universe\n",
                                    "GO:0009987 (cellular process) dropped; p < ", max.pvalue)),
          w = 25, h = 15, fmt = c("png", "pdf"),
          dir = mkpath(dir_proj, "plots/figures"))

## separate 'all' & 'treatmentDEG'
p.enrich.pvalue.treatmentAll <- ggplot(
    df.toplot %>% ## regenerate df.toplot for p.enrich.pvalue.full
    dplyr::filter(treatmentDEG == "all") %>%
    dplyr::mutate(GeneRatio = sapply(GeneRatio, function(x){eval(parse(text = x))}),
                  NnormRatio = sapply(1:nrow(.), function(i){
                      .[[i,"N.norm"]]/as.numeric(str_split(.[i,"GeneRatio"], '/')[[1]][[2]])
                  }))) +
    geom_point(aes(x = NnormRatio, y = tidytext::reorder_within(GO.name, N.norm,
                                                                list(universe, group, treatmentDEG)),
                   colour = p.adjust, size = N.norm), pch = 16) +
    scale_colour_gradient(high = "blue", low = "red") +
    tidytext::scale_y_reordered(labels = function(x){
        str_wrap(tidytext::reorder_func(x), width = axis.lab.width)}) +
    ## facet_wrap(~universe*has.blast*group, scales = "free", nrow = 3) +
    ggh4x::facet_nested(group~treatmentDEG*mapped.by*universe, independent = "all", scales = "free") +
    theme_bw() + ylab("GO term") +
    theme(axis.text.x = element_text(angle = 90))

save_plot("FigureSX_Fig7supp_map_method_overview_best10pvalue_treatmentAll",
          p.enrich.pvalue.treatmentAll +
          theme(axis.title.y = element_blank()) +
          ggtitle("GO enrichment with different mapping methods & universes",
                  subtitle = paste0("facets top to bottom: treatment DEG, map method, universe\n",
                                    "GO:0009987 (cellular process) dropped; p < ", max.pvalue)),
          w = 13, h = 15, fmt = c("png", "pdf"),
          dir = mkpath(dir_proj, "plots/figures"))

p.enrich.pvalue.treatmentDEG <- ggplot(
    df.toplot %>% ## regenerate df.toplot for p.enrich.pvalue.full
    dplyr::filter(treatmentDEG == "treatmentDEG") %>%
    dplyr::mutate(GeneRatio = sapply(GeneRatio, function(x){eval(parse(text = x))}),
                  NnormRatio = sapply(1:nrow(.), function(i){
                      .[[i,"N.norm"]]/as.numeric(str_split(.[i,"GeneRatio"], '/')[[1]][[2]])
                  }))) +
    geom_point(aes(x = NnormRatio, y = tidytext::reorder_within(GO.name, N.norm,
                                                                list(universe, group, treatmentDEG)),
                   colour = p.adjust, size = N.norm), pch = 16) +
    scale_colour_gradient(high = "blue", low = "red") +
    tidytext::scale_y_reordered(labels = function(x){
        str_wrap(tidytext::reorder_func(x), width = axis.lab.width)}) +
    ## facet_wrap(~universe*has.blast*group, scales = "free", nrow = 3) +
    ggh4x::facet_nested(group~treatmentDEG*mapped.by*universe, independent = "all", scales = "free") +
    theme_bw() + ylab("GO term") +
    theme(axis.text.x = element_text(angle = 90))

save_plot("FigureSX_Fig7supp_map_method_overview_best10pvalue_treatmentDEG",
          p.enrich.pvalue.treatmentDEG +
          theme(axis.title.y = element_blank()) +
          ggtitle("GO enrichment with different mapping methods & universes",
                  subtitle = paste0("facets top to bottom: treatment DEG, map method, universe\n",
                                    "GO:0009987 (cellular process) dropped; p < ", max.pvalue)),
          w = 13, h = 15, fmt = c("png", "pdf"),
          dir = mkpath(dir_proj, "plots/figures"))
