library(tidyverse)
library(ontologyIndex)
library(ontologyPlot)
## library(ggthemes) ## for scale_X_colorblind
## library(colorblindr) ## for cvd_grid (use: cvd_grid(p) to view plot under colour-blind conditions)

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
                      dir = mkpath(dir_proj, "plots", "figures"), ...){
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

f_go <- mkpath(dir_root, "/data/GO/go.obo")
go <- get_OBO(f_go)

figure_num <- "4-v3"

f_goterms <- mkpath(dir_proj, "/data/go_terms/", paste0("fig", figure_num, "_GO_full.tsv"))
df.goterms <- read.table(f_goterms, header = TRUE, sep = '\t')
df.goterms.terminal <- df.goterms %>% dplyr::filter(!is.na(order))
v.goterms.terminal <- df.goterms.terminal %>% dplyr::pull(id)

v.goterms.all <- get_ancestors(go, df.goterms$id)
v.goterms.prune <- remove_links(go, v.goterms.all)

############ PRUNED ##################
## prelim plot to figure out which overarching GO terms to choose to highlight (enriched GO terms in cyan)
par(mar = c(0.1, 0.1, 0.1, 0.1))
p.go.tree <- onto_plot(go, v.goterms.prune,
                       fillcolor = ifelse(v.goterms.prune %in% v.goterms.terminal, "cyan", "grey"))
p.go.tree

## GO terms to highlight (same as above)
v.goterms.highlight.name <- c(
    "response to external stimulus",
    "response to other organism",
    "response to abiotic stimulus",
    "response to chemical",
    "response to hormone",
    "regulation of cellular process",
    "metabolic process",
    ## "cellular component organization or biogenesis",
    "cellular component organization"
)
v.goterms.highlight <- go$name[go$name %in% v.goterms.highlight.name] %>% names

## overview plot to show enriched GO terms (cyan) AND selected higher level GO terms (black outline)
par(mar = c(0.1, 0.1, 0.1, 0.1))
p.go.tree.highlight <- onto_plot(go, v.goterms.prune,
                                 fillcolor = ifelse(v.goterms.prune %in% v.goterms.terminal, "cyan", "grey"),
                                 color = ifelse(v.goterms.prune %in% v.goterms.highlight, "black", NA))
plot_name <- paste0("fig", figure_num, "_GO_overview.pruned")
pdf(mkpath(dir_proj, "plots", "go_terms", paste0(plot_name, ".pdf")),
    width = 23.622, height = 13.7795) ## in inches, cuz apparently this doesn't come with a 'units' arg
p.go.tree.highlight
dev.off()
png(mkpath(dir_proj, "plots", "go_terms", paste0(plot_name, ".png")),
    width = 60, height = 35, units = "cm", res = 300)
p.go.tree.highlight
dev.off()


############ FULL ##################
## cuz some terminal GO terms seem to be missing from pruned :( ##

## prelim plot to figure out which overarching GO terms to choose to highlight (enriched GO terms in cyan)
par(mar = c(0.1, 0.1, 0.1, 0.1))
p.go.tree <- onto_plot(go, v.goterms.all,
                       fillcolor = ifelse(v.goterms.all %in% v.goterms.terminal, "cyan", "grey"),
                       fontsize = 40)
p.go.tree

## GO terms to highlight (same as above)
v.goterms.highlight.name <- c(
    "response to external stimulus",
    "response to other organism",
    "response to abiotic stimulus",
    "response to chemical",
    "response to hormone",
    "regulation of cellular process",
    "metabolic process",
    ## "cellular component organization or biogenesis",
    "cellular component organization"
)
v.goterms.highlight <- go$name[go$name %in% v.goterms.highlight.name] %>% names

## overview plot to show enriched GO terms (cyan) AND selected higher level GO terms (black outline)
par(mar = c(0.1, 0.1, 0.1, 0.1))
p.go.tree.highlight <- onto_plot(go, v.goterms.all,
                                 fillcolor = ifelse(v.goterms.all %in% v.goterms.terminal, "cyan", "grey"),
                                 color = ifelse(v.goterms.all %in% v.goterms.highlight, "black", NA))
plot_name <- paste0("fig", figure_num, "_GO_overview.all")
pdf(mkpath(dir_proj, "plots", "go_terms", paste0(plot_name, ".pdf")),
    width = 47.2441, height = 15.748) ## in inches, cuz apparently this doesn't come with a 'units' arg
p.go.tree.highlight
dev.off()
png(mkpath(dir_proj, "plots", "go_terms", paste0(plot_name, ".png")),
    width = 120, height = 40, units = "cm", res = 300)
p.go.tree.highlight
dev.off()

################## TILE #####################
## okay time to make the tile plot
df.goterms.plot <- df.goterms %>%
    dplyr::filter(!is.na(order)) %>%
    dplyr::select(id, name, ancestors, order) %>%
    tidyr::separate_longer_delim(ancestors, delim = ',') %>%
    dplyr::filter(ancestors %in% v.goterms.highlight) %>%
    dplyr::mutate(is.child = TRUE) %>%
    dplyr::right_join(df.goterms.terminal %>%
                      dplyr::select(id, name, order), by = c("id", "name", "order")) %>%
    dplyr::rename(ancestor.id = ancestors) %>%
    dplyr::left_join(df.goterms %>% dplyr::select(id, name) %>% dplyr::rename(ancestor.name = name),
                     by = c("ancestor.id" = "id")) %>%
    dplyr::mutate(name = factor(name, levels = c(df.goterms.terminal %>% dplyr::arrange(order) %>% dplyr::pull(name))),
                  ancestor.name = factor(ancestor.name, levels = v.goterms.highlight.name)) %>%
    tidyr::complete(nesting(id, name, order), nesting(ancestor.id, ancestor.name)) %>%
    dplyr::filter(!is.na(ancestor.id) & !is.na(name)) %>%
    tidyr::replace_na(list(is.child = FALSE))

p.goterms.tile <- ggplot(df.goterms.plot) +
    theme_bw() +
    geom_tile(aes(x = ancestor.name, y = name, fill = is.child), colour = "grey") +
    scale_fill_manual(values = c("TRUE" = "grey40", "FALSE" = "white"),
                      labels = c("TRUE" = "True", "FALSE" = "False"),
                      name = "inherits from") +
    scale_y_discrete(limits = rev, expand = c(0,0)) +
    scale_x_discrete(position = "bottom", expand = c(0,0)) +
    coord_fixed() +
    theme(
        legend.position = "top",
        ## axis.text.x = element_text(angle = 90, vjust = 0, hjust = 0),
        axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1),
        axis.title = element_blank(),
        axis.ticks = element_blank()
    ) +
    guides(fill = guide_legend(title.positio = "top"))

## if figure 3
save_plot("Figure3_subfig_ancestorGO", p.goterms.tile,
          w = 5.5, h = 8, units = "in", fmt = c("pdf", "png"),
          dir = mkpath(dir_proj, "plots", "go_terms"))

## if figure 3-v2
save_plot("Figure3-v2_subfig_ancestorGO", p.goterms.tile,
          w = 5.5, h = 8, units = "in", fmt = c("pdf", "png"),
          dir = mkpath(dir_proj, "plots", "go_terms"))

## if figure 4
save_plot("Figure4_subfig_ancestorGO", p.goterms.tile,
          w = 5.5, h = 6, units = "in", fmt = c("pdf", "png"),
          dir = mkpath(dir_proj, "plots", "go_terms"))

## if figure 4-v2
save_plot("Figure4-v2_subfig_ancestorGO", p.goterms.tile,
          w = 5.5, h = 6, units = "in", fmt = c("pdf", "png"),
          dir = mkpath(dir_proj, "plots", "go_terms"))

## if figure 4-v3
save_plot("Figure4-v3_subfig_ancestorGO", p.goterms.tile,
          w = 5.5, h = 6, units = "in", fmt = c("pdf", "png"),
          dir = mkpath(dir_proj, "plots", "go_terms"))
