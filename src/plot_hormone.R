library(tidyverse)
library(ggthemes) ## for scale_X_colorblind
library(ggh4x)
library(colorblindr) ## for cvd_grid (use: cvd_grid(p) to view plot under colour-blind conditions)

extrafont::font_import(pattern = "Arial", prompt = FALSE)
theme_set(theme_bw(base_family = "Arial", base_size = 8))
theme_set(theme_void(base_family = "Arial", base_size = 8))
update_geom_defaults("text", list(size = 8/.pt))

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
                      dir = mkpath(dir_proj, "plots", "hormone"), ...){
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

label_subplot_grob <- function(lab, p, fontsize = 10){
    gridExtra::arrangeGrob(grobs = list(plot = p),
                           left = grid::textGrob(rlang::expr(bold(!!lab)),
                                                 vjust = 2, x = unit(0.5, "npc"), y = unit(1, "npc"),
                                                 gp = grid::gpar(fontsize = fontsize)))
}

## applies layout matrix to grobs, then adds a white background
plot_fig <- function(grobs, layout_matrix, ...){
    p <- gridExtra::grid.arrange(grobs = grobs, layout_matrix = layout_matrix, ...)
    return(cowplot::ggdraw(p) + theme(plot.background = element_rect(fill = "white", colour = NA)))
}

###############################

f_aba_sa <- mkpath(dir_proj, "data", "hormone_SA_ABA.tsv")
df.abasa <- read.table(f_aba_sa, header = TRUE, sep = '\t') %>%
    tidyr::separate(condition, into = c("accession", "timept", "condition"), sep = '_') %>%
    dplyr::mutate(time = as.numeric(str_extract(timept, "\\d+"))) %>%
    dplyr::mutate(hormone = str_extract(subject, "^[^_]+"),
                  nested_condition = paste0(time, '&', condition),
                  metric_with_units = paste0(metric, " (", units, ")"))

p.overview <- ggplot(df.abasa %>% dplyr::filter(hormone != "Peak")) +
    geom_boxplot(aes(y = value, x = interaction(time, condition), fill = accession)) +
    geom_point(aes(y = value, x = interaction(time, condition), fill = accession),
               shape = 21, position = position_jitterdodge()) +
    scale_x_discrete(guide = "axis_nested") +
    theme_bw() +
    theme(axis.title.x = element_blank()) +
    ## ggh4x::facet_grid2(hormone~metric, scales = "free_y", independent = "y")
    ggh4x::facet_grid2(metric_with_units~hormone, scales = "free_y", independent = "y")

df.abasa.fw <- df.abasa %>% dplyr::filter(hormone != "Peak" & metric == "Absolute (0.4g FW)")
fw.units <- df.abasa.fw[1,"units"]
p.fw <- ggplot(df.abasa.fw) +
    geom_boxplot(aes(y = value, x = interaction(time, condition), fill = accession)) +
    geom_point(aes(y = value, x = interaction(time, condition), fill = accession),
               shape = 21, position = position_jitterdodge()) +
    scale_x_discrete(guide = "axis_nested") +
    theme_bw() +
    ylab(fw.units) +
    theme(axis.title.x = element_blank()) +
    ## ggh4x::facet_grid2(hormone~metric, scales = "free_y", independent = "y")
    ggh4x::facet_grid2(.~hormone)

save_plot("hormone_ABA-SA.overview", p.overview, h = 10, w = 8)
save_plot("hormone_ABA-SA.FW", p.fw, h = 5, w = 6)
