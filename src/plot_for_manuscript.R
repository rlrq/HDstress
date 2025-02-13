library(tidyverse)
library(ggthemes) ## for scale_X_colorblind
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

networks <- c("GRN" = "GRN_simplify", "keyGRN" = "keyGRN_simplify", "TFTF" = "TFTF_simplify")
## this list below is so we can easily make lists of plots; key and value are both network names
list_networks <- names(networks) %>% as.list()
names(list_networks) <- names(networks)

groups <- c("ECR", "C", "HD", "CHD")
groups.case <- c("ECR", "C", "HD", "cHD")
groups.levels.v2 <- c("E_CR", "C", "HD", "CHD")
groups_format_name_v1 <- function(x){ifelse(x=="ECR", "E_CR", x)}
groups_format_case <- function(x){ifelse(x=="CHD", "cHD", x)}

groups_format_name <- groups_format_name_v1
groups.levels <- groups.levels.v2

group_colours <- c(
    ## ECR = "darkgoldenrod",
    ## CR = "darkgoldenrod",
    ECR = "gold3",
    E_CR = "gold3",
    CR = "gold3",
    C = "green",
    CHD = "red",
    cHD = "red",
    HD = "blue"
)

## make a very small network of the TFs in the cyclic node that we're highlighting
## top-level cyclic nodes
l.nodes.highlight <- list(
    c("BraA01g023140.3.5C", "BraA03g055090.3.5C", "BraA04g001910.3.5C"),
    c("BraA04g012800.3.5C", "BraA05g010280.3.5C", "BraA06g037890.3.5C", "BraA07g012300.3.5C",
      "BraA07g031360.3.5C", "BraA08g036220.3.5C", "BraA10g000440.3.5C")
) %>%
    lapply(sort)

##########################
##  NETWORK DISPERSION  ##
##        PLOTS         ##
##########################

## adapted from plot_network_decay.R

gather_sizes <- function(df){
    df %>%
        tidyr::gather("group", "Nnodes", size, size_ECR, size_C, size_HD, size_CHD) %>%
        dplyr::mutate(group = sapply(group, function(s){
            if (s == "size"){"all"}else{str_extract(s, "(?<=size_).+")[[1]]}
        }))
}

## prune_algos <- c("maxOutAtInit", "maxOutAtInit2TF", "maxOutDynamic", "maxOutDynamic2TF")
prune_algos <- c("maxOutAtInit")

## read raw data
df.prune.raw <- data.frame()
df.prune.raw.long <- data.frame()
for (network_name in names(networks)){
    print(network_name)
    network <- networks[[network_name]]
    for (prune_algo in prune_algos){
        print(prune_algo)
        ## extract info from prune_algo string
        prune_order_determination_time <- ifelse(str_detect(prune_algo, "Init"), "Init", "Dynamic")
        prune_order_by_out_degree_to <- ifelse(str_detect(prune_algo, "2TF"), "out2TF", "out2Any")
        ## parse
        f_prune <- mkpath(dir_proj, "results", "network", paste0(network, ".prune.", prune_algo, ".tsv"))
        df.raw <- read.table(f_prune, header = TRUE, sep = '\t') %>%
            dplyr::mutate(graph_id = as.factor(graph_id),
                          parent_graph_id = as.factor(parent_graph_id),
                          prune_algorithm = as.factor(prune_algo),
                          network = network_name,
                          prune_order_determination_time = prune_order_determination_time,
                          prune_order_by_out_degree_to = prune_order_by_out_degree_to)
        df.raw.long <- df.raw %>% gather_sizes()
        ## merge with master list/df
        df.prune.raw <- rbind(df.prune.raw, df.raw)
        df.prune.raw.long <- rbind(df.prune.raw.long, df.raw.long)
    }
}

## reorder factors
group_order <- c("all", "ECR", "C", "HD", "CHD")
reorder_factors <- function(df.prune){
    df.prune <- df.prune %>%
        dplyr::mutate(prune_order_determination_time = factor(prune_order_determination_time,
                                                              levels = c("Init", "Dynamic")),
                      prune_grp = factor(prune_grp, levels = group_order))
    return(df.prune)
}
df.prune.raw <- df.prune.raw %>% reorder_factors()
df.prune.raw.long <- df.prune.raw.long %>% reorder_factors() %>%
    dplyr::mutate(group = factor(group, levels = group_order))

## get largest network at each iteration
df.prune.max_network <- df.prune.raw %>%
    dplyr::group_by(network, prune_algorithm, prune_grp, iteration) %>%
    dplyr::slice_max(order_by = size, n = 1) %>%
    gather_sizes() %>%
    dplyr::mutate(group = factor(group, levels = group_order))

## make df for plotting
df.prune.toplot <- df.prune.max_network %>%
    dplyr::filter(prune_order_determination_time == "Init",
                  prune_order_by_out_degree_to == "out2Any") %>%
    dplyr::mutate(prune_grp = as.character(prune_grp),
                  group = as.character(group)) %>%
    dplyr::mutate(prune_grp = factor(sapply(prune_grp, groups_format_name), levels = c("all", groups.levels.v2)),
                  group = factor(sapply(group, groups_format_name), levels = c("all", groups.levels.v2)))

## plot main
plot_prune_single_network <- function(network, nrow = 2){
    n <- network
    p.prune.withlegend <- ggplot(df.prune.toplot %>% dplyr::filter(network == n)) +
        geom_line(aes(x = iteration, y = Nnodes, colour = prune_grp)) +
        expand_limits(y = 0) +
        xlab("number of nodes removed") + ylab("number of X group nodes in largest network") +
        scale_colour_manual(values = c(c(all = "black"), group_colours)) +
        guides(colour = guide_legend(title = "group removed", override.aes = list(lwd = 4))) +
        theme_bw() # +
    ## theme(legend.position = "bottom") +
    ## facet_wrap(~group, scales = "free_y", nrow = 1)
    p.prune.legend.right <- ggpubr::get_legend(p.prune.withlegend + theme(legend.position = "right"))
    if (nrow == 2){
        layout_matrix <- rbind(c(1, 1, 1, 1, 1, 1, 1),
                               c(1, 1, 1, 1, 1, 2, 2))
    } else if (nrow == 3){
        layout_matrix <- rbind(c(1, 1, 1, 1, 1), c(1, 1, 1, 1, 1), c(1, 1, 1, 1, 1), c(1, 1, 1, 1, 1),
                               c(1, 1, 1, 1, 1), c(1, 1, 1, 1, 1), c(1, 1, 1, 1, 1), c(1, 1, 1, 1, 1),
                               c(1, 1, 1, 2, 2), c(1, 1, 1, 2, 2), c(1, 1, 1, 2, 2), c(1, 1, 1, 2, 2), c(1, 1, 1, 2, 2))
    }
    p.prune <- gridExtra::grid.arrange(p.prune.withlegend + theme(legend.position = "none") +
                                       facet_wrap(~group, nrow = nrow, scales = "free_y"),
                                       p.prune.legend.right,
                                       layout_matrix = layout_matrix)
    return(p.prune)
}

plots_prune.nrow3 <- lapply(list_networks, function(n){plot_prune_single_network(n, nrow = 3)})
plots_prune.nrow2 <- lapply(list_networks, function(n){plot_prune_single_network(n, nrow = 2)})

## plot supplementary
plot_prune_multiple_networks <- function(networks){
    p.prune <- ggplot(df.prune.toplot %>% dplyr::filter(network %in% networks)) +
        geom_line(aes(x = iteration, y = Nnodes, colour = prune_grp)) +
        expand_limits(y = 0) +
        xlab("number of nodes removed") + ylab("number of X group nodes in largest network") +
        scale_colour_manual(values = c(c(all = "black"), group_colours)) +
        guides(colour = guide_legend(title = "group removed", override.aes = list(lwd = 4))) +
        ## scale_color_colorblind() +
        theme_bw() +
        theme(legend.position = "bottom") +
        ggh4x::facet_nested(network~group,
                            scales = "free_y", independent = "y")
    return(p.prune)
}

plots_prune_excludeOne <- lapply(list_networks, function(n){
    plot_prune_multiple_networks(names(list_networks)[names(list_networks) != n])})


#################
##  MASTER TF  ##
##    PLOTS    ##
#################

## adapted from plot_master_TF.R

make_colnames_v1 <- function(grp){
    return (paste(c("child", "descendant", "descendant_fraction"), grp, sep = '_'))
}

gather_stats_v1 <- function(df, groups, make_colnames = function(grp){grp}){
    colnames_to_gather <- c(sapply(groups, make_colnames))
    df <- df %>%
        tidyr::gather("stat_name_full", "stat_value", all_of(colnames_to_gather)) %>%
        dplyr::mutate(descendants = as.integer(descendants),
                      stat_value = as.numeric(stat_value),
                      stat = paste0("stat_", str_extract(stat_name_full, ".+(?=_[^_]+$)")),
                      stat_group = str_extract(stat_name_full, "[^_]+$")) %>%
        dplyr::select(-c(stat_name_full)) %>%
        tidyr::spread("stat", "stat_value") %>%
        dplyr::mutate(stat_group = as.character(stat_group)) %>%
        dplyr::mutate(stat_group = factor(sapply(stat_group, groups_format_name), levels = groups.levels.v2))
    return(df)
}

simplify_multinode_group_v1 <- function(s){
    groups <- sort(unique(str_split(s, ',')[[1]]))
    return(paste(groups, collapse = ','))
}

gather_stats <- gather_stats_v1
make_colnames <- make_colnames_v1
simplify_multinode_group <- simplify_multinode_group_v1

## parse stats for downstream genes
df.descendants <- data.frame()
for (network_name in names(networks)){
    print(network_name)
    network <- networks[[network_name]]
    f_descendants <- mkpath(dir_proj, "results", "network", paste0(network, ".descendantStats.tsv"))
    df.tmp <- read.table(f_descendants, header = TRUE, sep = '\t') %>%
        gather_stats(groups, make_colnames) %>%
        dplyr::mutate(network = network_name,
                      group = sapply(group, simplify_multinode_group),
                      is_pseudonode = is_pseudonode == "True")
    df.descendants <- rbind(df.descendants, df.tmp)
}
df.descendants <- df.descendants %>%
    dplyr::mutate(group = as.factor(group),
                  nid = sapply(nid, function(s){paste(sort(str_split(s, ',')[[1]]), collapse=',')}))

## select nodes/pseudonodes to label
topn <- 2
df.descendants.tolabel <- df.descendants %>%
    dplyr::group_by(stat_group, network) %>%
    dplyr::slice_max(stat_descendant, n = topn) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(label = str_replace_all(nid, ',', ',\n'))

## make dummy legend. I just. cannot anymore with how difficult it is to do with the main plot
p.descendants.for_legend.colour <- ggplot(data = df.descendants,
                                          aes(x = stat_descendant, y = stat_descendant_fraction)) +
    geom_point(aes(colour = group)) +
    scale_colour_manual(values = c(list("C,CHD" = "darkgreen", "CHD,ECR" = "darkorange", "CHD,E_CR" = "darkorange",
                                        "CHD,HD" = "darkviolet", "ECR,HD" = "cyan3", "E_CR,HD" = "cyan3", "None" = "grey"),
                                   group_colours),
                        labels = c("C" = "C", "C,CHD" = "C,CHD", "CHD" = "CHD", "CHD,ECR" = "CHD,E_CR",
                                   "CHD,HD" = "CHD,HD", "ECR" = "E_CR", "ECR,HD" = "E_CR,HD", "HD" = "HD")) +
    theme_bw() + theme(legend.position = "bottom", legend.key.size = unit(0, "npc"))
p.descendants.for_legend.outline <- ggplot(df.descendants %>% dplyr::filter(is_pseudonode)) +
    geom_point(aes(x = stat_descendant, y = stat_descendant_fraction,
                   fill = group, colour = is_pseudonode)) +
    scale_colour_manual(values = c("black"), labels = c("cyclic subnetwork")) +
    guides(fill = FALSE,
           colour = guide_legend(title = "merged cyclic nodes", label = FALSE,
                                 override.aes = list(pch = 1, size = 3))) +
    theme_bw() + theme(legend.position = "right")
p.downstream.legend.colour.nrow1 <- ggpubr::get_legend(p.descendants.for_legend.colour +
                                                       guides(colour = guide_legend(override.aes = list(size = 3), nrow = 1)))
p.downstream.legend.colour.nrow2 <- ggpubr::get_legend(p.descendants.for_legend.colour +
                                                       guides(colour = guide_legend(override.aes = list(size = 3), nrow = 2)))
p.downstream.legend.colour.nrow3 <- ggpubr::get_legend(p.descendants.for_legend.colour +
                                                       guides(colour = guide_legend(override.aes = list(size = 3), nrow = 3)))
p.downstream.legend.outline <- ggpubr::get_legend(p.descendants.for_legend.outline)
p.downstream.legend.outline.nrow1 <- ggpubr::get_legend(p.descendants.for_legend.outline + theme(legend.position = "bottom"))
p.downstream.legend.nrow1 <- gridExtra::grid.arrange(p.downstream.legend.colour.nrow1, p.downstream.legend.outline.nrow1,
                                                     widths = c(3, 1))
p.downstream.legend.nrow2 <- gridExtra::grid.arrange(p.downstream.legend.colour.nrow2, p.downstream.legend.outline,
                                                     widths = c(5,1.5))
p.downstream.legend.nrow3 <- gridExtra::grid.arrange(p.downstream.legend.colour.nrow3, p.downstream.legend.outline,
                                                     widths = c(3,2))
## ...i can't get the legend title for merged cyclic nodes to NOT overlap with the label in PNG :( (was okay in PDF though...)

## ## plot
## p.downstream.full <- ggplot(data = df.descendants, aes(x = stat_descendant, y = stat_descendant_fraction)) +
##     geom_point(data = df.descendants %>% dplyr::filter(!is_pseudonode),
##                aes(colour = group, fill = group),
##                pch = 16, size = 2) +
##     scale_colour_manual(values = c(list("C,CHD" = "darkgreen", "CHD,ECR" = "darkorange",
##                                         "CHD,HD" = "darkviolet", "ECR,HD" = "cyan3", "None" = "grey"),
##                                    group_colours),
##                         drop = FALSE) +
##     scale_fill_discrete(drop = FALSE) +
##     ggnewscale::new_scale_colour() +
##     ggnewscale::new_scale_fill() +
##     geom_point(data = df.descendants %>% dplyr::filter(is_pseudonode),
##                aes(fill = group, colour = is_pseudonode),
##                pch = 21, size = 2) +
##     scale_fill_discrete(drop = FALSE) +
##     ggrepel::geom_text_repel(data = df.descendants.tolabel, aes(label = label), size = 3, force = 2, force_pull = 0.5) +
##     ggh4x::facet_nested(network~stat_group, scales = "free_x", independent = "x") +
##     ylab("fraction of downstream genes belonging to group") +
##     xlab("number of downstream genes belonging to group") +
##     scale_colour_manual(values = c("black"), labels = c("cyclic subnetwork")) +
##     theme_bw() +
##     theme(legend.position = "none")

plot_downstream_multiple_networks <- function(networks, p.legend, text_colour = "black"){
    df.descendants.filtered <- df.descendants %>% dplyr::filter(network %in% networks)
    df.tolabel.filtered <- df.descendants.tolabel %>% dplyr::filter(network %in% networks)
    p.downstream <- ggplot(data = df.descendants.filtered,
                           aes(x = stat_descendant, y = stat_descendant_fraction)) +
        geom_point(data = df.descendants.filtered %>% dplyr::filter(!is_pseudonode),
                   aes(colour = group, fill = group),
                   pch = 16, size = 2) +
        scale_colour_manual(values = c(list("C,CHD" = "darkgreen", "CHD,ECR" = "darkorange",
                                            "CHD,HD" = "darkviolet", "ECR,HD" = "cyan3", "None" = "grey"),
                                       group_colours),
                            drop = FALSE) +
        scale_fill_discrete(drop = FALSE) +
        ggnewscale::new_scale_colour() +
        ggnewscale::new_scale_fill() +
        geom_point(data = df.descendants.filtered %>% dplyr::filter(is_pseudonode),
                   aes(fill = group, colour = is_pseudonode),
                   pch = 21, size = 2) +
        scale_fill_manual(values = c(list("C,CHD" = "darkgreen", "CHD,ECR" = "darkorange",
                                          "CHD,HD" = "darkviolet", "ECR,HD" = "cyan3", "None" = "grey"),
                                     group_colours),
                          drop = FALSE) +
        ggrepel::geom_text_repel(data = df.tolabel.filtered, aes(label = label),
                                 size = 3, force = 2, force_pull = 0.5, colour = text_colour) +
        ggh4x::facet_nested(network~stat_group, scales = "free_x", independent = "x") +
        ylab("fraction of downstream genes belonging to group") +
        xlab("number of downstream genes belonging to group") +
        scale_colour_manual(values = c("black"), labels = c("cyclic subnetwork")) +
        theme_bw() +
        theme(legend.position = "none")
    return(gridExtra::grid.arrange(p.downstream, p.legend, heights = c(10,1)))
}

plot_downstream_single_network_nolegend <- function(network, nrow, text_colour = "black"){
    n <- network
    df.downstream.main <- df.descendants %>% dplyr::filter(network == n)
    df.downstream.main.tolabel <- df.descendants.tolabel %>% dplyr::filter(network == n)
    p.downstream.main.nolegend <- ggplot(data = df.downstream.main,
                                         aes(x = stat_descendant, y = stat_descendant_fraction)) +
        geom_point(data = df.downstream.main %>% dplyr::filter(!is_pseudonode),
                   aes(colour = group, fill = group),
                   pch = 16, size = 2) +
        scale_colour_manual(values = c(list("C,CHD" = "darkgreen", "CHD,ECR" = "darkorange",
                                            "CHD,HD" = "darkviolet", "ECR,HD" = "cyan3", "None" = "grey"),
                                       group_colours),
                            drop = FALSE) +
        scale_fill_discrete(drop = FALSE) +
        ggnewscale::new_scale_colour() +
        ggnewscale::new_scale_fill() +
        geom_point(data = df.downstream.main %>% dplyr::filter(is_pseudonode),
                   aes(fill = group, colour = is_pseudonode),
                   pch = 21, size = 2) +
        scale_fill_manual(values = c(list("C,CHD" = "darkgreen", "CHD,ECR" = "darkorange",
                                          "CHD,HD" = "darkviolet", "ECR,HD" = "cyan3", "None" = "grey"),
                                     group_colours),
                          drop = FALSE) +
        ggrepel::geom_text_repel(data = df.downstream.main.tolabel, aes(label = label),
                                 size = 3, force = 2, force_pull = 0.5, colour = text_colour) +
        facet_wrap(~stat_group, nrow = nrow, scales = "free_x") +
        ylab("fraction of downstream genes belonging to group") +
        xlab("number of downstream genes belonging to group") +
        scale_colour_manual(values = c("black"), labels = c("cyclic subnetwork")) +
        theme_bw() +
        theme(legend.position = "none")
    return(p.downstream.main.nolegend)
}

## df.downstream.main <- df.descendants %>% dplyr::filter(network == "GRN")
## df.downstream.main.tolabel <- df.descendants.tolabel %>% dplyr::filter(network == "GRN")
## p.downstream.main.nolegend <- ggplot(data = df.downstream.main,
##                                      aes(x = stat_descendant, y = stat_descendant_fraction)) +
##     geom_point(data = df.downstream.main %>% dplyr::filter(!is_pseudonode),
##                aes(colour = group, fill = group),
##                pch = 16, size = 2) +
##     scale_colour_manual(values = c(list("C,CHD" = "darkgreen", "CHD,ECR" = "darkorange",
##                                         "CHD,HD" = "darkviolet", "ECR,HD" = "cyan3", "None" = "grey"),
##                                    group_colours),
##                         drop = FALSE) +
##     scale_fill_discrete(drop = FALSE) +
##     ggnewscale::new_scale_colour() +
##     ggnewscale::new_scale_fill() +
##     geom_point(data = df.downstream.main %>% dplyr::filter(is_pseudonode),
##                aes(fill = group, colour = is_pseudonode),
##                pch = 21, size = 2) +
##     scale_fill_manual(values = c(list("C,CHD" = "darkgreen", "CHD,ECR" = "darkorange",
##                                       "CHD,HD" = "darkviolet", "ECR,HD" = "cyan3", "None" = "grey"),
##                                  group_colours),
##                       drop = FALSE) +
##     ggrepel::geom_text_repel(data = df.downstream.main.tolabel, aes(label = label), size = 3) +
##     facet_wrap(~stat_group, nrow = 2, scales = "free_x") +
##     ylab("fraction of downstream genes belonging to group") +
##     xlab("number of downstream genes belonging to group") +
##     scale_colour_manual(values = c("black"), labels = c("cyclic subnetwork")) +
##     theme_bw() +
##     theme(legend.position = "none")

## ## make subplot with legend
## p.downstream.main <- gridExtra::grid.arrange(p.downstream.main.nolegend,
##                                              p.downstream.legend.nrow3, heights = c(6, 1))

text_colour <- "grey30"
plots_downstream_nolegend.nrow2 <- lapply(
    list_networks, function(n){plot_downstream_single_network_nolegend(n, nrow = 2, text_colour = text_colour)})
plots_downstream.nrow2.nrow3legend <- lapply(
    plots_downstream_nolegend.nrow2, function(p){gridExtra::grid.arrange(p, p.downstream.legend.nrow3, heights = c(6,1))})
plots_downstream_nolegend.nrow1 <- lapply(
    list_networks, function(n){plot_downstream_single_network_nolegend(n, nrow = 1, text_colour = text_colour)})
plots_downstream.nrow1.nrow2legend <- lapply(
    plots_downstream_nolegend.nrow1, function(p){gridExtra::grid.arrange(p, p.downstream.legend.nrow2, heights = c(6,1))})
plots_downstream.nrow1.nrow1legend <- lapply(
    plots_downstream_nolegend.nrow1, function(p){gridExtra::grid.arrange(p, p.downstream.legend.nrow1, heights = c(6,1))})
plots_downstream.excludeOne <- lapply(
    list_networks, function(n){plot_downstream_multiple_networks(names(list_networks)[names(list_networks) != n],
                                                                 p.downstream.legend.nrow1, text_colour = text_colour)})

#####################
##  NETWORK PLOTS  ##
#####################

## adapted from plot_master_TF.R

library(tidygraph) ## for tbl_graph, activate
library(ggraph) ## for ggraph

## nodes.highlight <- sort(c("BraA01g023140.3.5C", "BraA03g055090.3.5C", "BraA04g001910.3.5C")) ## top-level cyclic nodes

parse_network_v1 <- function(network_name, networks = networks, nodes.highlight = nodes.highlight){
    ## parse names
    network <- networks[[network_name]]
    f_edges <- mkpath(dir_proj, "data", "network", paste0(network, ".edges"))
    f_nodes <- mkpath(dir_proj, "data", "network", paste0(network, ".nodes"))
    ## read files
    df.edges <- read.table(f_edges, header = TRUE, sep = '\t') %>%
        tibble::rowid_to_column("index")
    df.nodes <- read.table(f_nodes, header = TRUE, sep = '\t') %>%
        tibble::rowid_to_column("index") %>%
        tidyr::separate(group_module, into = c("group", "module"), sep = '_', remove = FALSE)
    ## nodes & edges descending from nodes to highlight
    df.edges.copy <- df.edges
    id.nodes.highlight_desc <- c()
    current_nodes <- nodes.highlight
    while(length(current_nodes) > 0){
        ## find child edges of current nodes
        df.filtered <- df.edges.copy %>%
            dplyr::filter(regulatoryGene %in% current_nodes)
        ## get child nodes of current nodes
        next_nodes <- df.filtered %>% dplyr::pull(targetGene)
        ## exclude already traversed nodes from child nodes
        next_nodes <- next_nodes[! next_nodes %in% c(id.nodes.highlight_desc, current_nodes)]
        ## remove child edges from df.edges.copy so we don't discover them again
        df.edges.copy <- df.edges.copy %>%
            dplyr::filter(! regulatoryGene %in% current_nodes)
        ## update variables
        id.nodes.highlight_desc <- c(id.nodes.highlight_desc, current_nodes)
        current_nodes <- next_nodes
    }
    id.nodes.highlight_desc <- id.nodes.highlight_desc %>% unique
    return(list(df.edges = df.edges, df.nodes = df.nodes, id.nodes.highlight_desc = id.nodes.highlight_desc))
}

make_df_tg_descendants <- function(network_name, networks, nodes.highlight = c()){
    dat.network <- parse_network_v1(network_name, networks = networks, nodes.highlight = nodes.highlight)
    df.edges.for_tg <- dat.network$df.edges %>%
    dplyr::mutate(from = regulatoryGene, to = targetGene,
                  to.highlight = regulatoryGene %in% dat.network$id.nodes.highlight_desc)
    df.tg.descendants <- tidygraph::tbl_graph(
        nodes = dat.network$df.nodes %>% dplyr::mutate(alpha = NA, colour = NA, fill = NA, size = NA, pch = NA),
        edges = df.edges.for_tg %>% dplyr::mutate(alpha = NA, colour = NA, fill = NA, size = NA, lty = NA),
        directed = TRUE, node_key = "Gene")
    return(list(df.tg.descendants = df.tg.descendants, dat.network = dat.network))
}

deg2rad <- function(d){return(d * pi / 180)}
## df needs columns named x and y; angle should be in degrees
rotate_xy_counter_clockwise <- function(df, angle){
    rad <- deg2rad(angle)
    df <- df %>%
        dplyr::mutate(x = (x*cos(rad)) + (y*sin(rad)),
                      y = (y*cos(rad)) + (x*sin(rad)))
    return(df)
}

## use cytoscape node coords (default prefused force directed layout, exported as cx2, converted to tsv)
cytoscape_algo <- "pfd"
network_names <- names(networks)
## network_names <- "keyGRN"
parsed_graphs <- list()
layouts <- list()
## for each set of nodes to highlight
for (i in 1:length(l.nodes.highlight)){
    parsed_graphs[[i]] <- list()
    layouts[[i]] <- list()
    ## parse graph + layout
    for (network_name in network_names){
        df.coords <- read.table(mkpath(dir_proj, "data", "network",
                                       paste0(networks[[network_name]], ".edges.cx2.cytoscape-",
                                              cytoscape_algo, ".coord")),
                                header = TRUE, sep = '\t')
        ## GRN network has a bunch of much tinier subgraphs that are WAY TOO RIGHT so we'll move them left
        if (network_name == "GRN"){
            bottom_right_min_x <- df.coords %>%
                dplyr::filter(x > 3000 & y < -3000) %>% pull(x) %>% min
            min_x <- df.coords %>% pull(x) %>% min
            df.coords <- df.coords %>%
                dplyr::mutate(x = sapply(1:nrow(.), function(i){
                    if(.[i,"x"] > 3000 && .[i,"y"] < -3000){.[i,"x"] - bottom_right_min_x + min_x}else{.[i,"x"]}
                }))
        } else if (network_name == "keyGRN"){ ## keyGRN needs to be rotated so it largely matches w/ GRN
            df.coords <- df.coords %>% rotate_xy_counter_clockwise(180)
        } else if (network_name == "TFTF"){
            ## nudge bottom right small networks a little left and flip y axis
            df.coords <- df.coords %>%
                dplyr::mutate(x = sapply(1:nrow(.), function(i){
                    if(.[i,"x"] > 1000 && .[i,"y"] < -1000){.[i,"x"] - 1000}else{.[i,"x"]}
                }),
                y = -y)
        }
        parsed_graph <- make_df_tg_descendants(network_name, networks, nodes.highlight = l.nodes.highlight[[i]])
        parsed_graphs[[i]][[network_name]] <- parsed_graph
        dat.network <- parsed_graph$dat.network
        df.tg.descendants <- parsed_graph$df.tg.descendants %>%
            activate(nodes) %>%
            dplyr::left_join(df.coords, by = c("Gene" = "node"))
        layout <- create_layout(df.tg.descendants, layout = "manual",
                                x = df.tg.descendants %>% activate(nodes) %>% pull(x),
                                y = df.tg.descendants %>% activate(nodes) %>% pull(y))
        ## set attributes for aesthetics that we didn't do before creating the layout object
        attributes(layout)$graph <- attributes(layout)$graph %>%
            tidygraph::activate(edges) %>%
            dplyr::mutate(to.highlight = regulatoryGene %in% dat.network$id.nodes.highlight_desc,
                          end_cap_size = sapply(
                              targetGene,
                              function(g){
                                  ifelse(g %in% l.nodes.highlight[[i]], "core",
                                  ifelse(g %in% dat.network$id.nodes.highlight_desc, "desc", "other"))
                              }
                          )) %>%
            tidygraph::activate(nodes) %>%
            dplyr::mutate(size = sapply(
                              Gene,
                              function(g){
                                  ifelse(g %in% l.nodes.highlight[[i]], "core",
                                  ifelse(g %in% dat.network$id.nodes.highlight_desc, "desc", "other"))
                              }
                          ))
        ## this update is required cuz apparently it doesn't propagate??
        layout$size <- attributes(layout)$graph %>% tidygraph::activate(nodes) %>% as_tibble %>% pull(size)
        layouts[[i]][[network_name]] <- layout
    }
}


## plot network
plots.network <- list()
for (network_name in network_names){
    ## plot (any highlighted nodes (i.e. i in layouts[[i]] doesn't matter) is fine since we're not highlighting nodes)
    plots.network[[network_name]] <- ggraph(layouts[[1]][[network_name]]) +
        ## geom_edge_link(show.legend = FALSE, width = 0.1) + ## without arrowheads
        geom_edge_link(show.legend = FALSE, width = 0.1,
                       end_cap = circle(ifelse(network_name == "TFTF", 2, 1), unit = "pt"),
                       arrow = arrow(type = "closed", angle = 20, length = unit(0.01, "npc"))) + ## with arrowheads
        ggnewscale::new_scale_colour() +
        geom_node_point(aes(colour = group), pch = 16, size = ifelse(network_name == "TFTF", 2, 1)) +
        scale_colour_manual(values = group_colours) +
        ## scale_size_manual(values = c("core" = 3, "desc" = 1, "other" = 1)) +
        ## scale_alpha_manual(values = c("core" = 1, "desc" = 1, "other" = 0.3)) +
        theme_void()
}

## plot master TF + downstream
plots.master_highlight <- list()
for (i in 1:length(l.nodes.highlight)){
    plots.master_highlight[[i]] <- list()
    for (network_name in network_names){
        ## plot
        plots.master_highlight[[i]][[network_name]] <- ggraph(layouts[[i]][[network_name]]) +
            ## geom_edge_link(aes(colour = to.highlight), show.legend = FALSE) + ## without arrowheads
            geom_edge_link(aes(colour = to.highlight,
                               end_cap = circle(2*ifelse(end_cap_size == "core", 3, 1), unit = "pt")),
                           show.legend = FALSE,
                           arrow = arrow(type = "closed", angle = 20, length = unit(0.01, "npc"))) + ## with arrowheads
            scale_edge_colour_manual(values = c("black", "grey"), breaks = c("TRUE", "FALSE"), labels = c("TRUE", "FALSE")) +
            ggnewscale::new_scale_colour() +
            geom_node_point(aes(colour = group, size = size, alpha = size),
                            pch = 16) +
            scale_colour_manual(values = group_colours) +
            scale_size_manual(values = c("core" = 3, "desc" = 1, "other" = 1)) +
            scale_alpha_manual(values = c("core" = 1, "desc" = 1, "other" = 0.3)) +
            theme_void()
    }
}

## update layouts with membership at prune iteration so we can plot it
for (network_name in network_names){
    ## plot (any highlighted nodes (i.e. i in layouts[[i]] doesn't matter) is fine since we're not highlighting nodes using those criteria)
    f_prune_membership <- mkpath(dir_proj, "results", "network",
                                 paste0(network_name, "_simplify.prune.membershipAtIter.tsv"))
    df.prune_membership <- read.table(f_prune_membership, header = TRUE, sep = '\t')
    prune_sample_ids <- colnames(df.prune_membership)[colnames(df.prune_membership) != "nid"]
    layout <- layouts[[1]][[network_name]]
    attributes(layout)$graph <- attributes(layout)$graph %>%
        tidygraph::activate(nodes) %>%
        dplyr::left_join(df.prune_membership, by = c("Gene" = "nid"))
    ## iterate through each prune algo/group/iteration combo
    for(prune_sample_id in prune_sample_ids){
        df.tmp <- df.prune_membership[,c("nid", prune_sample_id)]
        colnames(df.tmp) <- c("Gene", "tmp")
        ## add info for edge presence/absence at given iteration
        attributes(layout)$graph <- attributes(layout)$graph %>%
            tidygraph::activate(edges) %>%
            dplyr::left_join(df.tmp %>% dplyr::rename(tmp.A = tmp), by = c("regulatoryGene" = "Gene")) %>%
            dplyr::left_join(df.tmp %>% dplyr::rename(tmp.B = tmp), by = c("targetGene" = "Gene")) %>%
            ## node/link encoding: -1 = absent node/absent link
            ##    0 = node/link is in largest subnetwork
            ##    1 = node/link is not in largest subnetwork
            dplyr::mutate(!!prune_sample_id := ifelse(tmp.A == -1, -1,
                                               ifelse(tmp.B == -1, -1,
                                               ifelse(tmp.A != tmp.B, -1,
                                               ifelse(tmp.A == 0, 0, 1))))) %>%
            dplyr::select(-c(tmp.A, tmp.B))
        ## propagate node info to upper level
        layout[[prune_sample_id]] <- attributes(layout)$graph %>%
            tidygraph::activate(nodes) %>% as_tibble %>% pull(!!prune_sample_id)
    }
    ## update
    layouts[[1]][[network_name]] <- layout
}

## plot prune iter networks
plots.network_prune <- list()
prune_groups <- c("ECR", "CHD", "all")
prune_algos <- "maxOutAtInit"
prune_iters <- c(25, 50)
for (network_name in network_names){
    layout <- layouts[[1]][[network_name]]
    plots.network_prune[[network_name]] <- list()
    for(prune_algo in prune_algos){
        for(prune_group in prune_groups){
            for(prune_iter in prune_iters){
                prune_sample_id <- paste0(prune_algo, '_', prune_group, '_', prune_iter)
                sym_prune_sample_id <- sym(prune_sample_id)
                ## plots.network_prune[[network_name]][[prune_sample_id]] <- "tmp"
                plots.network_prune[[network_name]][[prune_sample_id]] <- ggraph(layout) +
                    ## without arrowheads
                    geom_edge_link(aes(alpha = as.factor(!!sym_prune_sample_id)),
                                   show.legend = FALSE, colour = "black") +
                    ## ## with arrowheads
                    ## geom_edge_link(aes(colour = !!prune_sample_id,
                    ##                    end_cap = circle(2*ifelse(end_cap_size == "core", 3, 1), unit="pt")),
                    ##                show.legend = FALSE,
                    ##                arrow = arrow(type = "closed", angle = 20, length = unit(0.01, "npc"))) +
                    scale_edge_alpha_manual(values = c(0, 1, 0.2),
                                            breaks = c(-1, 0, 1),
                                            labels = c("absent", "largest subnetwork", "other subnetwork")) +
                    ggnewscale::new_scale_colour() +
                    geom_node_point(data = . %>% dplyr::filter(!!sym_prune_sample_id != -1),
                                    aes(colour = group,
                                        alpha = as.factor(ifelse(!!sym_prune_sample_id > 0, 1,
                                                                 !!sym_prune_sample_id))),
                                        ## alpha = as.factor(!!sym_prune_sample_id)),
                                    pch = 16, size = ifelse(network_name == "TFTF", 2, 1)) +
                    geom_node_point(data = . %>% dplyr::filter(!!sym_prune_sample_id == -1),
                                    aes(fill = group),
                                    colour = "black", pch = 21, alpha = 0.5,
                                    size = ifelse(network_name == "TFTF", 2, 1)) +
                    scale_colour_manual(values = group_colours) +
                    scale_fill_manual(values = group_colours, guide = "none") +
                    scale_alpha_manual(values = c(0, 1, 0.3),
                                       breaks = c(-1, 0, 1),
                                       labels = c("absent", "largest subnetwork", "other subnetwork"),
                                       guide = "none") +
                    theme_void()
            }
        }
    }
}

## save to file
prune_groups <- c("ECR", "CHD", "all")
prune_algos <- "maxOutAtInit"
prune_iters <- c(25, 50)
## prune_groups <- c("all")
## prune_network_names <- c("TFTF")
prune_network_names <- network_names
for (network_name in prune_network_names){
    layout <- layouts[[1]][[network_name]]
    for(prune_algo in prune_algos){
        for(prune_group in prune_groups){
            for(prune_iter in prune_iters){
                prune_sample_id <- paste0(prune_algo, '_', prune_group, '_', prune_iter)
                print(prune_sample_id)
                p <- plots.network_prune[[network_name]][[prune_sample_id]]
                save_plot(paste0(network_name, ".prune.membershipAtIter.", prune_sample_id),
                          p + ggtitle(network_name,
                                      subtitle = paste0("group pruned: ", prune_group,
                                                        "; # pruned: ", prune_iter)),
                          h = 5, w = 5,
                          dir = mkpath(dir_proj, "plots", "network_dispersion_at_iter"))
            }
        }
    }
}

############################
##  PLOT MASTER CYCLIC TF
############################
library(tidygraph) ## for tbl_graph, activate
library(ggraph) ## for ggraph

make_cyclic_label_v1 <- function(gid, atg, gene_name){
    if(is.na(gene_name)){
        if(is.na(atg)){
            return(gid)
        } else {
            return(paste0(gid, "\n(", atg, ")"))
        }
    } else {
        if(atg == gene_name){
            return(paste0(gid, "\n(", atg, ")"))
        } else {
            return(paste0(gid, "\n(", atg, "/", gene_name, ")"))
        }
    }
}

## just remake the plots until the jittered positions look about right
layout_algo <- "fr"
## padding <- margin(0,0,1,0,"mm")
padding <- margin(1,1,1,-5,"mm")
make_cyclic_label <- make_cyclic_label_v1
plots.network.cyclic <- list()
i <- 1
for (i in 1:length(l.nodes.highlight)){
    nodes.highlight <- l.nodes.highlight[[i]]
    dat.network <- parse_network_v1("GRN", networks = networks, nodes.highlight = nodes.highlight)
    df.edges.highlight <- dat.network$df.edges %>%
        dplyr::filter(targetGene %in% nodes.highlight)
    df.nodes.highlight <- dat.network$df.nodes %>%
        dplyr::filter(Gene %in% nodes.highlight) %>%
        dplyr::mutate(label = sapply(1:nrow(.), function(i){make_cyclic_label(.[i,"Gene"], .[i,"AtLocus"], .[i,"gene_name"])}))
    df.edges.highlight.for_tg <- df.edges.highlight %>%
        dplyr::left_join(df.nodes.highlight %>% dplyr::select(Gene, label) %>% dplyr::rename(regulatoryGene.label = label),
                         by = c("regulatoryGene" = "Gene")) %>%
        dplyr::left_join(df.nodes.highlight %>% dplyr::select(Gene, label) %>% dplyr::rename(targetGene.label = label),
                         by = c("targetGene" = "Gene")) %>%
        dplyr::mutate(from = regulatoryGene.label, to = targetGene.label,
                      to.highlight = regulatoryGene %in% dat.network$id.nodes.highlight_desc)
    df.tg.highlight <- tidygraph::tbl_graph(
        nodes = df.nodes.highlight %>% dplyr::mutate(alpha = NA, colour = NA, fill = NA, size = NA, pch = NA),
        edges = df.edges.highlight.for_tg %>% dplyr::mutate(alpha = NA, colour = NA, fill = NA, size = NA, lty = NA),
        directed = TRUE, node_key = "label")
    plots.network.cyclic[[i]] <- ggraph(df.tg.highlight, layout = layout_algo) +
        geom_node_label(aes(label = label, colour = group), fontface = "bold", size = 3,
                        position = position_jitter(width = 0.2, height = 0)) +
        geom_edge_link(aes(start_cap = label_rect(regulatoryGene.label, padding = padding),
                           end_cap = label_rect(targetGene.label, padding = padding)),
                       arrow = arrow(type = "closed", angle = 20, length = unit(0.3, "cm"))) +
        scale_x_continuous(expand = c(0.2, 0.2)) +
        scale_y_continuous(expand = c(0.2, 0)) +
        scale_colour_manual(values = group_colours) +
        theme_void() +
        theme(legend.position = "none")
}
## preview plot to inspect jitter
gridExtra::grid.arrange(grobs = plots.network.cyclic, ncol = 1, heights = c(1,2))

## str_highlight <- paste(sort(nodes.highlight), collapse = '-')
## save_plot(paste0("network_cyclic_subgraph.GRN.", str_highlight),
##           p.highlight,
##           w = 5, h = 4, units = "in", fmt = NA, 
##           dir = mkpath(dir_proj, "plots", "network_cyclic_subgraph"))


#################
##  FIGURE 5
#################
label_fontsize <- 15

## MAIN
layout.fig5.v0 <- rbind(c(1, 1),
                        c(1, 1),
                        c(1, 1),
                        c(2, 3),
                        c(2, 3))
layout.fig5.v1 <- rbind(c(1, 1, 1, 1, 1, 2, 2, 2),
                        c(1, 1, 1, 1, 1, 2, 2, 2),
                        c(1, 1, 1, 1, 1, 2, 2, 2),
                        c(1, 1, 1, 1, 1, 4, 4, 4),
                        c(3, 3, 3, 3, 3, 4, 4, 4),
                        c(3, 3, 3, 3, 3, 4, 4, 4))
layout.fig5.v2 <- rbind(c(1, 1, 1, 1, 1, 2, 2, 2, 2),
                        c(1, 1, 1, 1, 1, 2, 2, 2, 2),
                        c(1, 1, 1, 1, 1, 2, 2, 2, 2),
                        c(3, 3, 3, 3, 3, 3, 3, 3, 3),
                        c(3, 3, 3, 3, 3, 3, 3, 3, 3))

grobs.fig5.v0.GRNcentric <- list(
    label_subplot_grob('A', fontsize = label_fontsize, plots.network$GRN + theme(legend.position = "none")),
    label_subplot_grob('B', fontsize = label_fontsize, plots_prune.nrow3[["GRN"]]),
    label_subplot_grob('C', fontsize = label_fontsize, plots_downstream.nrow2.nrow3legend[["GRN"]]))
grobs.fig5.v1.GRNcentric <- list(
    label_subplot_grob('A', fontsize = label_fontsize, plots.network$GRN + theme(legend.position = "none")),
    label_subplot_grob('B', fontsize = label_fontsize, plots.network$keyGRN + theme(legend.position = "none")),
    label_subplot_grob('C', fontsize = label_fontsize, plots_prune.nrow3[["GRN"]]),
    label_subplot_grob('D', fontsize = label_fontsize, plots_downstream.nrow2.nrow3legend[["GRN"]]))
grobs.fig5.v1.keyGRNcentric <- list(
    label_subplot_grob('A', fontsize = label_fontsize, plots.network$GRN + theme(legend.position = "none")),
    label_subplot_grob('B', fontsize = label_fontsize, plots.network$keyGRN + theme(legend.position = "none")),
    label_subplot_grob('C', fontsize = label_fontsize, plots_prune.nrow3[["keyGRN"]]),
    label_subplot_grob('D', fontsize = label_fontsize, plots_downstream.nrow2.nrow3legend[["keyGRN"]]))
grobs.fig5.v2.GRNkeyGRN <- list(
    label_subplot_grob('A', fontsize = label_fontsize, plots.network$GRN + theme(legend.position = "none")),
    label_subplot_grob('B', fontsize = label_fontsize, plots.network$keyGRN + theme(legend.position = "none")),
    label_subplot_grob('C', fontsize = label_fontsize, plots_prune_excludeOne[["TFTF"]] + theme(legend.position = "right")))

## p.fig5 <- plot_fig(grobs.fig5.v1.keyGRNcentric, layout.fig5.v1)
## p.fig5 %>% plot

## without arrowheads
save_plot("Figure5_v0-a-GRNcentric",
          plot_fig(grobs.fig5.v0.GRNcentric, layout.fig5.v0),
          h = 17, w = 10,)
save_plot("Figure5_v1-a-GRNcentric",
          plot_fig(grobs.fig5.v1.GRNcentric, layout.fig5.v1),
          h = 15, w = 15)
save_plot("Figure5_v1-b-keyGRNcentric",
          plot_fig(grobs.fig5.v1.keyGRNcentric, layout.fig5.v1),
          h = 15, w = 15)
save_plot("Figure5_v2-a-GRNkeyGRNdispersion",
          plot_fig(grobs.fig5.v2.GRNkeyGRN, layout.fig5.v2),
          h = 20, w = 20, units = "cm")

## v2-b is v2 with arrowheads
save_plot("Figure5_v2-b-GRNkeyGRNdispersion-arrowhead",
          plot_fig(grobs.fig5.v2.GRNkeyGRN, layout.fig5.v2),
          h = 20, w = 20, units = "cm")


## SUPPLEMENTARY
layout.fig5s.v0 <- rbind(c(1, 1, 3, 3, 3, 3),
                         c(1, 1, 3, 3, 3, 3),
                         c(1, 1, 3, 3, 3, 3),
                         c(1, 1, 4, 4, 4, 4),
                         c(1, 1, 4, 4, 4, 4),
                         c(2, 2, 4, 4, 4, 4),
                         c(2, 2, 4, 4, 4, 4),
                         c(2, 2, 4, 4, 4, 4),
                         c(5, 5, 7, 7, 7, 7),
                         c(5, 5, 7, 7, 7, 7),
                         c(6, 6, 7, 7, 7, 7),
                         c(6, 6, 7, 7, 7, 7),
                         c(6, 6, 7, 7, 7, 7),
                         c(6, 6, 7, 7, 7, 7))
layout.fig5s.v1 <- rbind(c(1, 1, 2, 2, 2, 2),
                         c(1, 1, 2, 2, 2, 2),
                         c(1, 1, 2, 2, 2, 2),
                         c(3, 3, 3, 3, 3, 3),
                         c(3, 3, 3, 3, 3, 3),
                         c(3, 3, 3, 3, 3, 3),
                         c(3, 3, 3, 3, 3, 3),
                         c(3, 3, 3, 3, 3, 3),
                         c(4, 4, 6, 6, 6, 6),
                         c(4, 4, 6, 6, 6, 6),
                         c(5, 5, 6, 6, 6, 6),
                         c(5, 5, 6, 6, 6, 6),
                         c(5, 5, 6, 6, 6, 6),
                         c(5, 5, 6, 6, 6, 6))
layout.fig5s.v2 <- rbind(c(1), c(2))

grobs.fig5s.v0.GRNcentric <- list(
    label_subplot_grob('A', fontsize = label_fontsize, plots.network$keyGRN + theme(legend.position = "none")),
    label_subplot_grob('B', fontsize = label_fontsize, plots.network$TFTF + theme(legend.position = "none")),
    label_subplot_grob('C', fontsize = label_fontsize, plots_prune_excludeOne[["GRN"]]),
    label_subplot_grob('D', fontsize = label_fontsize, plots_downstream.excludeOne[["GRN"]]),
    label_subplot_grob('E', fontsize = label_fontsize, plots.network.cyclic[[1]]),
    label_subplot_grob('F', fontsize = label_fontsize, plots.network.cyclic[[2]]),
    label_subplot_grob('G', fontsize = label_fontsize, plots.master_highlight[[1]][["GRN"]] + theme(legend.position = "none")))
grobs.fig5s.v1.GRNcentric <- list(
    label_subplot_grob('A', fontsize = label_fontsize, plots.network$TFTF + theme(legend.position = "none")),
    label_subplot_grob('B', fontsize = label_fontsize, plots_prune_excludeOne[["GRN"]]),
    label_subplot_grob('C', fontsize = label_fontsize, plots_downstream.excludeOne[["GRN"]]),
    label_subplot_grob('D', fontsize = label_fontsize, plots.network.cyclic[[1]]),
    label_subplot_grob('E', fontsize = label_fontsize, plots.network.cyclic[[2]]),
    label_subplot_grob('F', fontsize = label_fontsize, plots.master_highlight[[1]][["GRN"]] + theme(legend.position = "none")))
grobs.fig5s.v1.keyGRNcentric <- list(
    label_subplot_grob('A', fontsize = label_fontsize, plots.network$TFTF + theme(legend.position = "none")),
    label_subplot_grob('B', fontsize = label_fontsize, plots_prune_excludeOne[["keyGRN"]]),
    label_subplot_grob('C', fontsize = label_fontsize, plots_downstream.excludeOne[["keyGRN"]]),
    label_subplot_grob('D', fontsize = label_fontsize, plots.network.cyclic[[1]]),
    label_subplot_grob('E', fontsize = label_fontsize, plots.network.cyclic[[2]]),
    label_subplot_grob('F', fontsize = label_fontsize, plots.master_highlight[[1]][["keyGRN"]] + theme(legend.position = "none")))
grobs.fig5s.v2.GRNkeyGRN <- list(
    label_subplot_grob('A', fontsize = label_fontsize, plots.network$TFTF + theme(legend.position = "none")),
    label_subplot_grob('B', fontsize = label_fontsize, plots_prune.nrow2[["TFTF"]]))

## without arrowheads
save_plot("FigureS5_v0-a-GRNcentric",
          plot_fig(grobs.fig5s.v0.GRNcentric, layout.fig5s.v0),
          h = 20, w = 15)
save_plot("FigureS5_v1-a-GRNcentric",
          plot_fig(grobs.fig5s.v1.GRNcentric, layout.fig5s.v1),
          h = 20, w = 15)
save_plot("FigureS5_v1-b-keyGRNcentric",
          plot_fig(grobs.fig5s.v1.keyGRNcentric, layout.fig5s.v1),
          h = 20, w = 15)
save_plot("FigureS5_v2-a-GRNkeyGRNdispersion",
          plot_fig(grobs.fig5s.v2.GRNkeyGRN, layout.fig5s.v2),
          h = 21, w = 12, units = "cm")

## v2-b is v2 with arrowheads
save_plot("FigureS5_v2-b-GRNkeyGRNdispersion-arrowhead",
          plot_fig(grobs.fig5s.v2.GRNkeyGRN, layout.fig5s.v2),
          h = 21, w = 12, units = "cm")


#################
##  FIGURE 6
#################
label_fontsize <- 20

## MAIN
layout.fig6.v0 <- rbind(c(1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2),
                        c(1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2),
                        c(1, 1, 1, 1, 1, 1, 1, 1, 3, 3, 3, 3),
                        c(1, 1, 1, 1, 1, 1, 1, 1, 3, 3, 3, 3),
                        c(1, 1, 1, 1, 1, 1, 1, 1, 3, 3, 3, 3),
                        c(1, 1, 1, 1, 1, 1, 1, 1, 3, 3, 3, 3),
                        ## c(4, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 7),
                        c(4, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 7),
                        c(4, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 7),
                        c(4, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 7))

grobs.fig6.v0.GRNkeyGRN <- list(
    label_subplot_grob('A', fontsize = label_fontsize, plots_downstream.excludeOne[["TFTF"]]),
    label_subplot_grob('B', fontsize = label_fontsize, plots.network.cyclic[[1]]),
    label_subplot_grob('C', fontsize = label_fontsize, plots.network.cyclic[[2]]),
    label_subplot_grob('D', fontsize = label_fontsize, plots.master_highlight[[1]][["GRN"]] + theme(legend.position = "none")),
    label_subplot_grob('E', fontsize = label_fontsize, plots.master_highlight[[1]][["keyGRN"]] + theme(legend.position = "none")),
    label_subplot_grob('F', fontsize = label_fontsize, plots.master_highlight[[2]][["GRN"]] + theme(legend.position = "none")),
    label_subplot_grob('G', fontsize = label_fontsize, plots.master_highlight[[2]][["keyGRN"]] + theme(legend.position = "none")))

## ## replace some plots (they need to be regenerated if the jitter doesn't look good)
## grobs.fig6.v0.GRNkeyGRN[[2]] <- label_subplot_grob('B', fontsize = label_fontsize, plots.network.cyclic[[1]])
## grobs.fig6.v0.GRNkeyGRN[[3]] <- label_subplot_grob('C', fontsize = label_fontsize, plots.network.cyclic[[2]])

## ## i'm gonna save the good jitter
saveRDS(grobs.fig6.v0.GRNkeyGRN, mkpath(dir_proj, "data", "rds", "fig6_goodjitter.rds"))
## grobs.fig6.v0.GRNkeyGRN <- readRDS(mkpath(dir_proj, "data", "rds", "network_2024", "fig6_goodjitter.rds"))
grobs.fig6.v0.GRNkeyGRN <- readRDS(mkpath(dir_proj, "data", "rds", "fig6_goodjitter.rds"))
grobs.fig6.v0.GRNkeyGRN[[1]] <- label_subplot_grob('A', fontsize = label_fontsize, plots_downstream.excludeOne[["TFTF"]])
grobs.fig6.v0.GRNkeyGRN[[4]] <- label_subplot_grob('D', fontsize = label_fontsize, plots.master_highlight[[1]][["GRN"]] + theme(legend.position = "none"))
grobs.fig6.v0.GRNkeyGRN[[5]] <- label_subplot_grob('E', fontsize = label_fontsize, plots.master_highlight[[1]][["keyGRN"]] + theme(legend.position = "none"))
grobs.fig6.v0.GRNkeyGRN[[6]] <- label_subplot_grob('F', fontsize = label_fontsize, plots.master_highlight[[2]][["GRN"]] + theme(legend.position = "none"))
grobs.fig6.v0.GRNkeyGRN[[7]] <- label_subplot_grob('G', fontsize = label_fontsize, plots.master_highlight[[2]][["keyGRN"]] + theme(legend.position = "none"))

save_plot("Figure6_v0-a-GRNkeyGRNmaster",
          plot_fig(grobs.fig6.v0.GRNkeyGRN, layout.fig6.v0),
          h = 12, w = 15)

## v0-b is v0 with arrowheads
save_plot("Figure6_v0-b-GRNkeyGRNmaster-arrowhead",
          plot_fig(grobs.fig6.v0.GRNkeyGRN, layout.fig6.v0),
          h = 12, w = 15)


## SUPPLEMENTARY
layout.fig6s.v0 <- rbind(c(1, 1), c(1, 1), c(1, 1), c(1, 1),
                         c(2, 3), c(2, 3), c(2, 3))

grobs.fig6s.v0.GRNkeyGRN <- list(
    label_subplot_grob('A', fontsize = label_fontsize, plots_downstream.nrow1.nrow1legend[["TFTF"]]),
    label_subplot_grob('B', fontsize = label_fontsize, plots.master_highlight[[1]][["TFTF"]] + theme(legend.position = "none")),
    label_subplot_grob('C', fontsize = label_fontsize, plots.master_highlight[[2]][["TFTF"]] + theme(legend.position = "none")))

save_plot("FigureS6_v0-a-GRNkeyGRNmaster",
          plot_fig(grobs.fig6s.v0.GRNkeyGRN, layout.fig6s.v0),
          h = 10, w = 10)

## v0-b is v0 with arrowheads
save_plot("FigureS6_v0-b-GRNkeyGRNmaster-arrowhead",
          plot_fig(grobs.fig6s.v0.GRNkeyGRN, layout.fig6s.v0),
          h = 10, w = 10)


#################
##  FIGURE 7 (tree)
#################
library(ape)
library(ggtree)

f_bra_to_ath <- mkpath(dir_proj, "data/orthofam", "Anno_Orths.txt")

get_gene_v1 <- function(s){
    if(str_detect(s, "^AT\\dG")){unlist(str_extract(s, "^[^.]+"))}
    else if(str_detect(s, "^Brapa_")){s}
    else if(str_detect(s, "^BraA")){
        if(str_detect(s, "\\.3C$")){s}
        else if(str_detect(s, "\\.3\\.5C")){unlist(str_extract(s, "^.+C(?=\\.\\d+$)"))}
    }else{NA}
}

get_gene_v2 <- function(s){
    if(str_detect(s, "^AT\\dG")){unlist(str_extract(s, "^[^.]+"))}
    else if(str_detect(s, "^CAG")){unlist(str_extract(s, "^[^.]+"))}
    else if(str_detect(s, "_Bra")){unlist(str_extract(s, "^[^.]+"))}
    else if(str_detect(s, "^BraA")){
        if(str_detect(s, "\\.3\\.1C$")){s}
        else if(str_detect(s, "\\.3\\.5C")){unlist(str_extract(s, "^.+3.5C(?=\\.\\d+$)"))}
    }else{NA}
}

get_species_v1 <- function(s){
    if(str_detect(s, "^AT\\dG")){"Arabidopsis thaliana"}
    else if(str_detect(s, "^Brapa_")){"Brassica rapa"}
    else if(str_detect(s, "^BraA")){"Brassica rapa"}
    else{NA}
}

get_species_v2 <- function(s){
    if(str_detect(s, "^AT\\dG")){"Arabidopsis thaliana"}
    else if(str_detect(s, "_Bra")){"Brassica rapa"}
    else if(str_detect(s, "^BraA")){"Brassica rapa"}
    else if(str_detect(s, "^CAG")){"Brassica rapa"}
    else{NA}
}

get_accession_v1 <- function(s){
    if(str_detect(s, "^AT\\dG")){"Col-0"}
    else if(str_detect(s, "^Brapa_")){unlist(str_extract(s, "(?<=Brapa_).+(?=_NLR)"))}
    else if(str_detect(s, "^BraA")){"Chiifu"}
    else{NA}
}

get_accession_v2 <- function(s){
    if(str_detect(s, "^AT\\dG")){"Col-0"}
    else if(str_detect(s, "_Bra")){unlist(str_extract(s, "(?<=_Bra)[^/]+"))}
    else if(str_detect(s, "^CAG")){"Z1"}
    else if(str_detect(s, "^BraA")){"Chiifu"}
    else{NA}
}

get_source_v2 <- function(s){
    if(str_detect(s, "AT\\dG")){"TAIR10"}
    else if(str_detect(s, "_Bra")){"Cai2021"}
    else if(str_detect(s, "^CAG")){"Z1"}
    else if(str_detect(s, "^BraA")){
        if(str_detect(s, "\\.3\\.1C$")){"Brara_v3.1"}
        else if(str_detect(s, "\\.3\\.5C")){"Brara_v3.5"}
    }else{NA}
}

get_source_v1 <- function(s){
    if(str_detect(s, "^AT\\dG")){"TAIR10"}
    else if(str_detect(s, "^Brapa_")){"Cai2021"}
    else if(str_detect(s, "^BraA")){
        if(str_detect(s, "\\.3C$")){"Brara_v3.0"}
        else if(str_detect(s, "\\.3\\.5C")){"Brara_v3.5"}
    }else{NA}
}


################## NB-ARC tree ########################

get_gene <- get_gene_v1
get_species <- get_species_v1
get_accession <- get_accession_v1
get_source <- get_source_v1

nwk_nbarc <- mkpath(dir_root, "nlr_survey/results/tree/apt_TCRN", "NB-ARC.Col0-Brapa3-Brapa3.5-R500-Z1-Cai2021.pep.mafft.rooted.nwk")

t.nbarc <- read.tree(nwk_nbarc)
tip.labs <- t.nbarc$tip.label
## get tips to drop (anything with -Protein in it)
tips.drop <- tip.labs[str_detect(tip.labs, "-Protein")]
t.nbarc <- drop.tip(t.nbarc, tips.drop)
tip.labs <- t.nbarc$tip.label
## get NLR clade
mrca.nlr <- getMRCA(t.nbarc, c("Col-0_ref|AT1G59620|CDS|AT1G59620.1|NB-ARC|1|132-382", ## CNL
                               "Col-0_ref|AT5G58120|CDS|AT5G58120.1|NB-ARC|1|191-412")) ## TNL
t.nlr <- extract.clade(t.nbarc, mrca.nlr)

## parse seqids, make columns for source organisms
df.nlr.meta <- data.frame(seqid.nbarc = t.nbarc$tip.label) %>%
    dplyr::mutate(tid = unlist(str_extract(seqid.nbarc, "[^|]+(?=\\|revcomp\\|NB-ARC|\\|NB-ARC)")),
                  gid = sapply(tid, get_gene),
                  source = sapply(tid, get_source),
                  species = sapply(tid, get_species),
                  accession = sapply(tid, get_accession))

## identify which v3.5 genes were mapped by blast
df.bra.ath <- read.table(f_bra_to_ath, header = TRUE, sep = '\t', quote = '')
df.nlr.meta <- df.nlr.meta %>%
    dplyr::left_join(df.bra.ath %>%
                     dplyr::mutate(map.blast = !is.na(AtLocus)) %>%
                     dplyr::select(BrapaLocus, map.blast),
                     by = c("gid" = "BrapaLocus")) %>%
    dplyr::mutate(map.blast = sapply(1:nrow(.), function(i){
        if(.[[i,"source"]] == "Brara_v3.5" & is.na(.[[i,"map.blast"]])){FALSE}else{.[[i,"map.blast"]]}
    }))

## make df mapping seqid to presence of TCRN domains
df.nlr <- df.nlr.meta
domains <- c(TIR = "TIR", CC = "RX-CC_like", RPW8 = "RPW8", "NB-ARC" = "NB-ARC")
for (domain_name in names(domains)){
    domain <- domains[[domain_name]]
    df.seqids <- read.table(mkpath(dir_proj, "results/domain/NB-ARC",
                                   paste0(domain, ".Col0-Brapa3-Brapa3.5-R500-Z1-Cai2021.id.txt")),
                            header = FALSE, sep = '\t', col.names = c("tid", "seqid")) %>%
        dplyr::select(-c(seqid)) %>% dplyr::distinct()
    df.seqids[[domain_name]] <- TRUE
    df.nlr <- df.nlr %>%
        dplyr::left_join(df.seqids, by = "tid")
}
df.nlr <- df.nlr %>%
    tidyr::gather("domain", "present", all_of(names(domains))) %>%
    tidyr::replace_na(list(present = FALSE))

## plot i guess
rect.t.nlr <- ggtree(t.nlr)
rect.t.nlr$data <- rect.t.nlr$data %>%
    dplyr::left_join(df.nlr.meta, by = c("label" = "seqid.nbarc"))
circ.t.nlr <- ggtree(t.nlr, layout = "circular")
circ.t.nlr$data <- circ.t.nlr$data %>%
    dplyr::left_join(df.nlr.meta, by = c("label" = "seqid.nbarc"))
rect.t.nbarc <- ggtree(t.nbarc)
rect.t.nbarc$data <- rect.t.nbarc$data %>%
    dplyr::left_join(df.nlr.meta, by = c("label" = "seqid.nbarc"))

## data for gheatmap plotting of domains & mappability
df.domains.nlr <- df.nlr %>%
    dplyr::select(seqid.nbarc, domain, present) %>%
    dplyr::filter(domain != "NB-ARC" & present) %>%
    dplyr::mutate(present = domain) %>%
    tidyr::spread(domain, present) %>%
    ## dplyr::mutate(seqid.nbarc = factor(seqid.nbarc, levels = t.nbarc$tip.label)) %>%
    tibble::column_to_rownames(var = "seqid.nbarc")
df.mapped.nlr <- df.nlr %>%
    dplyr::select(seqid.nbarc, map.blast) %>%
    dplyr::distinct() %>%
    tibble::column_to_rownames(var = "seqid.nbarc")

## option 1, indicate mappability in ring
p.map.nlr <- gheatmap(
    circ.t.nlr,
    df.mapped.nlr,
    offset = -0.5, width = 0.1, color = NA, colnames = FALSE) +
    scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "blue"),
                      limits = c("TRUE", "FALSE"), na.value = NA,
                      name = "mappable to Ath by BLAST")
p.domain.nlr <- gheatmap(
    p.map.nlr + ggnewscale::new_scale_fill(),
    df.domains.nlr,
    offset = 0, width = 0.2, color = NA, colnames = FALSE,
    colnames_angle = 95, colnames_offset_y = 0.25) +
    scale_fill_discrete(limits = c("TIR", "CC", "RPW8"), na.value = "grey90", name = "domain")

## option 2, indicate mappability in tree
circ.t.nlr.ref.map <- circ.t.nlr +
    geom_tippoint(aes(subset = source %in% c("TAIR10", "Brara_v3.5"),
                      fill = map.blast, shape = source)) +
    scale_shape_manual(values = c("TAIR10" = 24, "Brara_v3.5" = 21)) +
    scale_fill_discrete(na.value = "white")
legend.circ.t.ref.map <- ggpubr::get_legend(
                                     circ.t.nlr.ref.map +
                                     scale_fill_discrete(
                                         na.value = "white",
                                         name = "mappable to Ath by BLAST",
                                         guide = guide_legend(override.aes = list(pch = 21))) +
                                 theme(legend.position = "bottom"))
p.domain.nlr <- gheatmap(
    circ.t.nlr.ref.map + ggnewscale::new_scale_fill(),
    df.domains.nlr,
    offset = -0.5, width = 0.2, color = NA, colnames = FALSE,
    colnames_angle = 95, colnames_offset_y = 0.25) +
    scale_fill_discrete(limits = c("TIR", "CC", "RPW8"), na.value = "grey90", name = "domain")

## option 3, rect tree (NB-ARC), indicate mappability in tree
rect.t.nbarc.ref.map <- rect.t.nbarc +
    geom_tippoint(aes(subset = source %in% c("TAIR10", "Brara_v3.5"),
                      fill = map.blast, shape = source)) +
    scale_shape_manual(values = c("TAIR10" = 24, "Brara_v3.5" = 21)) +
    scale_fill_discrete(na.value = "white")
legend.rect.t.ref.map <- ggpubr::get_legend(
                                     rect.t.nbarc.ref.map +
                                     scale_fill_discrete(
                                         na.value = "white",
                                         name = "mappable to Ath by BLAST",
                                         guide = guide_legend(override.aes = list(pch = 21))) +
                                 theme(legend.position = "bottom", legend.box = "vertical", legend.margin = margin()))
p.rect.domain.nbarc <- gheatmap(
    rect.t.nbarc.ref.map + ggnewscale::new_scale_fill() + theme_tree2(),
    df.domains.nlr,
    offset = 0.1, width = 0.3, color = NA, colnames = FALSE,
    colnames_angle = 95, colnames_offset_y = 0.25) +
    scale_fill_discrete(limits = c("TIR", "CC", "RPW8"), na.value = "grey90", name = "domain") +
    theme(legend.position = "none") +
    scale_x_ggtree()

save_plot("Figure7_subfig_tree_NB-ARC_v3",
          gridExtra::arrangeGrob(p.rect.domain.nbarc + vexpand(.01, direction = 1),
                                 legend.rect.t.ref.map, heights = c(10, 1)),
          h = 10, w = 5, fmt = c("pdf", "png"))

## option 4, rect tree (NB-ARC), indicate mappability in tile
p.rect.map.nbarc <- gheatmap(
    rect.t.nbarc,
    df.mapped.nlr,
    offset = -0.1, width = 0.05, color = NA, colnames = FALSE) +
    scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "blue"),
                      limits = c("TRUE", "FALSE"), na.value = NA,
                      name = "mappable to Ath by BLAST")
p.rect.domain.nbarc <- gheatmap(
    p.rect.map.nbarc + ggnewscale::new_scale_fill() + theme_tree2(),
    df.domains.nlr,
    offset = 0.1, width = 0.3, color = NA, colnames = FALSE,
    colnames_angle = 95, colnames_offset_y = 0.25) +
    scale_fill_discrete(limits = c("TIR", "CC", "RPW8"), na.value = "grey90",
                        name = "domain", guide = "none") +
    theme(legend.position = "bottom") +
    scale_x_ggtree()

save_plot("Figure7_subfig_tree_NB-ARC_v4", p.rect.domain.nbarc + vexpand(.01, direction = 1),
          h = 10, w = 5, fmt = c("pdf", "png"))

## p.t.nbarc <- ggtree(t.nbarc, layout = "circular")

################## PAN_2 tree ########################

get_gene <- get_gene_v2
get_species <- get_species_v2
get_accession <- get_accession_v2
get_source <- get_source_v2

nwk_pan2 <- mkpath(dir_proj, "results/domain/PAN_2/tree", "Col0-Brapa3.5-Z1-Cai2021.PAN_2.pep.mafft.nwk")

t.pan2 <- read.tree(nwk_pan2) %>%
    phytools::midpoint.root()
tip.labs <- t.pan2$tip.label

## ## this part here is just me figuring out which patterns to use for which seqs from which source
## ## these informed v2 of get_X functions
## tip.labs.cai <- tip.labs[str_detect(tip.labs, "_Bra")]
## tip.labs.ath <- tip.labs[str_detect(tip.labs, "AT.G")]
## tip.labs.v35 <- tip.labs[str_detect(tip.labs, "^BraA.+\\.3\\.5C")]
## tip.labs.v31 <- tip.labs[str_detect(tip.labs, "^BraA.+\\.3\\.1C")]
## tip.labs.z1 <- tip.labs[str_detect(tip.labs, "^CAG")]

df.pan2.meta <- data.frame(seqid.pan2 = t.pan2$tip.label) %>%
    dplyr::mutate(tid = unlist(str_extract(seqid.pan2, "^[^/]+")),
                  gid = sapply(tid, get_gene),
                  source = sapply(tid, get_source),
                  species = sapply(tid, get_species),
                  accession = sapply(tid, get_accession)) %>%
    dplyr::left_join(df.bra.ath %>%
                     dplyr::mutate(map.blast = !is.na(AtLocus)) %>%
                     dplyr::select(BrapaLocus, map.blast),
                     by = c("gid" = "BrapaLocus")) %>%
    dplyr::mutate(map.blast = sapply(1:nrow(.), function(i){
        if(.[[i,"source"]] == "Brara_v3.5" & is.na(.[[i,"map.blast"]])){FALSE}else{.[[i,"map.blast"]]}
    }))


## plot i guess
rect.t.pan2 <- ggtree(t.pan2)
rect.t.pan2$data <- rect.t.pan2$data %>%
    dplyr::left_join(df.pan2.meta, by = c("label" = "seqid.pan2"))
circ.t.pan2 <- ggtree(t.pan2, layout = "circular")
circ.t.pan2$data <- circ.t.pan2$data %>%
    dplyr::left_join(df.pan2.meta, by = c("label" = "seqid.pan2"))

## option 1, rect tree (PAN_2), indicate mappability in tree
rect.t.pan2.ref.map <- rect.t.pan2 +
    geom_tippoint(aes(subset = source %in% c("TAIR10", "Brara_v3.5"),
                      fill = map.blast, shape = source)) +
    scale_shape_manual(values = c("TAIR10" = 24, "Brara_v3.5" = 21)) +
    scale_fill_discrete(
        na.value = "white",
        name = "mappable to Ath by BLAST",
        guide = guide_legend(override.aes = list(pch = 21))) +
    theme(legend.position = "bottom", legend.box = "vertical", legend.margin = margin())

save_plot("Figure7_subfig_tree_PAN_2_v1",
          rect.t.pan2.ref.map + vexpand(.01, direction = 1),
          h = 5, w = 5, fmt = c("pdf", "png"))


################## S_locus_glycop tree ########################

get_gene <- get_gene_v2
get_species <- get_species_v2
get_accession <- get_accession_v2
get_source <- get_source_v2

domain_name <- "S_locus_glycop"
nwk_slg <- mkpath(dir_proj, "results/domain", domain_name, "tree",
                  paste0("Col0-Brapa3.5-Z1-Cai2021.", domain_name, ".pep.mafft.nwk"))

t.slg <- read.tree(nwk_slg) %>%
    phytools::midpoint.root()
tip.labs <- t.slg$tip.label

df.slg.meta <- data.frame(seqid.slg = t.slg$tip.label) %>%
    dplyr::mutate(tid = unlist(str_extract(seqid.slg, "^[^/]+")),
                  gid = sapply(tid, get_gene),
                  source = sapply(tid, get_source),
                  species = sapply(tid, get_species),
                  accession = sapply(tid, get_accession)) %>%
    dplyr::left_join(df.bra.ath %>%
                     dplyr::mutate(map.blast = !is.na(AtLocus)) %>%
                     dplyr::select(BrapaLocus, map.blast),
                     by = c("gid" = "BrapaLocus")) %>%
    dplyr::mutate(map.blast = sapply(1:nrow(.), function(i){
        if(.[[i,"source"]] == "Brara_v3.5" & is.na(.[[i,"map.blast"]])){FALSE}else{.[[i,"map.blast"]]}
    }))


## plot i guess
rect.t.slg <- ggtree(t.slg)
rect.t.slg$data <- rect.t.slg$data %>%
    dplyr::left_join(df.slg.meta, by = c("label" = "seqid.slg"))
circ.t.slg <- ggtree(t.slg, layout = "circular")
circ.t.slg$data <- circ.t.slg$data %>%
    dplyr::left_join(df.slg.meta, by = c("label" = "seqid.slg"))

## option 1, rect tree (S_locus_glycop), indicate mappability in tree
rect.t.slg.ref.map <- rect.t.slg +
    geom_tippoint(aes(subset = source %in% c("TAIR10", "Brara_v3.5"),
                      fill = map.blast, shape = source)) +
    scale_shape_manual(values = c("TAIR10" = 24, "Brara_v3.5" = 21)) +
    scale_fill_discrete(
        na.value = "white",
        name = "mappable to Ath by BLAST",
        guide = guide_legend(override.aes = list(pch = 21))) +
    theme(legend.position = "bottom", legend.box = "vertical", legend.margin = margin())

save_plot(paste0("Figure7_subfig_tree_", domain_name, "_v1"),
          rect.t.slg.ref.map + vexpand(.01, direction = 1),
          h = 5, w = 5, fmt = c("pdf", "png"))



