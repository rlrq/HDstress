library(tidyverse)

mkpath <- function(...){paste(..., sep = '/')}

machine_name <- Sys.info()[["nodename"]]
if (machine_name == "chaelab-rlrq"){
    dir_root <- "/home/rachelle/OneDrive/NUS/CEY"
} else if (machine_name == "rlrq-home") {
    dir_root <- "/mnt/d/OneDrive_doysd/OneDrive - Default Directory/NUS/CEY"
} else if (machine_name == "chaelab-ws.nus.edu.sg") {
    dir_root <- "/mnt/chaelab/rachelle"
}

dir_proj <- mkpath(dir_root, "/zzOtherzz/XiaoMei/HDstress")


############################
##    PLOT NETWORK W/     ##
##  CANDIDATE MASTER TFs  ##
############################

save_object <- function(obj, fname){
    saveRDS(obj, file = mkpath(dir_proj, "data", "rds", fname))
}

group_colours <- c(ECR = "darkgoldenrod",
                   CR = "darkgoldenrod",
                   C = "green",
                   CHD = "red",
                   cHD = "red",
                   HD = "blue")

networks <- c("GRN" = "GRN_simplify", "keyGRN" = "keyGRN_simplify", "TFTF" = "TFTF_simplify")
groups <- c("ECR", "C", "HD", "CHD")

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

## (we'll start simple-ish with the keyGRN network)
network_name <- "keyGRN"
nodes.highlight <- c("BraA03g055090.3.5C", "BraA01g023140.3.5C", "BraA04g001910.3.5C") ## top-level cyclic nodes
## nodes.highlight <- "BraA03g030860.3.5C" ## this is a child node of BraA01g023140.3.5C

dat.network <- parse_network_v1(network_name, networks = networks, nodes.highlight = nodes.highlight)

######################
##  igraph
######################
library(igraph)

df.edges <- dat.network$df.edges

vector.edges <- df.edges %>%
    split(.$index) %>%
    lapply(function(df){c(df[1,"regulatoryGene"],df[1,"targetGene"])}) %>%
    unlist()

## make graph
g <- make_graph(edges = vector.edges,
                directed = TRUE)

## ## test layout algorithms
## coords <- layout_with_fr(g) ## this is so random :( (but it's real fast)
## coords <- layout_with_kk(g) ## more reproducible (maybe it's deterministic?)
## coords <- layout_with_drl(g) ## this is so skinny lol
## coords <- layout_with_gem(g) ## wow this takes forever; I terminated it for keyGRN
## coords <- layout_with_graphopt(g) ## lol this is just a pom pom or a blob; it's so homogeneous
## coords <- layout_with_lgl(g) ## ...this looks weird, like a human epidermal cell with its cytoskeleton
## coords <- layout_with_mds(g) ## looks...weird-ish? like a daddy long legs spider
## coords <- layout_with_sugiyama(g) ## oh this is cool in the sense that it's clearly hierarchical but uh. too clearly hierarchical. (it is real fast though)

## after testing, we'll go with kk

## for keyGRN
igraph_layout_algo <- "kk"
coords <- layout_with_kk(g) ## good for keyGRN (...apparently less good for GRN. looks like a ball)
igraph_layout_algo <- "mds"
coords <- layout_with_mds(g)

## for GRN
igraph_layout_algo <- "mds"
coords <- layout_with_mds(g) ## takes a while to run decently separates out cHD & ECR for GRN

## for TFTF
igraph_layout_algo <- "kk"
coords <- layout_with_kk(g) ## good for keyGRN (...apparently less good for GRN. looks like a ball)
igraph_layout_algo <- "mds"
coords <- layout_with_mds(g)

save_object(coords, paste0(network, ".igraph.", igraph_layout_algo, ".rds"))


######################
##  tidygraph & ggraph
######################
library(tidygraph) ## for tbl_graph, activate
library(ggraph) ## for ggraph

make_df_tg_descendants <- function(network_name){
    dat.network <- parse_network_v1(network_name, networks = networks, nodes.highlight = nodes.highlight)    
    df.edges.for_tg <- dat.network$df.edges %>%
        dplyr::mutate(from = regulatoryGene, to = targetGene,
                      to.highlight = regulatoryGene %in% dat.network$id.nodes.highlight_desc)
    ## make a bunch of dummy columns for aesthetics (and some extras for misc variables)
    ## that we can hopefully modify even after creating layout
    df.tg.descendants <- tidygraph::tbl_graph(
        nodes = dat.network$df.nodes %>% dplyr::mutate(alpha = NA, colour = NA, fill = NA, size = NA, pch = NA,
                                                       var1 = NA, var2 = NA, var3 = NA),
        edges = df.edges.for_tg %>% dplyr::mutate(alpha = NA, colour = NA, fill = NA, size = NA, lty = NA,
                                                  var1 = NA, var2 = NA, var3 = NA),
        directed = TRUE, node_key = "Gene")
    return(df.tg.descendants)
}

## function to make layout and save
make_ggraph_layout_and_save <- function(df.tg.descendants, network, layout_algo){
    layout <- ggraph::create_layout(df.tg.descendants, layout = layout_algo)
    save_object(layout, paste0(network, ".ggraph.", layout_algo, ".rds"))
}

## keyGRN
network <- "keyGRN"
df.tg.descendants <- make_df_tg_descendants(network)
make_ggraph_layout_and_save(df.tg.descendants, network, "fr")
make_ggraph_layout_and_save(df.tg.descendants, network, "kk")
make_ggraph_layout_and_save(df.tg.descendants, network, "mds")

## GRN
network <- "GRN"
df.tg.descendants <- make_df_tg_descendants(network)
make_ggraph_layout_and_save(df.tg.descendants, network, "fr") ## try a couple and see
make_ggraph_layout_and_save(df.tg.descendants, network, "mds")

## TFTF
network <- "TFTF"
df.tg.descendants <- make_df_tg_descendants(network)
make_ggraph_layout_and_save(df.tg.descendants, network, "fr")
make_ggraph_layout_and_save(df.tg.descendants, network, "kk")
make_ggraph_layout_and_save(df.tg.descendants, network, "mds")


