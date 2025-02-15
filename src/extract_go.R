library(tidyverse)
library(ontologyIndex)

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

f_go <- "/mnt/chaelab/rachelle/data/GO/2023_04/go.obo"
go <- get_OBO(f_go)

f_go_terms <- mkpath(dir_proj, "data", "go_terms", "GOterms_for_Figure3.txt")
df.goterms.terminal <- read.table(f_go_terms, header = FALSE, sep = '\t', col.names = c("name", "id"), stringsAsFactors = FALSE) %>%
    dplyr::mutate(order = 1:nrow(.))

f_go_terms <- mkpath(dir_proj, "data", "go_terms", "GOterms_for_Figure3-v2.txt")
df.goterms.terminal <- read.table(f_go_terms, header = FALSE, sep = '\t', col.names = c("name", "id"), stringsAsFactors = FALSE) %>%
    dplyr::mutate(order = 1:nrow(.))

f_go_terms <- mkpath(dir_proj, "data", "go_terms", "GOterms_for_Figure4.txt")
df.goterms.terminal <- read.table(f_go_terms, header = TRUE, sep = '\t', stringsAsFactors = FALSE) %>%
    dplyr::rename(id = ID, name = Description) %>%
    dplyr::mutate(order = 1:nrow(.))

f_go_terms <- mkpath(dir_proj, "data", "go_terms", "GOterms_for_Figure4-v2.txt")
df.goterms.terminal <- read.table(f_go_terms, header = TRUE, sep = '\t', stringsAsFactors = FALSE) %>%
    dplyr::rename(id = ID, name = Description) %>%
    dplyr::mutate(order = 1:nrow(.))

f_go_terms <- mkpath(dir_proj, "data", "go_terms", "GOterms_for_Figure4-v3.txt")
df.goterms.terminal <- read.table(f_go_terms, header = FALSE, sep = '\t', stringsAsFactors = FALSE,
                                  col.names = c("id", "name")) %>%
    dplyr::mutate(order = 1:nrow(.))

## get all relevant GO terms (id of GO terms themselves and all of their parents)
v.goterms <- df.goterms.terminal %>% pull(id)
v.goterms.ancestors <- lapply(v.goterms, function(x){unlist(go$ancestors[[x]])}) %>% unlist() %>% unique()
v.goterms <- c(v.goterms, v.goterms.ancestors) %>% unique()

## make df of GO terms
df.goterms <- data.frame(id = v.goterms) %>%
    dplyr::mutate(name = go$name[v.goterms],
                  parents = lapply(go$parents[v.goterms], function(x){paste(x,collapse = ',')}) %>% unlist(),
                  children = lapply(go$children[v.goterms], function(x){paste(x,collapse = ',')}) %>% unlist(),
                  ancestors = lapply(go$ancestors[v.goterms], function(x){paste(unlist(x),collapse = ',')}) %>% unlist(),
                  obsolete = go$obsolete[v.goterms]
                  ) %>%
    dplyr::left_join(df.goterms.terminal %>% dplyr::select(id, order), by = c("id"))

## write
write.table(df.goterms, mkpath(dir_proj, "data", "go_terms", "fig3_GO_full.tsv"),
            col.names = TRUE, sep = '\t', row.names = FALSE, quote= FALSE)
write.table(df.goterms, mkpath(dir_proj, "data", "go_terms", "fig3-v2_GO_full.tsv"),
            col.names = TRUE, sep = '\t', row.names = FALSE, quote= FALSE)
write.table(df.goterms, mkpath(dir_proj, "data", "go_terms", "fig4_GO_full.tsv"),
            col.names = TRUE, sep = '\t', row.names = FALSE, quote= FALSE)
write.table(df.goterms, mkpath(dir_proj, "data", "go_terms", "fig4-v2_GO_full.tsv"),
            col.names = TRUE, sep = '\t', row.names = FALSE, quote= FALSE)
write.table(df.goterms, mkpath(dir_proj, "data", "go_terms", "fig4-v3_GO_full.tsv"),
            col.names = TRUE, sep = '\t', row.names = FALSE, quote= FALSE)
