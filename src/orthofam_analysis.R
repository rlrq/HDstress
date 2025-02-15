library(tidyverse)

## library(clusterProfiler)

library(AnnotationDbi)
library(GO.db)
library(org.At.tair.db)

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

f_deg <- mkpath(dir_proj, "/data/orthofam", "7555_genes.txt")
f_hddeg <- mkpath(dir_proj, "/data/orthofam", "DEGs_HDvsC.list")
f_accdeg <- mkpath(dir_proj, "/data/orthofam", "accDEGs_4groups.txt")
f_all <- mkpath(dir_proj, "/data/orthofam", "Anno_Orths.txt")
f_z1 <- mkpath(dir_proj, "/data/orthofam", "Brachiifu_orth_Z1V1.1stBlast.ids")
f_orthofam <- mkpath(dir_proj, "/data/orthofam", "genefamily_data.ORTHOFAM.ath-bra.tsv")
## preview a GO term:
## AnnotationDbi::select(GO.db, keys = "GO:0008150", keytype = "GOID", columns = "TERM") %>% head
## AnnotationDbi::select(GO.db, keys = "response to cold", keytype = "TERM", columns = "GOID") %>% head

## parse chiifu-z1 mapping
df.chiifu_z1 <- read.table(f_z1, header = TRUE, sep = '\t')

## parse brassica genes
df.deg <- read.table(f_deg, header = TRUE, sep = '\t')
df.hddeg <- read.table(f_hddeg, header = TRUE, sep = '\t') %>%
    dplyr::left_join()
df.accdeg <- read.table(f_accdeg, header = TRUE, sep = '\t')
df.bra <- read.table(f_all, header = TRUE, sep = '\t', quote = '') %>%
    dplyr::rename(Gene = BrapaLocus) %>%
    dplyr::full_join(df.chiifu_z1 %>% dplyr::rename(z1.id = "Z1V1_ID"), by = c("Gene" = "Chiifu_ID")) %>%
    dplyr::full_join(df.accdeg %>% dplyr::select(Gene, group2) %>% dplyr::rename(group = group2),
                     by = "Gene") %>%
    dplyr::full_join(df.hddeg %>% dplyr::mutate(treatmentDEG = TRUE) %>%
                     dplyr::select(Gene, AtLocus, treatmentDEG), by = c("Gene", "AtLocus")) %>%
    tidyr::replace_na(list(treatmentDEG = FALSE)) %>%
    dplyr::mutate(group = ifelse(group == "non-accDEG", "ECR", group))
v.ath.id <- df.bra %>% dplyr::filter(!is.na(AtLocus)) %>% dplyr::pull(AtLocus)

## read orthofam groups for brassica and arabidopsis
df.of <- read.table(f_orthofam, header = FALSE, sep = '\t', col.names = c("gf_id", "species", "gene_id")) %>%
    dplyr::rename(orthofam = gf_id, gid = gene_id)

## get all relevant ath gene ids
v.ath.id <- c(
    v.ath.id,
    df.of %>% dplyr::filter(species == "ath") %>% dplyr::pull(gid) %>% as.character
) %>% unique()

## get which orthofams are found in which groups
df.of.bygroup <- df.bra %>%
    dplyr::select(Gene, z1.id, group) %>%
    dplyr::left_join(df.of %>% dplyr::filter(species == "bra"), by = c("z1.id" = "gid")) %>%
    dplyr::filter(!is.na(group)) %>%
    dplyr::select(orthofam, group) %>%
    dplyr::distinct()

########################

## try to reproduce Xiaomei's results
get_ath_genes_in_orthofams <- function(orthofams){
    df.of %>%
        dplyr::filter(orthofam %in% orthofams & species == "ath") %>%
        dplyr::pull(gid) %>%
        unique
}

x2 <- clusterProfiler::enrichGO(get_ath_genes_in_orthofams(df.of.bygroup %>% dplyr::filter(group == "ECR") %>% pull(orthofam)), org.At.tair.db, keyType = "TAIR", ont = "BP")
x2 %>% head(10) %>% dplyr::select(-c(geneID))

## ??? can't reproduce all the sesquiterpene and terpene metabolism pathways. need to ask for script.

#########################

## get stats for all relevant ath gene ids


#########################
## GET ARATH GO TERMS  ##
#########################

## get biological process GO terms for ath genes
df.go.arath <- AnnotationDbi::select(org.At.tair.db, keys = v.ath.id, columns = "GOALL", keytype = "TAIR") %>%
    dplyr::filter(ONTOLOGYALL == "BP") %>%
    dplyr::select(-c(EVIDENCEALL)) %>%
    dplyr::distinct() %>%
    dplyr::filter(GOALL != "GO:0008150") ## remove biological_process GO term. we don't need that. it's given

######### ATH ORTHOFAM GO TERMS ############

## map orthofam group to ath genes
df.go.arath.of <- df.of %>%
    dplyr::filter(species == "ath") %>%
    dplyr::select(-c(species)) %>%
    dplyr::left_join(df.go.arath, by = c("gid" = "TAIR")) %>%
    dplyr::filter(!is.na(ONTOLOGYALL))

## get number of ath/bra genes per orthofam group
N.per.orthofam <- df.of %>%
    dplyr::group_by(orthofam, species) %>%
    dplyr::summarise(N.genes = n())

## get number of annotated ath genes per orthofam group
N.annotated.ath.per.orthofam <- df.go.arath.of %>%
    dplyr::select(orthofam, gid) %>%
    dplyr::distinct() %>%
    dplyr::group_by(orthofam) %>%
    dplyr::summarise(N.ath = n()) %>%
    dplyr::ungroup()

## get frequency of each GO term within each orthofam
N.go.ath.per.orthofam <- df.go.arath.of %>%
    dplyr::select(orthofam, GOALL) %>%
    dplyr::group_by(orthofam, GOALL) %>%
    dplyr::summarise(N.ath.go = n()) %>%
    dplyr::ungroup() %>%
    dplyr::left_join(N.annotated.ath.per.orthofam, by = c("orthofam")) %>%
    dplyr::mutate(fraction.ath.go = N.ath.go/N.ath)

## write.table(N.go.ath.per.orthofam,
##             mkpath(dir_proj, "/data/orthofam", "orthofam.ath.db3.16.0.tsv"),
##             col.names = TRUE, sep = '\t', row.names = FALSE, quote = FALSE)

## get orthofams which have annotated ath genes AND contain bra genes
df.orthofam.go.filtered <- N.per.orthofam %>%
    tidyr::spread(species, N.genes, fill = 0) %>%
    dplyr::filter(ath > 0 & bra > 0) %>%
    dplyr::left_join(N.go.ath.per.orthofam, by = c("orthofam")) %>%
    dplyr::filter(N.ath.go > 0)

## vector of orthofams that contain at least 1 ath gene and 1 bra gene,
## AND where at least 1 ath gene is annotated
v.shared.annotated.orthofam <- df.orthofam.go.filtered %>%
    dplyr::pull(orthofam) %>%
    unique

## write.table(mkpath(dir_proj, "/data/orthofam", "orthofam.ath.summary.db3.16.0.tsv"),
##             N.go.ath.per.orthofam)

write.table(N.per.orthofam %>% dplyr::mutate(shared.annotated.orthofam = orthofam %in% v.shared.annotated.orthofam),
            mkpath(dir_proj, "results/enrichment", "orthofam.counts.ath-bra.tsv"),
            col.names = TRUE, sep = '\t', row.names = FALSE, quote = FALSE)

## get normalised freq of GO terms for Z1 genes in orthofam
df.orthofam.go.z1 <- df.orthofam.go.filtered %>%
    dplyr::mutate(N.z1.norm = fraction.ath.go * bra)

################################
##  MAP GO TERMS TO BRASSICA  ##
################################

############## BLASTABLE BRA GENES #################

## match GO terms to brassica genes
df.go.bra.blast <- df.bra %>%
    dplyr::filter(!is.na(AtLocus)) %>%
    dplyr::left_join(df.go.arath, by = c("AtLocus" = "TAIR"), relationship = "many-to-many") %>%
    dplyr::select(Gene, AtLocus, group, GOALL, treatmentDEG)

## get number of blastable bra genes in each group
N.bra.blast.annotated <- df.go.bra.blast %>%
    dplyr::select(Gene, group, treatmentDEG) %>%
    dplyr::distinct() %>%
    dplyr::group_by(group, treatmentDEG) %>%
    dplyr::summarise(N.bra = n()) %>%
    dplyr::ungroup()

## get number of bra genes from each group mapped to each GO term
df.go.bra.blast.summ <- df.go.bra.blast %>%
    dplyr::filter(!is.na(GOALL)) %>%
    dplyr::group_by(group, treatmentDEG, GOALL) %>%
    dplyr::summarise(N.bra = n()) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(N.bra.norm = N.bra)

## write.table(df.go.bra.blast.summ,
##             mkpath(dir_proj, "results/enrichment", "go_summary.brassica.blast.tsv"),
##             col.names = TRUE, sep = '\t', row.names = FALSE, quote = FALSE)

#################### GO TERMS VIA ORTHOFAM ###################

## map bra genes to orthofams
df.bra.of <- df.bra %>%
    dplyr::mutate(has.blast = !is.na(AtLocus)) %>%
    dplyr::select(Gene, group, treatmentDEG, has.blast) %>%
    dplyr::left_join(df.chiifu_z1, by = c("Gene" = "Chiifu_ID")) %>%
    dplyr::rename(z1.id = Z1V1_ID) %>%
    dplyr::left_join(df.of, by = c("z1.id" = "gid")) %>%
    dplyr::select(-c(species)) %>%
    dplyr::mutate(has.z1 = !is.na(z1.id),
                  has.orthofam = !is.na(orthofam),
                  annotated.orthofam = orthofam %in% v.shared.annotated.orthofam)

write.table(df.bra.of %>%
            dplyr::group_by(orthofam, treatmentDEG, group, has.z1, has.orthofam, annotated.orthofam, has.blast) %>%
            dplyr::summarise(N.bra = n()),
            mkpath(dir_proj, "results/enrichment", "orthofam.counts.brassica.tsv"),
            col.names = TRUE, sep = '\t', row.names = FALSE, quote = FALSE)

## get number of brassia genes in each combination of treatmentDEG, group, has.z1, has.orthofam, annotated.orthofam, has.blast
N.bra.of <- df.bra.of %>%
    dplyr::group_by(treatmentDEG, group, has.z1, has.orthofam, annotated.orthofam, has.blast) %>%
    dplyr::summarise(N.bra = n()) %>%
    dplyr::ungroup()

## map orthogroups & GO terms to bra genes (filter for bra genes with annotated orthofams first)
df.go.bra.of <- df.bra.of %>%
    dplyr::filter(annotated.orthofam) %>%
    dplyr::select(Gene, group, orthofam, treatmentDEG, has.blast) %>%
    dplyr::left_join(df.orthofam.go.filtered %>%
                     dplyr::select(orthofam, GOALL, N.ath.go, N.ath, fraction.ath.go),
                     by = c("orthofam"),
                     relationship = "many-to-many")

## get normalised GO term frequency for each GO term within each group
## for unmapped bra genes mappable to orthofams
df.go.bra.of.summ <- df.go.bra.of %>%
    dplyr::select(group, treatmentDEG, has.blast, GOALL, fraction.ath.go) %>%
    dplyr::group_by(group, treatmentDEG, has.blast, GOALL) %>%
    dplyr::summarise(N.bra.norm = sum(fraction.ath.go),
                     N.bra = n()) %>%
    dplyr::ungroup()

## write.table(df.go.bra.of.summ,
##             mkpath(dir_proj, "results/enrichment", "go_summary.brassica.orthofam.tsv"),
##             col.names = TRUE, sep = '\t', row.names = FALSE, quote = FALSE)

## ## map unmappable bra genes to orthofams
## df.bra.unmapped <- df.bra %>%
##     dplyr::filter(is.na(AtLocus)) %>%
##     dplyr::select(Gene, group, treatmentDEG) %>%
##     dplyr::left_join(df.chiifu_z1, by = c("Gene" = "Chiifu_ID")) %>%
##     dplyr::rename(z1.id = Z1V1_ID) %>%
##     dplyr::left_join(df.of, by = c("z1.id" = "gid")) %>%
##     dplyr::select(-c(species)) %>%
##     dplyr::mutate(
##                has.z1 = !is.na(z1.id),
##                annotated.orthofam = orthofam %in% v.shared.annotated.othofam)

## ## get number of unmappable bra genes in each group, and also how many are in orthofams w/ annotated ath genes
## N.bra.unmapped <- df.bra.unmapped %>%
##     dplyr::group_by(treatmentDEG, group, has.z1, annotated.orthofam) %>%
##     dplyr::summarise(N.bra = n()) %>%
##     dplyr::ungroup()
## ## ~87% of unmapped bra genes are mappable to annotated orthofams


################################
##  SUMMARISE ANNOTATABILITY  ##
################################

df.annotated.blast <- df.go.bra.blast %>%
    dplyr::filter(!is.na(GOALL)) %>%
    dplyr::select(Gene) %>%
    dplyr::distinct() %>%
    dplyr::mutate(ann.blast = TRUE)
df.annotated.of <- df.go.bra.of %>%
    dplyr::filter(!is.na(GOALL)) %>%
    dplyr::select(Gene) %>%
    dplyr::distinct() %>%
    dplyr::mutate(ann.of = TRUE)
df.bra.annotated <- df.bra %>%
    dplyr::mutate(has.blast = !is.na(AtLocus)) %>%
    dplyr::left_join(df.bra.of %>% dplyr::select(Gene, orthofam, has.z1, has.orthofam, annotated.orthofam), by = "Gene") %>%
    dplyr::left_join(df.annotated.blast, by = "Gene") %>%
    dplyr::left_join(df.annotated.of, by = "Gene") %>%
    tidyr::replace_na(list(ann.blast = FALSE, ann.of = FALSE))
remove(df.annotated.blast, df.annotated.of)

## get number of annotated bra genes (whether via initial mapping to arabidopsis or via orthofam) in universe by group
N.bra.annotated <- df.bra.annotated %>%
    dplyr::group_by(group, treatmentDEG, has.blast, has.orthofam, ann.blast, ann.of) %>%
    dplyr::summarise(N.bra = n()) %>%
    dplyr::ungroup()

## N.bra.annotated <- rbind(
##     N.bra.mapped.annotated %>% dplyr::mutate(initial.mapped = TRUE),
##     N.bra.unmapped %>% dplyr::filter(annotated.orthofam) %>% dplyr::mutate(initial.mapped = FALSE) %>%
##     dplyr::select(group, N.bra, initial.mapped, treatmentDEG)
## )


## ## map orthogroups & GO terms to unmapped bra genes
## df.go.bra.unmapped <- df.bra.unmapped %>%
##     dplyr::filter(annotated.orthofam) %>%
##     dplyr::select(Gene, group, orthofam, treatmentDEG) %>%
##     dplyr::left_join(df.orthofam.go.filtered %>%
##                      dplyr::select(orthofam, GOALL, N.ath.go, N.ath, fraction.ath.go),
##                      by = c("orthofam"),
##                      relationship = "many-to-many")

## ## get normalised GO term frequency for each GO term within each group
## ## for unmapped bra genes mappable to orthofams
## df.go.bra.unmapped.summ <- df.go.bra.unmapped %>%
##     dplyr::select(group, treatmentDEG, GOALL, fraction.ath.go) %>%
##     dplyr::group_by(group, treatmentDEG, GOALL) %>%
##     dplyr::summarise(N.bra.norm = sum(fraction.ath.go),
##                      N.bra = n())

## ## make the 'universe' from the annotated bra genes (mapped to annotated ath or mapped to orthofams with annotated ath) among these 7555 genes
## df.universe.go <- rbind(
##     df.go.bra.mapped.summ %>% dplyr::mutate(mapped = TRUE, N.bra.norm = N.bra),
##     df.go.bra.unmapped.summ %>% dplyr::mutate(mapped = FALSE)
## ) %>%
##     dplyr::group_by(group, mapped, GOALL, treatmentDEG) %>%
##     dplyr::summarise(universe.N.bra.norm = sum(N.bra.norm),
##                      universe.N.bra = sum(N.bra)) %>%
##     dplyr::ungroup()
## df.universe.go <- rbind(
##     df.universe.go %>% dplyr::mutate(universe = ifelse(mapped, "mapped", "unmapped")) %>%
##     dplyr::select(-c(mapped)),
##     df.universe.go %>% dplyr::group_by(group, GOALL, treatmentDEG,
##                                        universe.N.bra.norm, universe.N.bra) %>%
##     dplyr::summarise(universe.N.bra.norm = sum(universe.N.bra.norm),
##                      universe.N.bra = sum(universe.N.bra)) %>%
##     dplyr::ungroup() %>%
##     dplyr::mutate(universe = "all")
## ) %>%
##     dplyr::filter(!is.na(GOALL))

## df.universe.go.summ <- df.universe.go %>%
##     dplyr::group_by(universe, treatmentDEG, GOALL) %>%
##     dplyr::summarise(universe.N.bra.norm = sum(universe.N.bra.norm),
##                      universe.N.bra = sum(universe.N.bra))



############################
##  GET ENRICHMENT STATS  ##
############################

## df.counts should have columns: group, GO, N, N.norm
## (where N is the number of genes of interest annotated with a given GO term)
## (where N.norm is the normalised number of genes of interest annotated with a given GO term)
## df.N should have columns: group, N (joined w/ df.counts with "group", N renamed to N.total)
## (where N.total is the number of genes of interest that are annotated)
## df.universe.counts should have columns: GO, N, N.norm
## (joined w/ df.counts with "GO"; additional columns optional)
## (N renamed to N.universe, N.norm renamed to N.norm.universe)
## (where N is the number of genes in the background that are annotated)
## (where N.norm is the normalised number of genes in the background that are annotated)
custom_enrich_GO <- function(df.counts, df.N, df.universe.counts, N.universe.total,
                             pAdjustMethod = "BH"){
    ## N.universe.total <- df.universe.counts %>% dplyr::pull(N.norm) %>% sum()
    df.universe.counts <- df.universe.counts %>%
        dplyr::rename(N.universe = N, N.norm.universe = N.norm)
    df.enrich <- df.counts %>%
        dplyr::left_join(df.N %>% dplyr::rename(N.total = N), by = c("group")) %>%
        dplyr::left_join(df.universe.counts, by = c("GO")) %>%
        dplyr::mutate(GeneRatio = paste0(N.norm, "/", N.total),
                      BgRatio = paste0(N.norm.universe, "/", N.universe.total),
                      pvalue = phyper(N.norm-1,
                                      N.norm.universe,
                                      N.universe.total-N.norm.universe,
                                      N.total,
                                      lower.tail = FALSE)) %>%
    dplyr::group_by(group) %>%
    dplyr::mutate(p.adjust = p.adjust(pvalue, method = pAdjustMethod)) %>%
    dplyr::ungroup()
    return(df.enrich)
}

## calculate GeneRatio, BgRatio, and pvalue for these GO terms
## N_annotated_bra_treatment_initial_mapped <- N.bra.annotated %>%
##     dplyr::filter(initial.mapped & treatmentDEG) %>% dplyr::pull(N.bra) %>% sum()
## N_annotated_bra_treatment_initial_unmapped <- N.bra.annotated %>%
##     dplyr::filter(!initial.mapped & treatmentDEG) %>% dplyr::pull(N.bra) %>% sum()
## N_annotated_bra_treatment <- sum(N_annotated_bra_treatment_initial_mapped,
##                                  N_annotated_bra_treatment_initial_unmapped)
## N_annotated_bra_initial_mapped <- N.bra.annotated %>%
##     dplyr::filter(initial.mapped) %>% dplyr::pull(N.bra) %>% sum()
## N_annotated_bra_initial_unmapped <- N.bra.annotated %>%
##     dplyr::filter(!initial.mapped) %>% dplyr::pull(N.bra) %>% sum()
## N_annotated_bra <- sum(N_annotated_bra_initial_mapped, N_annotated_bra_initial_unmapped)

## subject dataframes
df.go.bra.blast.treatment <- df.go.bra.blast.summ %>% dplyr::filter(treatmentDEG) %>%
    dplyr::rename(GO = GOALL, N = N.bra, N.norm = N.bra.norm) %>% dplyr::mutate(has.blast = TRUE)
df.go.bra.blast.all <- df.go.bra.blast.summ %>%
    dplyr::group_by(group, GOALL) %>% dplyr::summarise(N.bra.norm = sum(N.bra.norm), N.bra = sum(N.bra)) %>% dplyr::ungroup() %>%
    dplyr::rename(GO = GOALL, N = N.bra, N.norm = N.bra.norm) %>%
    dplyr::mutate(treatmentDEG = NA, has.blast = TRUE)
df.go.bra.of.treatment <- df.go.bra.of.summ %>% dplyr::filter(treatmentDEG) %>%
    dplyr::group_by(group, GOALL) %>% dplyr::summarise(N.bra.norm = sum(N.bra.norm), N.bra = sum(N.bra)) %>% dplyr::ungroup() %>%
    dplyr::rename(GO = GOALL, N = N.bra, N.norm = N.bra.norm) %>%
    dplyr::mutate(treatmentDEG = TRUE, has.blast = NA)
df.go.bra.of.treatment.merged <- df.go.bra.of.summ %>% dplyr::filter(treatmentDEG) %>%
    dplyr::group_by(GOALL) %>% dplyr::summarise(N.bra.norm = sum(N.bra.norm), N.bra = sum(N.bra)) %>% dplyr::ungroup() %>%
    dplyr::rename(GO = GOALL, N = N.bra, N.norm = N.bra.norm) %>%
    dplyr::mutate(treatmentDEG = TRUE, group = "merged", has.blast = NA)
df.go.bra.of.all <- df.go.bra.of.summ %>%
    dplyr::group_by(group, GOALL) %>% dplyr::summarise(N.bra.norm = sum(N.bra.norm), N.bra = sum(N.bra)) %>% dplyr::ungroup() %>%
    dplyr::rename(GO = GOALL, N = N.bra, N.norm = N.bra.norm) %>%
    dplyr::mutate(treatmentDEG = NA, has.blast = NA)
df.go.bra.of.all.merged <- df.go.bra.of.summ %>%
    dplyr::group_by(GOALL) %>% dplyr::summarise(N.bra.norm = sum(N.bra.norm), N.bra = sum(N.bra)) %>% dplyr::ungroup() %>%
    dplyr::rename(GO = GOALL, N = N.bra, N.norm = N.bra.norm) %>%
    dplyr::mutate(treatmentDEG = NA, group = "merged", has.blast = NA)
df.go.bra.of.initiallyunmapped.merged <- df.go.bra.of.summ %>% dplyr::filter(!has.blast) %>%
    dplyr::group_by(GOALL) %>% dplyr::summarise(N.bra.norm = sum(N.bra.norm), N.bra = sum(N.bra)) %>% dplyr::ungroup() %>%
    dplyr::rename(GO = GOALL, N = N.bra, N.norm = N.bra.norm) %>%
    dplyr::mutate(treatmentDEG = NA, group = "noblast_merged", has.blast = FALSE)
df.go.bra.of.treatment.initiallyunmapped.merged <- df.go.bra.of.summ %>% dplyr::filter(!has.blast & treatmentDEG) %>%
    dplyr::group_by(GOALL) %>% dplyr::summarise(N.bra.norm = sum(N.bra.norm), N.bra = sum(N.bra)) %>% dplyr::ungroup() %>%
    dplyr::rename(GO = GOALL, N = N.bra, N.norm = N.bra.norm) %>%
    dplyr::mutate(treatmentDEG = TRUE, group = "noblast_merged", has.blast = FALSE)

## N dataframes
df.N.bra.blast.annotated.treatment <- N.bra.blast.annotated %>% dplyr::filter(treatmentDEG) %>%
    dplyr::rename(N = N.bra) %>% dplyr::select(-c(treatmentDEG))
df.N.bra.blast.annotated.all <- N.bra.blast.annotated %>%
    dplyr::group_by(group) %>% dplyr::summarise(N.bra = sum(N.bra)) %>% dplyr::ungroup() %>%
    dplyr::rename(N = N.bra)
df.N.bra.of.annotated.treatment <- N.bra.of %>% dplyr::filter(annotated.orthofam & treatmentDEG) %>%
    dplyr::group_by(group) %>% dplyr::summarise(N.bra = sum(N.bra)) %>% dplyr::ungroup() %>%
    dplyr::rename(N = N.bra)
df.N.bra.of.annotated.treatment.merged <- N.bra.of %>% dplyr::filter(annotated.orthofam & treatmentDEG) %>%
    dplyr::group_by(annotated.orthofam) %>% dplyr::summarise(N.bra = sum(N.bra)) %>% dplyr::ungroup() %>%
    dplyr::rename(N = N.bra) %>% dplyr::select(N) %>% dplyr::mutate(group = "merged")
df.N.bra.of.annotated.all <- N.bra.of %>% dplyr::filter(annotated.orthofam) %>%
    dplyr::group_by(group) %>% dplyr::summarise(N.bra = sum(N.bra)) %>% dplyr::ungroup() %>%
    dplyr::rename(N = N.bra)
df.N.bra.of.annotated.all.merged <- N.bra.of %>% dplyr::filter(annotated.orthofam) %>%
    dplyr::group_by(annotated.orthofam) %>% dplyr::summarise(N.bra = sum(N.bra)) %>% dplyr::ungroup() %>%
    dplyr::rename(N = N.bra) %>% dplyr::select(N) %>% dplyr::mutate(group = "merged")
df.N.bra.of.initiallyunmapped.merged <- N.bra.of %>% dplyr::filter(annotated.orthofam & !has.blast) %>%
    dplyr::group_by(annotated.orthofam) %>% dplyr::summarise(N.bra = sum(N.bra)) %>% dplyr::ungroup() %>%
    dplyr::rename(N = N.bra) %>% dplyr::select(N) %>% dplyr::mutate(group = "noblast_merged")
df.N.bra.of.treatment.initiallyunmapped.merged <- N.bra.of %>% dplyr::filter(annotated.orthofam & !has.blast & treatmentDEG) %>%
    dplyr::group_by(annotated.orthofam) %>% dplyr::summarise(N.bra = sum(N.bra)) %>% dplyr::ungroup() %>%
    dplyr::rename(N = N.bra) %>% dplyr::select(N) %>% dplyr::mutate(group = "noblast_merged")

## universe GO dataframes & Ns
N.universe.blast.treatment <- df.bra.annotated %>% dplyr::filter(treatmentDEG & ann.blast) %>% nrow()
df.go.universe.blast.treatment <- df.go.bra.blast.summ %>%
    dplyr::filter(treatmentDEG) %>%
    dplyr::select(-c(treatmentDEG)) %>%
    dplyr::group_by(GOALL) %>%
    dplyr::summarise(N.bra = sum(N.bra), N.bra.norm = sum(N.bra.norm)) %>%
    dplyr::rename(GO = GOALL, N = N.bra, N.norm = N.bra.norm) %>%
    dplyr::mutate(universe = "blast_HDresponsive")
N.universe.blast <- df.bra.annotated %>% dplyr::filter(ann.blast) %>% nrow()
df.go.universe.blast <- df.go.bra.blast.summ %>%
    dplyr::group_by(GOALL) %>%
    dplyr::summarise(N.bra = sum(N.bra), N.bra.norm = sum(N.bra.norm)) %>%
    dplyr::rename(GO = GOALL, N = N.bra, N.norm = N.bra.norm) %>%
    dplyr::mutate(universe = "blast_all")
## N.universe.blast.merged <- N.universe.blast
## df.go.universe.blast.merged <- df.go.bra.blast.summ %>%
##     dplyr::group_by(GOALL) %>%
##     dplyr::summarise(N.bra = sum(N.bra), N.bra.norm = sum(N.bra.norm)) %>%
##     dplyr::rename(GO = GOALL, N = N.bra, N.norm = N.bra.norm) %>%
##     dplyr::mutate(group = "merged", universe = "blast_all_merged")
N.universe.of.treatment <- df.bra.annotated %>% dplyr::filter(treatmentDEG & ann.of) %>% nrow()
df.go.universe.of.treatment <- df.go.bra.of.summ %>%
    dplyr::filter(treatmentDEG) %>%
    dplyr::select(-c(treatmentDEG)) %>%
    dplyr::group_by(GOALL) %>%
    dplyr::summarise(N.bra = sum(N.bra), N.bra.norm = sum(N.bra.norm)) %>%
    dplyr::rename(GO = GOALL, N = N.bra, N.norm = N.bra.norm) %>%
    dplyr::mutate(universe = "orthofam_HDresponsive")
N.universe.of <- df.bra.annotated %>% dplyr::filter(ann.of) %>% nrow()
df.go.universe.of <- df.go.bra.of.summ %>%
    dplyr::group_by(GOALL) %>%
    dplyr::summarise(N.bra = sum(N.bra), N.bra.norm = sum(N.bra.norm)) %>%
    dplyr::rename(GO = GOALL, N = N.bra, N.norm = N.bra.norm) %>%
    dplyr::mutate(universe = "orthofam_all")
## N.universe.of.merged <- N.universe.of
## df.go.universe.of.merged <- df.go.bra.of.summ %>%
##     dplyr::group_by(GOALL) %>%
##     dplyr::summarise(N.bra = sum(N.bra), N.bra.norm = sum(N.bra.norm)) %>%
##     dplyr::rename(GO = GOALL, N = N.bra, N.norm = N.bra.norm) %>%
##     dplyr::mutate(group = "merged", universe = "orthofam_all_merged")
N.universe.of.z1 <- df.orthofam.go.z1 %>% dplyr::select(orthofam, bra) %>% dplyr::distinct() %>% dplyr::pull(bra) %>% sum()
df.go.universe.of.z1 <- df.orthofam.go.z1 %>%
    dplyr::group_by(GOALL) %>%
    dplyr::summarise(N.bra = sum(bra), N.bra.norm = sum(N.z1.norm)) %>%
    dplyr::rename(GO = GOALL, N = N.bra, N.norm = N.bra.norm) %>%
    dplyr::mutate(universe = "orthofam_Z1")

## blast GxE groups, universe = annotated blast HD-responsive bra genes
df.enrich.blast.treatmentDEG.universeHD <- custom_enrich_GO(
    df.go.bra.blast.treatment,
    df.N.bra.blast.annotated.treatment,
    df.go.universe.blast.treatment,
    N.universe.blast.treatment)
## blast GxE groups, universe = annotated blast bra genes
df.enrich.blast.treatmentDEG.universeAll <- custom_enrich_GO(
    df.go.bra.blast.treatment,
    df.N.bra.blast.annotated.treatment,
    df.go.universe.blast,
    N.universe.blast)
## blast G groups, universe = annotated blast bra genes
df.enrich.blast.all.universeAll <- custom_enrich_GO(
    df.go.bra.blast.all,
    df.N.bra.blast.annotated.all,
    df.go.universe.blast,
    N.universe.blast)
## orthofam GxE groups, universe = annotated orthofam HD-responsive bra genes
df.enrich.of.treatmentDEG.universeHD <- custom_enrich_GO(
    df.go.bra.of.treatment,
    df.N.bra.of.annotated.treatment,
    df.go.universe.of.treatment,
    N.universe.of.treatment)
## orthofam GxE groups, universe = annotated orthofam bra genes
df.enrich.of.treatmentDEG.universeAll <- custom_enrich_GO(
    df.go.bra.of.treatment,
    df.N.bra.of.annotated.treatment,
    df.go.universe.of,
    N.universe.of)
## orthofam GxE groups, universe = annotated orthofam Z1 genes
df.enrich.of.treatmentDEG.universeZ1 <- custom_enrich_GO(
    df.go.bra.of.treatment,
    df.N.bra.of.annotated.treatment,
    df.go.universe.of.z1,
    N.universe.of.z1)
## orthofam groups, universe = annotated orthofam bra genes
df.enrich.of.all.universeAll <- custom_enrich_GO(
    df.go.bra.of.all,
    df.N.bra.of.annotated.all,
    df.go.universe.of,
    N.universe.of)
## orthofam groups, universe = annotated orthofam Z1 genes
df.enrich.of.all.universeZ1 <- custom_enrich_GO(
    df.go.bra.of.all,
    df.N.bra.of.annotated.all,
    df.go.universe.of.z1,
    N.universe.of.z1)
## orthofam initially unmapped all groups, universe = annotated orthofam bra genes
df.enrich.of.initiallyunmapped.universeAll <- custom_enrich_GO(
    df.go.bra.of.initiallyunmapped.all,
    df.N.bra.of.initiallyunmapped.all,
    df.go.universe.of,
    N.universe.of
)
## orthofam initially unmapped all groups, universe = annotated Z1 genes
df.enrich.of.initiallyunmapped.universeZ1 <- custom_enrich_GO(
    df.go.bra.of.initiallyunmapped.all,
    df.N.bra.of.initiallyunmapped.all,
    df.go.universe.of.z1,
    N.universe.of.z1
)
## orthofam initially unmapped all HD repsonsive genes, universe = annotated orthofam bra genes
df.enrich.of.treatmentDEG.initiallyunmapped.universeAll <- custom_enrich_GO(
    df.go.bra.of.treatment.initiallyunmapped.merged,
    df.N.bra.of.treatment.initiallyunmapped.merged,
    df.go.universe.of,
    N.universe.of
)
## orthofam initially unmapped all HD repsonsive genes, universe = annotated Z1 genes
df.enrich.of.treatmentDEG.initiallyunmapped.universeZ1 <- custom_enrich_GO(
    df.go.bra.of.treatment.initiallyunmapped.merged,
    df.N.bra.of.treatment.initiallyunmapped.merged,
    df.go.universe.of.z1,
    N.universe.of.z1
)


df.enrich <- rbind(
    df.enrich.blast.treatmentDEG.universeHD,
    df.enrich.blast.treatmentDEG.universeAll,
    df.enrich.blast.all.universeAll,
    df.enrich.of.treatmentDEG.universeHD,
    df.enrich.of.treatmentDEG.universeAll,
    df.enrich.of.treatmentDEG.universeZ1,
    df.enrich.of.all.universeAll,
    df.enrich.of.all.universeZ1,
    df.enrich.of.initiallyunmapped.universeAll,
    df.enrich.of.initiallyunmapped.universeZ1,
    df.enrich.of.treatmentDEG.initiallyunmapped.universeAll,
    df.enrich.of.treatmentDEG.initiallyunmapped.universeZ1
)

write.table(df.enrich %>% dplyr::rename(GO.id = GO) %>%
            ## reorder columns
            dplyr::select(universe, has.blast, treatmentDEG, group, GO.id,
                          N, N.norm, N.total, N.universe, N.norm.universe,
                          GeneRatio, BgRatio, pvalue, p.adjust),
            mkpath(dir_proj, "/results/enrichment", "enrichment.txt"),
            col.names = TRUE, sep = '\t', quote = FALSE, row.names = FALSE)

## df.oi <- df.enrich.unmapped.treatmentDEG.universeTreatment %>% dplyr::filter(group == "CHD")
## x <- df.oi %>% arrange(pvalue) %>% dplyr::select(group, GO, GeneRatio, pvalue, p.adjust, N.norm) %>% head(20)
## x <- df.oi %>% arrange(N.norm) %>% dplyr::select(group, GO, GeneRatio, pvalue, p.adjust, N.norm) %>% tail(20)
## AnnotationDbi::select(GO.db, keys = x %>% pull(GO), keytype = "GOID", columns = "TERM")

## fisher.test(matrix(c(7, 286, 88, 4523), nrow = 2))

## x3 <- clusterProfiler::enrichGO(df.deg %>% dplyr::filter(group == "C" & !is.na(AtLocus)) %>% dplyr::pull(AtLocus), org.At.tair.db, ont="ALL", pvalueCutoff=20, qvalueCutoff=20, keyType = "TAIR")
## p <- enrichplot::dotplot(x3)
## ggsave(mkpath(dir_proj, "tmp.png"), p)


## 627 ath genes were mapped to 2 different brassica genes etc.
## enrichGO only takes the unique set of input gene ids,
## so these 627 will be half as represented as they should be
## also, using enrichGO as-is, the 'universe' is assumed to be annotated ath genes.
## however, it's clear that the ath universe is pretty different from the bra universe
## # A tibble: 4 × 2
##   N.bra N.ath
##   <int> <int>
## 1     1  4133
## 2     2   627
## 3     3    68
## 4     4     3

