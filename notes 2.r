library(tidyverse)
library(readxl)
library(clipr)

#importa file e riordina le variabili
#importante: 
    # 123 CTRL
    # 456 TRT
####123VS456
df <- as_tibble(read_excel(r"(Original Input/1.1,2,3.vs.4,5,6.xlsx)")) %>%
    transmute(
        name = `gene_id`,
        logFC = as.double(`log2(fold_change)`),
        Pval = as.double(`p_value`),
        logPval = -log(Pval, 10),
        FDR = as.double(`q_value(FDR)`),
        avgCTRL = as.double(value_1),
        avgTRT = as.double(value_2)
                )

#siccome ci sono dei valori indefiniti
# ho assegnato +Inf = +3
# e -Inf = -3
#l'ho fatto guardando i valori massimi e minimi nel df
df$logFC[df$logFC == Inf] <- 3
df$logFC[df$logFC == -Inf] <- -3

####DESEQ
dfdown <- as_tibble(read_excel(r"(Original Input/deseq2_3a8b6bf4670112_0.xlsx)", sheet = 2))
dfup <- as_tibble(read_excel(r"(Original Input/deseq2_3a8b6bf4670112_0.xlsx)", sheet = 1))
df <- rbind(dfup, dfdown) %>%
    transmute(
        name = `Gene name`,
        logFC = logFC,
        Pval = `P-value`,
        logPval = -log(Pval, 10),
        FDR = FDR,
        avgCTRL = `Mean FPKM reference`,
        avgTRT = `Mean FPKM query`
                )


#background di tutti i geni
    #devo prendere solo i geni "accesi" in almeno uno dei due fenotipi
    #per trovare i geni accesi devo fare una media di FPKM, che deve essere > 1
    #ricorda che FPKM è la "quantità" di RNA trovato nella cellula.

#   questo background 

bck <- read_excel(r"(Original Input/20220222_Aurogene_DanieleCarlini_FPKM.xlsx)") %>%
    transmute(
        name = name,
        sample1 = as.double(`1`),
        sample2 = as.double(`2`),
        sample3 = as.double(`3`),
        sample4 = as.double(`4`),
        sample5 = as.double(`5`),
        sample6 = as.double(`6`),
        avgCONTROL = (sample1 + sample2 + sample3)/3,
        avgTRT = (sample4 + sample5 + sample6)/3
            ) %>%
        filter(avgCONTROL > 1 | avgTRT > 1) %>%
        view() %>%
        write_tsv(r"(data_frames/aurogene/background.tsv)")

#solo geni significativi
df %>%
    filter(FDR <= 0.01) %>%
    write_tsv(r"(data_frames/deseq/geni_significativi.tsv)")

df %>%
    filter(logFC < -1, FDR <= 0.01) %>%
    write_tsv(r"(data_frames/deseq/geni_downregulated_extreme.tsv)")



# funzione per scrivere 4 file diversi con liste di geni
separateLists <- function(df) {
    filtered <- filter(df, avgCTRL > 1 | avgTRT > 1)
    filtered %>%
        filter(logFC < -1 & FDR <= 0.01) %>%
            write_tsv(r"(data_frames/123_vs_456/geni_downregulated_extreme.tsv)")

    filtered %>%
        filter(logFC < 0 & FDR <= 0.01) %>%
            write_tsv(r"(data_frames/123_vs_456/geni_downregulated.tsv)")

    filtered %>%
        filter(logFC > 0 & FDR <= 0.01) %>%
            write_tsv(r"(data_frames/123_vs_456/geni_upregulated.tsv)")

    filtered %>%
        filter(logFC > +1 & FDR <= 0.01) %>%
            write_tsv(r"(data_frames/123_vs_456/geni_upregulated_extreme.tsv)")
}

separateLists(df)


#   VOLCANO PLOT: Sotto c'è la Funzione
# df %>% filter(
#     avgCONTROL >= 1 | avgTRT >= 1) %>%
#     mutate(
#     significant = ifelse(`FDR` <= 0.01, "yes", "no")) %>%
#     ggplot(aes(
#     x = `logFC`,
#     y = `logPval`,
#     color = `significant`)) +
#     geom_point(
#         alpha = 0.1
#     ) +
#     scale_color_manual(
#         values = c("black", "red")
#     ) +
#     coord_cartesian(
#         xlim = c(-5,5),
#         ylim = c(0,4)
#     )


volcano_plot <- function(df, xlim = c(-5,5), ylim = c(0,5)){
    df %>% filter(
        avgCTRL >= 1 | avgTRT >= 1) %>%
            mutate(
            significant = ifelse(`FDR` <= 0.01, "yes", "no")) %>%
                ggplot(aes(
                x = `logFC`,
                y = `logPval`,
                color = `significant`)) +
                geom_point(
                    alpha = 0.1
                ) +
                scale_color_manual(
                    values = c("black", "red")
                ) +
                coord_cartesian(
                    xlim = xlim,
                    ylim = ylim
                )
}

df <- read_tsv(r"(data_frames/deseq/geni_upregulated_extreme.tsv)")
write_clip(df$name)


files <- list.files(r"(data_frames/deseq/Molecular Function/files)",full.names=TRUE, recursive=FALSE)
files

lapply(files, function(x) {
    df <- as_tibble(read_delim(x, skip = 11)) %>%
        transmute(
            name = `GO molecular function complete`,
            sign = ifelse(`upload_1 (over/under)` == "+", as.double(+1), as.double(-1)),
            fold_enrichment = ifelse(is.numeric(`upload_1 (fold Enrichment)`) == TRUE, 
                                    as.double(`upload_1 (fold Enrichment)`, 
                                    ifelse(`upload_1 (fold Enrichment)` == ">100", 
                                        as.double(100),
                                        as.double(0))),
            fold_enrichment = fold_enrichment * sign,
            FDR = as.double(`upload_1 (FDR)`),
            Pval = `upload_1 (raw P-value)`
        ) %>%
        select(-sign) %>%
        filter(FDR <= 0.01) %>%
        ggplot(
                mapping = aes(
                    x = name,
                    y = fold_enrichment)) +
                geom_col(
                    aes(fill = fold_enrichment)) +
                scale_fill_gradient(
                    low = "#5e8319",
                    high = "#f32121",
                    space = "Lab",
                    #na.value = "grey50",
                    guide = "colourbar",
                    aesthetics = "fill") +
                coord_flip() +
                theme(legend.position = "none") +
                labs(title = "Downregulated (Molecular Function)") %>%


    png("Downregulated_MolecularFunx.png")
    print(x)
    dev.off()
    })





#problema:
#risolto

df <- as_tibble(
    read_delim(
        r"(data_frames/deseq/Cellular Component/GOEA_upregulated_extreme_deseq.txt)", 
        skip = 11)) %>%
    transmute(
            name = `GO cellular component complete`,
            sign = ifelse(`upload_1 (over/under)` == "+", as.double(+1), as.double(-1)),
            fold_enrichment = as.double(`upload_1 (fold Enrichment)`) * sign,
            FDR = as.double(`upload_1 (FDR)`),
            Pval = `upload_1 (raw P-value)`
        ) %>%
        select(-sign) %>%
        arrange(fold_enrichment) %>%
        mutate(name = fct(name))


p <- df %>%
    filter(FDR <= 0.01) %>%
    ggplot(
            mapping = aes(
                x = `name`,
                y = fold_enrichment)) +
            geom_col(
                aes(fill = fold_enrichment)) +
            scale_fill_gradient(
                low = "#5e8319",
                high = "#f32121",
                space = "Lab",
                #na.value = "grey50",
                guide = "colourbar",
                aesthetics = "fill") +
            coord_flip() +
            theme(legend.position = "none") +
            labs(title = "Extreme Upregulated (Cellular Comp)")
p


png(filename = r"(data_frames/deseq/Cellular Component/Upregulated_Extreme_CellularComp.png)", width = 1080, height = 1080)
print(last_plot())
dev.off()






#

df <- transmute(df,
    name = gene,
    avgCONTROL = as.double(value_1),
    avgTRT = as.double(value_2),
    logFC = as.double(`log2(fold_change)`),
    Pval = as.double(p_value),
    logPval = -log(Pval, 10),
    FDR = as.double(`q_value(FDR)`))




#GOEA programmatic access

df <- read_tsv(r"(data_frames/123_vs_456/geni_downregulated_extreme.tsv)")


enriched <- df %>%
                transmute(
                    gene = name,
                    `0/1` = 1
                ) %>%
                as.data.frame() %>%
                go_enrich()


Go_Enrich_Out<- go_enrich(df)

input <- tibble(
    gene = bck$name,
    `0/1` = ifelse(bck$name %in% df$gene, 1, 0)
    )

enriched <- go_enrich(as.data.frame(input))

df <- enriched$results
df <- df %>%
        transmute(
            category = ontology,
            
        )
