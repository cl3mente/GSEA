library(tidyverse)
library(readxl)
library(clipr)

fig <- function(width, heigth){
     options(repr.plot.width = width, repr.plot.height = heigth)
}

#importa file e riordina le variabili
df <- as_tibble(read_excel(r"(1.1,2,3.vs.4,5,6.xlsx)")) %>%
    transmute(
        name = `gene_id`,
        logFC = as.double(`log2(fold_change)`),
        Pval = as.double(`p_value`),
        logPval = -log(Pval, 10),
        FDR = as.double(`q_value(FDR)`),
                )

dfdown <- as_tibble(read_excel(r"(deseq2_3a8b6bf4670112_0.xlsx)", sheet = 2)) %>%
    transmute(
        name = `Gene name`,
        logFC = logFC,
        Pval = `P-value`,
        logPval = -log(Pval, 10),
        FDR = FDR
                )

dfup <- as_tibble(read_excel(r"(deseq2_3a8b6bf4670112_0.xlsx)", sheet = 1)) %>%
    transmute(
        name = `Gene name`,
        logFC = logFC,
        Pval = `P-value`,
        logPval = -log(Pval, 10),
        FDR = FDR
                ) 
df <- rbind(dfup, dfdown)


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
    write_tsv(r"(data_frames\deseq\geni_downregulated_extreme.tsv)")




#volcano plot
df %>% mutate(
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
        xlim = c(-5,5)
    )


df <- read_tsv(r"(data_frames/deseq/geni_upregulated_extreme.tsv)")
write_clip(df$name)



df <- as_tibble(
    read_delim(
        r"(data_frames\deseq\Cellular Component\GOEA_downregulated_extreme_deseq.txt)", skip = 11)) %>%
    transmute(
        name = `GO cellular component complete`,
        fold_enrichment = as.double(`upload_1 (fold Enrichment)`),
        FDR = as.double(`upload_1 (FDR)`),
        Pval = `upload_1 (raw P-value)`
    )

df[is.na(df)] <- 0
df <- as_tibble(df)

p <- df %>%
    filter(FDR <= 0.01) %>%
    ggplot(df,
            mapping = aes(
                x = name,
                y = fold_enrichment),
                fig(1, 20)) +
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
            labs(title = "Extreme Downregulated (Cellular Comp)")
p

png("Extreme Downregulated_CellularComp.png")
print(last_plot())
dev.off()


