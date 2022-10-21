library(tidyverse)
library(readxl)
library(clipr)

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
bck <- df %>%
    filter(average of control FPKM > 1 | average of treated FPKM > 1)



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
    color = `significant`,
    alpha = .2)) +
    geom_point() +
    coord_cartesian(
        xlim = c(-5,5)
    )


df <- read_tsv(r"(data_frames/deseq/geni_downregulated_extreme.tsv)")
write_clip(df$name)



df <- as_tibble(read_delim(r"(data_frames\deseq\Biological Process\GOEA_upregulated_extreme_deseq.txt)"))

df <- df %>%
    transmute(
        name = `GO biological process complete`,
        fold_enrichment = ifelse(`upload_1 (over/under)` == "+", as.double(`upload_1 (fold Enrichment)`), -1 * (as.double(`upload_1 (fold Enrichment)`))),
        FDR = as.double(`upload_1 (FDR)`),
        Pval = `upload_1 (raw P-value)`
    )


p <- df %>%
    filter(FDR <= 0.01) %>%
    ggplot(df,
            mapping = aes(
                x = name,
                y = fold_enrichment
                )) +
        geom_col(
            aes(fill = fold_enrichment)
            ) +
        scale_fill_gradient(
            low = "#5e8319",
            high = "#f32121",
            space = "Lab",
            #na.value = "grey50",
            guide = "colourbar",
            aesthetics = "fill"
            ) +
        coord_flip() +
        theme(legend.position = "none") +
        labs(title = "Extreme Upregulated (Bio Process)")
p

png("Extreme_Upregulated_BioProcess.png")
print(last_plot())
dev.off()
