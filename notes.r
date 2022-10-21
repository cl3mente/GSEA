library(tidyverse)
library(readxl)

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

#background di tutti i geni
    #devo prendere solo i geni "accesi" in almeno uno dei due fenotipi
    #per trovare i geni accesi devo fare una media di FPKM, che deve essere > 1
    #ricorda che FPKM è la quantità di RNA trovato nella cellula.
bck <- df %>%
    filter(average of control FPKM > 1 | average of treated FPKM > 1)



df <- rbind(dfup, dfdown)
#solo geni significativi
df %>%
    filter(FDR <= 0.01) %>%
    write_tsv(r"(geni_significativi.tsv)")

df %>%
    filter(logFC < 0) %>%
    write_tsv(r"(data_frames\deseq\geni_downregulated.tsv)")




#volcano plot
df %>%
    transmute(
    significant = ifelse(`FDR` <= 0.01, "yes", "no")
    ) %>%

    ggplot(aes(
    x = `logFC`,
    y = `logPval`,
    color = `significant`,
    alpha = .2)) +
    geom_point()
