#AUROGENE - GSEA
    #C2 Collection (report for FB23-2)
    df <- read_tsv(r"(GSEA/Aurogene/FB23vsControl_Aurogene_C2.Gsea.1659172720632/gsea_report_for_FB23-2_1659172720632.tsv)")

    #Hallmark Collection (report for FB23-2)
    df <- read_tsv(r"(GSEA/Aurogene/FB23vsControl_Aurogene_HALLMARK/gsea_report_for_FB23-2_1659167997837.tsv)")

    #Hallmark Collection (report for Control)
    df <- read_tsv(r"(GSEA/Aurogene/FB23vsControl_Aurogene_HALLMARK/gsea_report_for_Control_1659167997837.tsv)")

#FPKMlibrary - GSEA
    #C2 Collection (report for FB23-2)
    df <- read_tsv(r"(GSEA/fpkm/FB23vsControl_FPKM_C2.Gsea.1659176504362/gsea_report_for_FB23-2_1659176504362.tsv)")

    #Hallmark Collection (report for FB23-2)
    df <- read_tsv(r"(GSEA/fpkm/FB23vsControl_FPKM_HALLMARK.Gsea.1659182157144/gsea_report_for_FB23-2_1659182157144.tsv)")


#123 - GO
#Gene Ontology disclaimer:
    # nelle cartelle GO ci sono diverse
    # analisi (upregulated, downregulated, completa).
    # significa che io ho messo nella query i geni
    # corrispondenti: nell'analisi "upregulated" GO ha calcolato
    # i pathway a partire SOLO dai geni upregulated ecc.
    # nell'analisi completa ci sono tutti i geni enriched.

    ##   in queste directory ho messo solo le analisi complete

    #Biological Process completo
    df <-  read_tsv(r"(GO/123/Gene_Ontology_123_biological_function/analisi_completa_123_BP.txt)")

    #Cellular Component completo 
    df <-  read_tsv(r"(GO/123/Gene_Ontology_123_cellular_component/analisi_completa_123_cellular_component.txt)")

    #Molecular Funx completo
    df <-  read_tsv(r"(GO/123/Gene_Ontology_123_molecular_function/analisi_completa_123_MF.txt)") # nolint

#deSeq - GO

    # qua ci sono solo up e downregulated
    # non mi ricordo più perchè
    #Biological Process downregulated
    df <- read_tsv(r"(GO/deseq/GO_deseq_BP/GO_deseq_downregulated_BP.txt)")

    #Biological Process upregulated
    df <- read_tsv(r"(GO/deseq/GO_deseq_BP/GO_deseq_upregulated_BP.txt)")


#### SCRIPT

#aggiungi variabile "significant - nonsignificant"
df <- df %>%
    select(c(NAME, SIZE, NES, `NOM p-val`, `FDR q-val`)) %>%
        arrange(NES) %>%
            mutate(
                adjPvalue = ifelse(
                df$`NOM p-val` <= 0.05,
                "significant",
                "non-significant"),
                ) %>%
                    view()

#filtra i pvalue significativi
df <- filter(df, `NOM p-val` <= 0.05)

#barplot con Enrichment Score; specifica la significatività del risultato
ggplot(df, aes(NAME, abs(NES), fill = NES)) +
    geom_col() +
    #scale_fill_manual(
    #    values = c(
    #        "significant" = "red",
    #        "non-significant" = "gray")) 
    scale_fill_gradient(
        low = "#5e8319",
        high = "#f32121",
        space = "Lab",
        #na.value = "grey50",
        guide = "colourbar",
        aesthetics = "fill"
        ) +

    coord_flip() +
    labs(
        x = "Pathway",
        y = "Normalized Enrichment Score",
        title = "C2 pathways Enrichment Score from GSEA")

#rinomina e filtra i dati usciti da GO
df <- df %>%
    rename(
        fold_enrichment = `upload_1 (fold Enrichment)`,
     name = `GO biological process complete`
        ) %>%
    arrange(`fold_enrichment`) %>%
    filter(
        `upload_1 (raw P-value)` <= 0.05,
        #`upload_1 (FDR)` <= 0.01
        ) %>%
    select(fold_enrichment, name)


#barplot con Fold Enrichment e gradiente di colore
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
    labs(title = "Enriched Pathways")

#questo funziona


df <- read.table(file = "clipboard", 
    sep = "\t", header=TRUE)
df1 <- read.table(file = "clipboard", 
    sep = "\t", header=TRUE)
df <- rbind(df, df1) %>%
    as_tibble()


write_tsv(df, "deseq2_data.tsv")

df <- read_tsv(r"(deseq2_data.tsv)")

df <- df %>%
    mutate(P.value = as.numeric(gsub(",", ".", gsub("\\.", "", P.value))))


flt <- df %>%
    select(Gene.ID,Gene.name,logFC,P.value,FDR) %>%
    filter(`P.value` <= 0.05) %>%
    filter(FDR <= 0.01) %>%
    view()

arrange(flt, desc(FDR))

flt <- flt %>%
    mutate(regulation = ifelse(logFC > 0, "upregulated", "downregulated"))



up <- flt %>% filter(regulation == "upregulated")
down <- flt %>% filter(regulation == "downregulated")


write_clip(down$Gene.name)

df %>%
    ggplot(aes(x = Gene.ID, y = logFC)) +
    geom_col(fill = logFC) +
    scale_fill_gradient(
    low = "#5e8319",
    high = "#f32121",
    space = "Lab",
    na.value = "grey50",
    guide = "colourbar",
    aesthetics = "fill"
    ) +
    coord_flip() +
    theme(legend.position="none") +
    labs(title = "Fold Change")
