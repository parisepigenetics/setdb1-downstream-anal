# Load packages.
library(magrittr)
library(DOSE)
library(clusterProfiler)
library(org.Mm.eg.db)
library(msigdbr)
library(ReactomePA)
library(enrichplot)
library(RColorBrewer)
library(gplots)
library(factoextra)
library(cluster)
library(mclust)

### Load/preprocess data ----------------------------------
dKK <- read.table("deSeq2_DiffAnal_WT_KK.tab", header = TRUE)
dCA <- read.table("deSeq2_DiffAnal_WT_CA.tab", header = TRUE)

# Keep only the rows with complete data.
dKK <- dKK[complete.cases(dKK),]
dCA <- dCA[complete.cases(dCA),]

# Keep only the significantly differentially expressed.
dKKs <- dKK[dKK$padj <= 0.05,]
dCAs <- dCA[dCA$padj <= 0.05,]

# Keep only the ones with logFC more (or less) than 1.
dKKsf1 <- dKKs[abs(dKKs$log2FoldChange) >= 1,]
dCAsf1 <- dCAs[abs(dCAs$log2FoldChange) >= 1,]
# Keep only the ones with logFC more (or less) than 0.5.
dKKsf05 <- dKKs[abs(dKKs$log2FoldChange) >= 0.5,]
dCAsf05 <- dCAs[abs(dCAs$log2FoldChange) >= 0.5,]

# Create name/logFC objects keep only the significant genes
geneListCA <- dCAs[["log2FoldChange"]]
geneListKK <- dKKs[["log2FoldChange"]]
namesCA <- as.character(dCAs[["Gene_ID"]])
namesKK <- as.character(dKKs[["Gene_ID"]])
names(geneListCA) <- namesCA
names(geneListKK) <- namesKK
geneListCA1 <- dCAsf1[["log2FoldChange"]]
geneListKK1 <- dKKsf1[["log2FoldChange"]]
namesCA1 <- as.character(dCAsf1[["Gene_ID"]])
namesKK1 <- as.character(dKKsf1[["Gene_ID"]])
names(geneListCA1) <- namesCA1
names(geneListKK1) <- namesKK1
length(namesCA1)
length(namesKK1)
geneListCA05 <- dCAsf05[["log2FoldChange"]]
geneListKK05 <- dKKsf05[["log2FoldChange"]]
namesCA05 <- as.character(dCAsf05[["Gene_ID"]])
namesKK05 <- as.character(dKKsf05[["Gene_ID"]])
names(geneListCA05) <- namesCA05
names(geneListKK05) <- namesKK05
length(namesCA05)
length(namesKK05)

# Write files to the disk
write(namesCA1, file = "deSeq2_DiffAnal_WT_CA_p05_FC1_geneNames.txt")
write(namesKK1, file = "deSeq2_DiffAnal_WT_KK_p05_FC1_geneNames.txt")
write(namesCA05, file = "deSeq2_DiffAnal_WT_CA_p05_FC05_geneNames.txt")
write(namesKK05, file = "deSeq2_DiffAnal_WT_KK_p05_FC05_geneNames.txt")



### Enrichments -------------------------------------------
# Get wikipathID names.
# Add column with the WikiGene ID.
convt <- read.table("wikiID_geneID.txt", sep = "\t")  # tabledownloaded from ENSEMBL BioMart
colnames(convt) <- c("WikiGeneID", "GeneID")
namesCAwiki1 <- convt[convt$GeneID %in% namesCA1,]$WikiGeneID
namesKKwiki1 <- convt[convt$GeneID %in% namesKK1,]$WikiGeneID
namesCAwiki05 <- convt[convt$GeneID %in% namesCA05,]$WikiGeneID
namesKKwiki05 <- convt[convt$GeneID %in% namesKK05,]$WikiGeneID
#namesCAwiki <- as.character(namesCAwiki[complete.cases(namesCAwiki)])
#namesKKwiki <- as.character(namesKKwiki[complete.cases(namesKKwiki)])
length(namesCAwiki1)
length(namesKKwiki1)
length(namesCAwiki05)
length(namesKKwiki05)

# Prepare the reference pathways and the TERM2 objects.
wp2genes <- read.gmt("../wikiPatways_Mm/wikipathways-20190110-gmt-Mus_musculus.gmt")
wp2genes <- wp2genes %>% tidyr::separate("ont", c("name","version","wpid","org"), "%")
wpid2gene <- wp2genes %>% dplyr::select("wpid", "gene") #TERM2GENE
wpid2name <- wp2genes %>% dplyr::select("wpid", "name") #TERM2NAME


## Perform WikiPaths enrichment analysis.
ewpCA1 <- enricher(namesCAwiki1, minGSSize = 40, TERM2GENE = wpid2gene, TERM2NAME = wpid2name)
barplot(ewpCA1, showCategory = 30, title = "WikiPaths enrichment CA1")
ewpKK1 <- enricher(namesKKwiki1, minGSSize = 20, TERM2GENE = wpid2gene, TERM2NAME = wpid2name)
barplot(ewpKK1, showCategory = 30,  title = "WikiPaths enrichment KK1")
ewpCA05 <- enricher(namesCAwiki05, minGSSize = 50, TERM2GENE = wpid2gene, TERM2NAME = wpid2name)
barplot(ewpCA05, showCategory = 30,  title = "WikiPaths enrichment CA05")
ewpKK05 <- enricher(namesKKwiki05, minGSSize = 10, TERM2GENE = wpid2gene, TERM2NAME = wpid2name)
barplot(ewpKK05, showCategory = 30,  title = "WikiPaths enrichment KK05")

# Prepare two new gene lists without filtering for logFC.
namesCAs <- as.character(dCAs[["Gene_ID"]])
namesCAwikiS <- convt[convt$GeneID %in% namesCAs,]$WikiGeneID
namesKKs <- as.character(dKKs[["Gene_ID"]])
namesKKwikiS <- convt[convt$GeneID %in% namesKKs,]$WikiGeneID
length(namesCAwikiS)
length(namesKKwikiS)


## Perform WikiPaths enrichment experiment WITHOUT logFC.
ewpCA0 <- enricher(namesCAwikiS, minGSSize = 50, TERM2GENE = wpid2gene, TERM2NAME = wpid2name)
barplot(ewpCA0, showCategory = 30,  title = "WikiPaths enrichment ALL DEGs")
ewpKK0 <- enricher(namesKKwikiS, minGSSize = 50, TERM2GENE = wpid2gene, TERM2NAME = wpid2name)
barplot(ewpKK0, showCategory = 30,  title = "WikiPaths enrichment ALL KK DEGs")


# Prepare the geneLists (ordered gene list)
# These are the gene names that we have Wiki IDs.
nn <- unique(convt[convt$WikiGene.ID %in% namesCAwikiS,]$GeneID)
mm <- unique(convt[convt$WikiGene.ID %in% namesKKwikiS,]$GeneID)

glCA <- geneListCA[names(geneListCA) %in% levels(nn)]
glKK <- geneListKK[names(geneListKK) %in% levels(mm)]

length(unique(convt[convt$GeneID %in% names(glCA),]$WikiGeneID))
length(names(glCA))
length(convt[convt$WikiGeneID %in% wpid2gene$gene,]$GeneID)
subset(convt, convt$WikiGeneID == wpid2gene[1,]$gene)[1,]$GeneID  # Just to check

# Conversion from wikipath to geneID
wpid2geneID <- data.frame();
for (i in 1:nrow(wpid2gene)) {
    wp <- wpid2gene[i,]$wpid
    gn <- wpid2gene[i,]$gene
    gid <- as.character(subset(convt, convt$WikiGeneID == gn)[1,]$GeneID)
    wpid2geneID <- rbind(wpid2geneID, data.frame(wp, gid))
}
wpid2geneID <- wpid2geneID[complete.cases(wpid2geneID),]
dim(wpid2geneID)
head(wpid2geneID)

# Sort the gene lists.
geneListCA1 <- sort(geneListCA1, decreasing = TRUE)
geneListKK1 <- sort(geneListKK1, decreasing = TRUE)
geneListCA05 <- sort(geneListCA05, decreasing = TRUE)
geneListKK05 <- sort(geneListKK05, decreasing = TRUE)

## Perform Gene Set Enrichment of WikiPaths.
ewpgsCA1 <- GSEA(geneListCA1, TERM2GENE = wpid2geneID, TERM2NAME = wpid2name, verbose = FALSE)
barplot(ewpgsCA1, title = "WikiPaths GSEA CA1")
ewpgsKK1 <- GSEA(geneListKK1, TERM2GENE = wpid2geneID, TERM2NAME = wpid2name, verbose = FALSE)
barplot(ewpgsKK1, title = "WikiPaths GSEA KK1")
ewpgsCA05 <- GSEA(geneListCA05, TERM2GENE = wpid2geneID, TERM2NAME = wpid2name, verbose = FALSE)
barplot(ewpgsCA05, showCategory = 30, title = "WikiPaths GSEA CA05")
ewpgsKK05 <- GSEA(geneListKK05, TERM2GENE = wpid2geneID, TERM2NAME = wpid2name, verbose = FALSE)
barplot(ewpgsKK05, showCategory = 30, title = "WikiPaths GSEA KK05")


# MSIGDB enrichment
m_df <- msigdbr(species = "Mus musculus")
m_df$gs_id <- m_df$gene_symbol # BIG TRICK TO SWAP THE COLUMN NAMES!!!!!!


## Perform GSIGDB gene enrichment.
esigCA1 <- enricher(namesCA1, minGSSize = 50, TERM2GENE = m_df)
barplot(esigCA1, showCategory = 50, title = "MsiGDB enrichment CA1")
esigKK1 <- enricher(namesKK1, minGSSize = 50, TERM2GENE = m_df)
barplot(esigKK1, showCategory = 50, title = "MsiGDB enrichment KK1")
esigCA05 <- enricher(namesCA05, minGSSize = 50, TERM2GENE = m_df)
barplot(esigCA05, showCategory = 50, title = "MsiGDB enrichment CA05")
esigKK05 <- enricher(namesKK05, minGSSize = 50, TERM2GENE = m_df)
barplot(esigKK05, showCategory = 50, title = "MsiGDB enrichment KK05")


## Perform GSIGDB Gene Set Enrichment.
esigsCA1 <- GSEA(geneListCA1, minGSSize = 40, TERM2GENE = m_df)
dotplot(esigsCA1, showCategory = 50, title = "MsiGDB GSEA CA1")
esigsKK1 <- GSEA(geneListKK1, minGSSize = 20, TERM2GENE = m_df)
dotplot(esigsKK1, showCategory = 50, title = "MsiGDB GSEA KK1")
esigsCA05 <- GSEA(geneListCA05, minGSSize = 10, TERM2GENE = m_df)
dotplot(esigsCA05, showCategory = 50, title = "MsiGDB GSEA CA05")
esigsKK05 <- GSEA(geneListKK05, minGSSize = 20, TERM2GENE = m_df)
dotplot(esigsKK05, showCategory = 50, title = "MsiGDB GSEA KK05")

# Transform the common gene names to ENTREZIDs
nidCA1 <- bitr(namesCA1, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
nidKK1 <- bitr(namesKK1, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
nidCA05 <- bitr(namesCA05, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
nidKK05 <- bitr(namesKK05, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)

## Perform GO group enrichment analysis.
ggoCA1_MF2 <- groupGO(gene = nidCA1$ENTREZID, OrgDb = org.Mm.eg.db, ont = "MF", level = 2, readable = TRUE)
barplot(ggoCA1_MF2, showCategory = 40,  title = "GroupGO CA1 MF2")
ggoKK1_MF2 <- groupGO(gene = nidKK1$ENTREZID, OrgDb = org.Mm.eg.db, ont = "MF", level = 2, readable = TRUE)
barplot(ggoKK1_MF2, showCategory = 40,  title = "GroupGO KK1 MF2")
ggoCA05_MF2 <- groupGO(gene = nidCA05$ENTREZID, OrgDb = org.Mm.eg.db, ont = "MF", level = 2, readable = TRUE)
barplot(ggoCA05_MF2, showCategory = 40, title = "GroupGO CA05 MF2")
ggoKK05_MF2 <- groupGO(gene = nidKK05$ENTREZID, OrgDb = org.Mm.eg.db, ont = "MF", level = 2, readable = TRUE)
barplot(ggoKK05_MF2, showCategory = 40, title = "GroupGO KK05 MF2")
# Here one needs to play a lot with the level.


## Perform the enrichment in GO Molecular Functions analysis.
egoCA1_MF <- enrichGO(gene = nidCA1$ENTREZID, OrgDb = org.Mm.eg.db, minGSSize = 20, ont = "MF", pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)
dotplot(egoCA1_MF, title = "EnrichGO CA1 MF")
egoKK1_MF <- enrichGO(gene = nidKK1$ENTREZID, OrgDb = org.Mm.eg.db, minGSSize = 50, ont = "MF", pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)
dotplot(egoKK1_MF, title = "EnrichGO KK1 MF")
egoCA05_MF <- enrichGO(gene = nidCA05$ENTREZID, OrgDb = org.Mm.eg.db, minGSSize = 30, ont = "MF", pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)
dotplot(egoCA05_MF, title = "EnrichGO CA05 MF")
egoKK05_MF <- enrichGO(gene = nidKK05$ENTREZID, OrgDb = org.Mm.eg.db, minGSSize = 20, ont = "MF", pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)
dotplot(egoKK05_MF, title = "EnrichGO KK05 MF")


## Perform enrichment in GO Biological Processes analysis.
egoCA1_BP <- enrichGO(gene = nidCA1$ENTREZID, OrgDb = org.Mm.eg.db, minGSSize = 20, ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)
dotplot(egoCA1_BP, title = "EnrichGO CA1 BP")
egoKK1_BP <- enrichGO(gene = nidKK1$ENTREZID, OrgDb = org.Mm.eg.db, minGSSize = 20, ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)
dotplot(egoKK1_BP, title = "EnrichGO KK1 BP")
egoCA05_BP <- enrichGO(gene = nidCA05$ENTREZID, OrgDb = org.Mm.eg.db, minGSSize = 20, ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)
dotplot(egoCA05_BP, title =  "EnrichGO CA05 BP")
egoKK05_BP <- enrichGO(gene = nidKK05$ENTREZID, OrgDb = org.Mm.eg.db, minGSSize = 20, ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)
dotplot(egoKK05_BP, title = "EnrichGO KK05 BP")  # Interesting enrichments for all.


## Perform gene set enrichment in GO Molecular Functions analysis.
egogsCA1_MF <- gseGO(geneList = geneListCA1, OrgDb = org.Mm.eg.db, ont = "MF", nPerm = 1000, minGSSize = 50, maxGSSize = 500, pvalueCutoff = 0.05, verbose = FALSE, keyType = "SYMBOL")
dotplot(egogsCA1_MF, title = "GSEA GO MF CA1")
egogsKK1_MF <- gseGO(geneList = geneListKK1, OrgDb = org.Mm.eg.db, ont = "MF", nPerm = 1000, minGSSize = 50, maxGSSize = 500, pvalueCutoff = 0.05, verbose = FALSE, keyType = "SYMBOL")
dotplot(egogsKK1_MF, title = "GSEA GO MF KK1")
egogsCA05_MF <- gseGO(geneList = geneListCA05, OrgDb = org.Mm.eg.db, ont = "MF", nPerm = 1000, minGSSize = 50, maxGSSize = 500, pvalueCutoff = 0.05, verbose = FALSE, keyType = "SYMBOL")
dotplot(egogsKK1_MF, title = "GSEA GO MF CA05")
egogsKK05_MF <- gseGO(geneList = geneListKK05, OrgDb = org.Mm.eg.db, ont = "MF", nPerm = 1000, minGSSize = 50, maxGSSize = 500, pvalueCutoff = 0.05, verbose = FALSE, keyType = "SYMBOL")
dotplot(egogsKK05_MF, title = "GSEA GO MF KK05")
# NONE enriched


## Perform the gene set enrichment in GO biological processes analysis.
egogsCA1_BP <- gseGO(geneList = geneListCA1, OrgDb = org.Mm.eg.db, ont = "BP", nPerm = 1000, minGSSize = 50, maxGSSize = 500, pvalueCutoff = 0.05, verbose = FALSE, keyType = "SYMBOL")
dotplot(egogsCA1_BP, title = "GSEA GO BP CA1")
egogsKK1_BP <- gseGO(geneList = geneListKK1, OrgDb = org.Mm.eg.db, ont = "BP", nPerm = 1000, minGSSize = 50, maxGSSize = 500, pvalueCutoff = 0.05, verbose = FALSE, keyType = "SYMBOL")
dotplot(egogsKK1_BP, title = "GSEA GO BP KK1")
egogsCA05_BP <- gseGO(geneList = geneListCA05, OrgDb = org.Mm.eg.db, ont = "BP", nPerm = 1000, minGSSize = 50, maxGSSize = 500, pvalueCutoff = 0.05, verbose = FALSE, keyType = "SYMBOL")
dotplot(egogsCA05_BP, title = "GSEA GO BP CA05")
egogsKK05_BP <- gseGO(geneList = geneListKK05, OrgDb = org.Mm.eg.db, ont = "BP", nPerm = 1000, minGSSize = 50, maxGSSize = 500, pvalueCutoff = 0.05, verbose = FALSE, keyType = "SYMBOL")
dotplot(egogsKK05_BP, title = "GSEA GO BP KK05")
# No enrichment found.


## Perform KEEG pathway enrichment.
ekeCA1 <- enrichKEGG(gene = nidCA1$ENTREZID, organism = "mmu", pvalueCutoff = 0.05)
barplot(ekeCA1, title = "KEGG enrichment CA1")
ekeKK1 <- enrichKEGG(gene = nidKK1$ENTREZID, organism = "mmu", pvalueCutoff = 0.05)
barplot(ekeKK1, title = "KEGG enrichment KK1")
ekeCA05 <- enrichKEGG(gene = nidCA05$ENTREZID, organism = "mmu", pvalueCutoff = 0.05)
barplot(ekeCA05, showCategory = 60, title = "KEGG enrichment CA05") ##
ekeKK05 <- enrichKEGG(gene = nidKK05$ENTREZID, organism = "mmu", pvalueCutoff = 0.05)
barplot(ekeKK05, title = "KEGG enrichment KK05")

ekmCA1 <- enrichMKEGG(gene = nidCA1$ENTREZID, organism = "mmu")
barplot(ekmCA1, title = "KEGG modules enrichment CA1")
ekmKK1 <- enrichMKEGG(gene = nidKK1$ENTREZID, organism = "mmu")
barplot(ekmKK1, showCategory = 30, title = "KEGG enrichment CA05")
ekmCA05 <- enrichMKEGG(gene = nidCA05$ENTREZID, organism = "mmu")
barplot(ekmCA05, title = "KEGG modules enrichment KK1")
ekmKK05 <- enrichMKEGG(gene = nidKK05$ENTREZID, organism = "mmu")
barplot(ekmKK05, title = "KEGG enrichment KK05")

epaCA1 <- enrichPathway(gene = nidCA1$ENTREZID, organism = "mouse", pvalueCutoff = 0.05)
barplot(epaCA1, showCategory = 30, title = "Pathways enrichment CA1")
epaKK1 <- enrichPathway(gene = nidKK1$ENTREZID, organism = "mouse", pvalueCutoff = 0.05)
barplot(epaKK1, title = "Pathways enrichment KK1")
epaCA05 <- enrichPathway(gene = nidCA05$ENTREZID, organism = "mouse", pvalueCutoff = 0.05)
barplot(epaCA05, title = "Pathways enrichment CA05")
epaKK05 <- enrichPathway(gene = nidKK05$ENTREZID, organism = "mouse", pvalueCutoff = 0.05)
barplot(epaKK05, showCategory = 20, title = "Pathways enrichment KK05")



### Visualisations ----------------------------------------
edoCA1 <- enrichDGN(nidCA1$ENTREZID)
edoKK1 <- enrichDGN(nidKK1$ENTREZID)
# We can possibly do that if we transform to human orthologous genes.

## BArplots - Dotplots.
# for Jupyter use the following to produce better figures
#options(repr.plot.width=10,repr.plot.height=14) )
barplot(ggoCA1_MF2, drop = TRUE, showCategory = 200, title = "")
barplot(ggoKK1_MF2, drop = TRUE, showCategory = 200, title = "")
dotplot(egoCA05_MF, showCategory = 30) + ggtitle("Dotplot for MF enrichment for CA")
dotplot(egoCA1_BP, showCategory = 50) + ggtitle("Dotplot for BP enrichment for CA")
dotplot(egoKK05_MF, showCategory = 50) + ggtitle("Dotplot for MF enrichment for KK")
dotplot(egoKK05_BP, showCategory = 50) + ggtitle("Dotplot for BP enrichment for KK")


## Enrichment Map GO plots
# BF
emapplot(egoCA1_BP) + ggtitle("EMAplot BP CA1")
emapplot(egoKK1_BP) + ggtitle("EMAplot BP KK1")
emapplot(egoCA05_BP) + ggtitle("EMAplot BP CA05")
emapplot(egoKK05_BP) + ggtitle("EMAplot BP KK05")
# MF
emapplot(egoCA1_MF) + ggtitle("EMAplot MF CA1")
emapplot(egoKK1_MF) + ggtitle("EMAplot MF KK1")
emapplot(egoCA05_MF) + ggtitle("EMAplot MF CA05")
emapplot(egoKK05_MF) + ggtitle("EMAplot MF KK05")


## Category Network (CNET) plots (perhaps the most usefull!)
cnetplot(egoCA1_MF, foldChange = geneListCA1, colorEdge = TRUE) + ggtitle("CNETplot MF CA1")
cnetplot(egoKK1_MF, foldChange = geneListKK1, colorEdge = TRUE) + ggtitle("CNETplot MF KK1")
cnetplot(egoCA05_MF, foldChange = geneListCA05, colorEdge = TRUE) + ggtitle("CNETplot MF CA05")
cnetplot(egoKK05_MF, foldChange = geneListKK05, colorEdge = TRUE) + ggtitle("CNETplot MF KK05")
cnetplot(egoCA1_BP, foldChange = geneListCA1, colorEdge = TRUE) + ggtitle("CNETplot BP CA1")
cnetplot(egoKK1_BP, foldChange = geneListKK1, colorEdge = TRUE) + ggtitle("CNETplot BP KK1")
cnetplot(egoCA05_BP, foldChange = geneListCA05, colorEdge = TRUE) + ggtitle("CNETplot BP CA05")
cnetplot(egoKK05_BP, foldChange = geneListKK05, colorEdge = TRUE) + ggtitle("CNETplot BP KK05")


## GOplots (extended EMA plots, perhaps confusing)
goplot(egoCA05_MF, showCategory = 30, geom = "text", alph = 0.50) + ggtitle("GOAplot MF CA05")
goplot(egoKK05_MF, showCategory = 10, geom = "text", alph = 0.50) + ggtitle("EMAplot MF KK05")
goplot(egoCA05_BP, showCategory = 10, geom = "text", alph = 0.50) + ggtitle("GOplot BP CA05")
goplot(egoKK05_BP, showCategory = 10, geom = "text", alph = 0.50) + ggtitle("EMAplot BP KK05")


## Perform GO gene set enrichment analysis.
# MF
edoCA1_MF <- gseGO(geneListCA1, nPerm = 10000, OrgDb = org.Mm.eg.db, ont = "MF", keyType = "SYMBOL")
dotplot(edoCA1_MF, title = "GSEA GO MF CA1") #
edoKK1_MF <- gseGO(geneListKK1, nPerm = 10000, OrgDb = org.Mm.eg.db, ont = "MF", keyType = "SYMBOL")
dotplot(edoKK1_MF, title = "GSEA GO MF KK1")
edoCA05_MF <- gseGO(geneListCA05, nPerm = 10000, OrgDb = org.Mm.eg.db, ont = "MF", keyType = "SYMBOL")
dotplot(edoCA05_MF,  title = "GSEA GO MF CA05")
edoKK05_MF <- gseGO(geneListKK05, nPerm = 10000, OrgDb = org.Mm.eg.db, ont = "MF", keyType = "SYMBOL")
dotplot(edoKK05_MF,  title = "GSEA GO MF KK05")
# BP
edoCA1_BP <- gseGO(geneListCA1, nPerm = 10000, OrgDb = org.Mm.eg.db, ont = "BP", keyType = "SYMBOL")
dotplot(edoCA1_BP, title = "GSEA GO BP CA1") #
edoKK1_BP <- gseGO(geneListKK1, nPerm = 10000, OrgDb = org.Mm.eg.db, ont = "BP", keyType = "SYMBOL")
dotplot(edoKK1_BP, title = "GSEA GO BP KK1")
edoCA05_BP <- gseGO(geneListCA05, nPerm = 10000, OrgDb = org.Mm.eg.db, ont = "BP", keyType = "SYMBOL")
dotplot(edoCA05_BP, title = "GSEA GO BP CA05")
edoKK05_BP <- gseGO(geneListKK05, nPerm = 10000, OrgDb = org.Mm.eg.db, ont = "BP", keyType = "SYMBOL")
dotplot(edoKK05_BP, title = "GSEA GO BP KK05")


# Heat plots, grouping and gene expression!
heatplot(edoCA1_BP, foldChange = geneListCA1)
heatplot(edoKK05_BP, foldChange = geneListKK1)

# Ridge plot (funcy)
ridgeplot(edoCA1_BP)

# GSEA plots, USEFULL more difficult to interpret.
gseaplot2(edoCA1_BP, geneSetID = 1, title = edoCA1_BP$Description[1])
gseaplot2(edoCA1_BP, geneSetID = 2, title = edoCA1_BP$Description[2])
gseaplot2(edoCA1_BP, geneSetID = 1:5, title = "GSES of the top 5")  # This is usefull as we can superimpose many different enrichments!



### Clustering --------------------------------------------
# Collect data.
ddCA1 <- read.table("tpm_WT_CA_KK_NLSB_degsCA_p05_fc1.tab", header = TRUE)
ddKK1 <- read.table("tpm_WT_CA_KK_NLSB_degsKK_p05_fc1.tab", header = TRUE)
ddCA05 <- read.table("tpm_WT_CA_KK_NLSB_degsCA_p05_fc05.tab", header = TRUE)
ddKK05 <- read.table("tpm_WT_CA_KK_NLSB_degsKK_p05_fc05.tab", header = TRUE)
#Clean the zeros
ddCA1 <- ddCA1[apply(ddCA1, 1, function(row) all(row != 0 )),]
ddKK1 <- ddKK1[apply(ddKK1, 1, function(row) all(row != 0 )),]
ddCA05 <- ddCA05[apply(ddCA05, 1, function(row) all(row != 0 )),]
ddKK05 <- ddKK05[apply(ddKK05, 1, function(row) all(row != 0 )),]
summary(ddCA05)

# Get the logs. (natural)
dCA1 <- log(ddCA1[,c(1,2,3,4,5,6)])
dKK1 <- log(ddKK1[,c(1,2,3,7,8,9)])
dCA05 <- log(ddCA05[,c(1,2,3,4,5,6)])
dKK05 <- log(ddKK05[,c(1,2,3,7,8,9)])
colnames(dCA1) <- c("WT1", "WT2", "WT3", "CA1", "CA2", "CA3")
colnames(dCA05) <- c("WT1", "WT2", "WT3", "CA1", "CA2", "CA3")
colnames(dKK1) <- c("WT1", "WT2", "WT3", "KK1", "KK2", "KK3")
colnames(dKK05) <- c("WT1", "WT2", "WT3", "KK1", "KK2", "KK3")

### For the CA DEGs
## Hierarhical clustering of CA1
my_palette <- brewer.pal(n = 11, name = "RdYlGn")
# Do the clustering.
clust_CA1 <- hclust(dist(dCA1), method = "single")
# define clusters (hard thresold)
CA1_Clusts <- cutree(clust_CA1, h = max(clust_CA1$height/1.5))
# Colour vector for clusters side bar.
myClustCols <- rainbow(length(unique(CA1_Clusts)))
myClusters_CA1 <- myClustCols[CA1_Clusts]
heatmap.2(as.matrix(dCA1), main = "DEGs CA_1 hierarchical clustering", Rowv = as.dendrogram(clust_CA1), Colv = FALSE, dendrogram = "row", scale = "row", col = my_palette, cexCol = 2, cexRow = 0.4, key.title = NA, keysize = 0.8, key.xlab = NA, ylab = "Genes", RowSideColors = myClusters_CA1)

m## Hierarhical clustering of CA05
clust_CA05 <- hclust(dist(dCA05), method = "single")
# define clusters (hard thresold)
CA05_Clusts <- cutree(clust_CA05, h = max(clust_CA05$height/2))
#CA05_Clusts <- cutree(clust_CA05, k = 3)
# Colour vector for clusters side bar.
myClustCols <- rainbow(length(unique(CA05_Clusts)))
myClusters_CA05 <- myClustCols[CA05_Clusts]
heatmap.2(as.matrix(dCA05), main = "DEGs CA_05 hierarchical clustering", Rowv = as.dendrogram(clust_CA05), Colv = FALSE, dendrogram = "row", scale = "row", col = my_palette, cexCol = 2, cexRow = 0.4, key.title = NA, keysize = 0.8, key.xlab = NA, ylab = "Genes", RowSideColors = myClusters_CA05)


## Determine number of clusters by plotting sum of squares whitnin groups
wssCA <- (nrow(dCA1))*sum(apply(dCA1,2,var))
for (i in 2:15) wssCA[i] <- sum(kmeans(dCA1, centers = i, iter.max = 100, nstart = 250)$withinss)
my_palette <- brewer.pal(n = 11, name = "RdYlGn")
plot(1:15, wssCA, type = "b", xlab = "Number of Clusters", ylab = "Within groups sum of squares CA")

## Perform k-means clustering
#set.seed(1)  # Always set the same random seed.
kmCA1 <- kmeans(as.matrix(dCA1), 6)
kmCA05 <- kmeans(as.matrix(dCA05), 6)

# Append id and cluster
dfcallCA1 <- cbind(dCA1, id = seq(nrow(dCA1)), cluster = kmCA1$cluster)
# Add idsort, the id number ordered by cluster
dfcallCA1$idsort <- dfcallCA1$id[order(dfcallCA1$cluster)]
dfcallCA1$idsort <- order(dfcallCA1$idsort)
clusterColsCA1 <- as.character(sort(kmCA1$cluster))

## Plot k-means clustering.
heatmap(as.matrix(dCA1)[order(kmCA1$cluster),], Rowv = NA, col = my_palette, Colv = NA, cexCol = 1.5, cexRow = 0.4, RowSideColors = clusterColsCA1, ylab = "Genes", main = "k-means clustering of CA_1 TPMs")


## For the KK DEGs.
## Hierarchical clustering of KK1
clust_KK1 <- hclust(dist(dKK1), method = "single")
# define clusters (hard thresold)
KK1_Clusts <- cutree(clust_KK1, h = max(clust_KK1$height/1.3))
#kk1_Clusts <- cutree(clust_KK1, k = 3)
# Colour vector for clusters side bar.
myClustCols <- rainbow(length(unique(KK1_Clusts)))
myClusters_KK1 <- myClustCols[KK1_Clusts]
heatmap.2(as.matrix(dKK1), main = "DEGs KK_1 hierarchical clustering", Rowv = as.dendrogram(clust_KK1), Colv = FALSE, dendrogram = "row", scale = "row", col = my_palette, cexCol = 2, cexRow = 0.6, key.title = NA, keysize = 0.8, key.xlab = NA, ylab = "Genes", RowSideColors = myClusters_KK1)

## Hierarchical clustering of KK05
clust_KK05 <- hclust(dist(dKK05), method = "single")
# define clusters (hard thresold)
KK05_Clusts <- cutree(clust_KK05, h = max(clust_KK05$height/3))
#KK05_Clusts <- cutree(clust_KK05, k = 20)
# Colour vector for clusters side bar.
myClustCols <- rainbow(length(unique(KK05_Clusts)))
myClusters_KK05 <- myClustCols[KK05_Clusts]
heatmap.2(as.matrix(dKK05), main = "DEGs KK_05 hierarchical clustering", Rowv = as.dendrogram(clust_KK05), dendrogram = "row", Colv = FALSE, scale = "row", col = my_palette, cexCol = 2, cexRow = 0.6, key.title = NA, keysize = 0.8, key.xlab = NA, ylab = "Genes", RowSideColors = myClusters_KK05)

# determine number of cluster by plotting sum of squares whitnin groups
wssKK <- (nrow(dKK05))*sum(apply(dKK05,2,var))
for (i in 2:15) wssKK[i] <- sum(kmeans(dKK05, centers = i, iter.max = 100, nstart = 250)$withinss)
plot(1:15, wssKK, type = "b", xlab = "Number of Clusters", ylab = "Within groups sum of squares KK")

## Perform k-means clustering
#set.seed(1)  # Always set the same random seed.
kmKK1 <- kmeans(dKK1, 4, iter.max = 500, nstart = 200)
kmKK05 <- kmeans(dKK05, 4, iter.max = 500, nstart = 200)

# Append id and cluster
dfcallKK05 <- cbind(dKK05, id = seq(nrow(dKK05)), cluster = kmKK05$cluster)
# Add idsort, the id number ordered by cluster
dfcallKK05$idsort <- dfcallKK05$id[order(dfcallKK05$cluster)]
dfcallKK05$idsort <- order(dfcallKK05$idsort)
clusterColsKK05 <- as.character(sort(kmKK05$cluster))

# Plot k-means clustering.
heatmap(as.matrix(dKK05)[order(kmKK05$cluster),], Rowv = NA, Colv = NA, col = my_palette, cexCol = 1.5, cexRow = 0.7, RowSideColors = clusterColsKK05, ylab = "Genes", main = "k-means clustering of KK_05 TPMs")


# The PAM method.
pam_dKK05 <- pam(dKK05, 4)

## Perform some clustering tests and visulaisations.
km4_KK05.res <- kmeans(dKK05, 4, nstart = 1000)
km2_KK05.res <- kmeans(dKK05, 4, nstart = 1000)
kmCorr_KK05.res <- kmeans(cor(t(dKK05)), 2, nstart = 250)
clust_KK05 <- hclust(dist(dKK05), method = "single")

fviz_cluster(kmCorr_KK05.res, data = dKK05, ellipse.type = "convex") + theme_minimal()


# Lastly try also the MClust package.
mclust_KK05 <- Mclust(dKK05)
summary(mclust_KK05)

fviz_cluster(mclust_KK05, data = dKK05, ellipse.type = "convex") + theme_minimal()



## Clustering raw data of DEGs is not so informative.

## We proceed by clustering the up and dow regulated genes independently.

### CLUSTERING OF UP AND DOW REGULATED GENES IN CA1 and KK05.

# Prepare the data frames.
# Keep only the ones with logFC more (or less) than 1.
dCAsf1u <- dCAs[dCAs$log2FoldChange >= 1,]
dCAsf1d <- dCAs[dCAs$log2FoldChange <= -1,]
# Keep only the ones with logFC more (or less) than 0.5.
dKKsf05u <- dKKs[dKKs$log2FoldChange >= 0.5,]
dKKsf05d <- dKKs[dKKs$log2FoldChange <= -0.5,]

# Create name/logFC objects keep only the significant genes
geneListCA1u <- dCAsf1u[["log2FoldChange"]]
geneListCA1d <- dCAsf1d[["log2FoldChange"]]
namesCA1u <- as.character(dCAsf1u[["Gene_ID"]])
namesCA1d <- as.character(dCAsf1d[["Gene_ID"]])
names(geneListCA1u) <- namesCA1u
names(geneListCA1d) <- namesCA1d
length(namesCA1u)
length(namesCA1d)
geneListKK05u <- dKKsf05u[["log2FoldChange"]]
geneListKK05d <- dKKsf05d[["log2FoldChange"]]
namesKK05u <- as.character(dKKsf05u[["Gene_ID"]])
namesKK05d <- as.character(dKKsf05d[["Gene_ID"]])
names(geneListKK05u) <- namesKK05u
names(geneListKK05d) <- namesKK05d
length(namesKK05u)
length(namesKK05d)

## Select the data
ddCA1u <- ddCA05[namesCA1u,]
ddCA1d <- ddCA05[namesCA1d,]
ddKK05u <- ddKK05[namesKK05u,]
ddKK05d <- ddKK05[namesKK05d,]

#Clean zeros and NAs.
ddCA1u <- ddCA1u[apply(ddCA1u, 1, function(row) all(row != 0 )),]
ddCA1d <- ddCA1d[apply(ddCA1d, 1, function(row) all(row != 0 )),]
ddKK05u <- ddKK05u[apply(ddKK05u, 1, function(row) all(row != 0 )),]
ddKK05d <- ddKK05u[apply(ddKK05d, 1, function(row) all(row != 0 )),]

# Get the logs. (natural)
dCA1u <- log(ddCA1u[,c(1,2,3,4,5,6)])
dCA1d <- log(ddCA1d[,c(1,2,3,4,5,6)])
dKK05u <- log(ddKK05[,c(1,2,3,7,8,9)])
dKK05d <- log(ddKK05[,c(1,2,3,7,8,9)])
colnames(dCA1u) <- c("WT1", "WT2", "WT3", "CA1", "CA2", "CA3")
colnames(dCA1d) <- c("WT1", "WT2", "WT3", "CA1", "CA2", "CA3")
colnames(dKK05u) <- c("WT1", "WT2", "WT3", "KK1", "KK2", "KK3")
colnames(dKK05d) <- c("WT1", "WT2", "WT3", "KK1", "KK2", "KK3")

## Hierarhical clustering of CA1-up
my_palette <- brewer.pal(n = 11, name = "RdYlGn")
# Do the clustering.
clust_CA1u <- hclust(dist(dCA1u), method = "single")
# define clusters (hard thresold)
CA1_Clusts <- cutree(clust_CA1, h = max(clust_CA1$height/1.5))
# Colour vector for clusters side bar.
myClustCols <- rainbow(length(unique(CA1_Clusts)))
myClusters_CA1 <- myClustCols[CA1_Clusts]
heatmap.2(as.matrix(dCA1u), main = "DEGs CA_1 hierarchical clustering", Rowv = as.dendrogram(clust_CA1u), Colv = FALSE, dendrogram = "row", scale = "row", col = my_palette, cexCol = 2, cexRow = 0.4, key.title = NA, keysize = 0.8, key.xlab = NA, ylab = "Genes", RowSideColors = myClusters_CA1)
