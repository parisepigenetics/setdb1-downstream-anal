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
# We have preselecteed for logFC1 and significant genes.
dnls <- read.table("DESeq2_Genes_NLS_vs_WT.tab", header = TRUE)

# Create name/logFC objects keep only the significant genes
geneListNLS <- dnls[["log2FoldChange"]]
namesNLS <- as.character(rownames(dnls))
names(geneListNLS) <- namesNLS
geneListNLS <- sort(geneListNLS, decreasing = TRUE)

# Create also two separated lists one for UP one for DOWN regulated.
geneListNLSUP <- geneListNLS[geneListNLS > 0]
geneListNLSDOWN <- geneListNLS[geneListNLS < 0]
namesNLSUP <- names(geneListNLSUP)
namesNLSDOWN <- names(geneListNLSDOWN)


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
wp2genes <- read.gmt("../../wikiPatways_Mm/wikipathways-20190110-gmt-Mus_musculus.gmt")
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
barplot(ewpCA1, title = "WikiPaths GSEA CA1")
ewpgsKK1 <- GSEA(geneListKK1, TERM2GENE = wpid2geneID, TERM2NAME = wpid2name, verbose = FALSE)
barplot(ewpKK1, title = "WikiPaths GSEA KK1")
ewpgsCA05 <- GSEA(geneListCA05, TERM2GENE = wpid2geneID, TERM2NAME = wpid2name, verbose = FALSE)
barplot(ewpCA05, showCategory = 30, title = "WikiPaths GSEA CA05")
ewpgsKK05 <- GSEA(geneListKK05, TERM2GENE = wpid2geneID, TERM2NAME = wpid2name, verbose = FALSE)
barplot(ewpKK05, showCategory = 30, title = "WikiPaths GSEA KK05")


# MSIGDB enrichment
m_df <- msigdbr(species = "Mus musculus")
m_df$gs_id <- m_df$gene_symbol # BIG TRICK TO SWAP THE COLUMN NAMES!!!!!!


## Perform GSIGDB gene enrichment.
esigALL <- enricher(namesNLS, minGSSize = 50, TERM2GENE = m_df)
barplot(esigALL, showCategory = 50, title = "MsiGDB enrichment ALL")
esigUP <- enricher(namesNLSUP, minGSSize = 50, TERM2GENE = m_df)
barplot(esigUP, showCategory = 50, title = "MsiGDB enrichment UP")
esigDOWN <- enricher(namesNLSDOWN, minGSSize = 50, TERM2GENE = m_df)
barplot(esigDOWN, showCategory = 50, title = "MsiGDB enrichment DOWN")

## Perform GSIGDB Gene Set Enrichment.
esigsALL <- GSEA(geneListNLS, minGSSize = 40, TERM2GENE = m_df)
dotplot(esigsALL, showCategory = 50, title = "MsiGDB GSEA ALL")
esigsUP <- GSEA(geneListNLSUP, minGSSize = 20, TERM2GENE = m_df)
dotplot(esigsUP, showCategory = 50, title = "MsiGDB GSEA UP")
esigsDOWN <- GSEA(geneListNLSDOWN, minGSSize = 10, TERM2GENE = m_df)
dotplot(esigsDOWN, showCategory = 50, title = "MsiGDB GSEA DOWN")

# Transform the common gene names to ENTREZIDs
nidALL <- bitr(namesNLS, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
nidUP <- bitr(namesNLSUP, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
nidDOWN <- bitr(namesNLSDOWN, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)

## Perform GO group enrichment analysis.
ggoALL_MF3 <- groupGO(gene = nidALL$ENTREZID, OrgDb = org.Mm.eg.db, ont = "MF", level = 3, readable = TRUE)
barplot(ggoALL_MF3, showCategory = 40,  title = "GroupGO CA1 MF3")
ggoUP_MF3 <- groupGO(gene = nidUP$ENTREZID, OrgDb = org.Mm.eg.db, ont = "MF", level = 3, readable = TRUE)
barplot(ggoUP_MF3, showCategory = 40,  title = "GroupGO KK1 MF3")
ggoDOWN_MF3 <- groupGO(gene = nidDOWN$ENTREZID, OrgDb = org.Mm.eg.db, ont = "MF", level = 3, readable = TRUE)
barplot(ggoDOWN_MF3, showCategory = 40, title = "GroupGO CA05 MF2")
# Here one needs to play a lot with the level.


## Perform the enrichment in GO Molecular Functions analysis.
egoALL_MF <- enrichGO(gene = nidALL$ENTREZID, OrgDb = org.Mm.eg.db, minGSSize = 20, ont = "MF", pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)
dotplot(egoALL_MF, title = "EnrichGO ALL MF")
egoUP_MF <- enrichGO(gene = nidUP$ENTREZID, OrgDb = org.Mm.eg.db, minGSSize = 50, ont = "MF", pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)
dotplot(egoUP_MF, title = "EnrichGO UP MF")
egoDOWN_MF <- enrichGO(gene = nidDOWN$ENTREZID, OrgDb = org.Mm.eg.db, minGSSize = 30, ont = "MF", pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)
dotplot(egoDOWN_MF, title = "EnrichGO DOWN MF")



## Perform enrichment in GO Biological Processes analysis.
egoALL_BP <- enrichGO(gene = nidALL$ENTREZID, OrgDb = org.Mm.eg.db, minGSSize = 20, ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)
dotplot(egoALL_BP, title = "EnrichGO ALL BP")
egoUP_BP <- enrichGO(gene = nidUP$ENTREZID, OrgDb = org.Mm.eg.db, minGSSize = 50, ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)
dotplot(egoUP_BP, title = "EnrichGO UP BP")
egoDOWN_BP <- enrichGO(gene = nidDOWN$ENTREZID, OrgDb = org.Mm.eg.db, minGSSize = 30, ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)
dotplot(egoDOWN_BP, title = "EnrichGO DOWN BP")


## Perform gene set enrichment in GO Molecular Functions analysis.
egogsALL_MF <- gseGO(geneList = geneListNLS, OrgDb = org.Mm.eg.db, ont = "MF", nPerm = 1000, minGSSize = 50, maxGSSize = 500, pvalueCutoff = 0.05, verbose = FALSE, keyType = "SYMBOL")
dotplot(egogsALL_MF, title = "GSEA GO ALL MF")
egogsUP_MF <- gseGO(geneList = geneListNLSUP, OrgDb = org.Mm.eg.db, ont = "MF", nPerm = 1000, minGSSize = 50, maxGSSize = 500, pvalueCutoff = 0.05, verbose = FALSE, keyType = "SYMBOL")
dotplot(egogsUP_MF, title = "GSEA GO UP MF")
egogsDOWN_MF <- gseGO(geneList = geneListNLSDOWN, OrgDb = org.Mm.eg.db, ont = "MF", nPerm = 1000, minGSSize = 50, maxGSSize = 500, pvalueCutoff = 0.05, verbose = FALSE, keyType = "SYMBOL")
dotplot(egogsDOWN_MF, title = "GSEA GO DOWN MF")


## Perform the gene set enrichment in GO biological processes analysis.
egogsALL_BP <- gseGO(geneList = geneListNLS, OrgDb = org.Mm.eg.db, ont = "BP", nPerm = 1000, minGSSize = 50, maxGSSize = 500, pvalueCutoff = 0.05, verbose = FALSE, keyType = "SYMBOL")
dotplot(egogsALL_BP, title = "GSEA GO ALL BP")
egogsUP_BP <- gseGO(geneList = geneListNLSUP, OrgDb = org.Mm.eg.db, ont = "BP", nPerm = 1000, minGSSize = 50, maxGSSize = 500, pvalueCutoff = 0.05, verbose = FALSE, keyType = "SYMBOL")
dotplot(egogsUP_BP, title = "GSEA GO UP BP")
egogsDOWN_BP <- gseGO(geneList = geneListNLSDOWN, OrgDb = org.Mm.eg.db, ont = "BP", nPerm = 1000, minGSSize = 50, maxGSSize = 500, pvalueCutoff = 0.05, verbose = FALSE, keyType = "SYMBOL")
dotplot(egogsDOWN_BP, title = "GSEA GO DOWN BP")


## Perform KEEG pathway enrichment.
ekeALL <- enrichKEGG(gene = nidALL$ENTREZID, organism = "mmu", pvalueCutoff = 0.05)
barplot(ekeALL, title = "KEGG enrichment NLS_ALL")
ekeUP <- enrichKEGG(gene = nidUP$ENTREZID, organism = "mmu", pvalueCutoff = 0.05)
barplot(ekeUP, title = "KEGG enrichment NLS_UP")
ekeDOWN <- enrichKEGG(gene = nidDOWN$ENTREZID, organism = "mmu", pvalueCutoff = 0.05)
barplot(ekeDOWN, showCategory = 60, title = "KEGG enrichment NLS_DOWN")

ekmALL <- enrichMKEGG(gene = nidALL$ENTREZID, organism = "mmu")
barplot(ekmALL, title = "KEGG modules enrichment NLS_ALL")
ekmUP <- enrichMKEGG(gene = nidUP$ENTREZID, organism = "mmu")
barplot(ekmUP, showCategory = 30, title = "KEGG enrichment NLS_UP")
ekmDOWN <- enrichMKEGG(gene = nidDOWN$ENTREZID, organism = "mmu")
barplot(ekmDOWN, title = "KEGG modules enrichment NLS_DOWN")

epaALL <- enrichPathway(gene = nidALL$ENTREZID, organism = "mouse", pvalueCutoff = 0.05)
barplot(epaALL, showCategory = 30, title = "Pathways enrichment NLS_ALL")
epaUP <- enrichPathway(gene = nidUP$ENTREZID, organism = "mouse", pvalueCutoff = 0.05)
barplot(epaUP, title = "Pathways enrichment NLS_UP")
epaDOWN <- enrichPathway(gene = nidDOWN$ENTREZID, organism = "mouse", pvalueCutoff = 0.05)
barplot(epaDOWN, title = "Pathways enrichment CA05")



### Visualisations ----------------------------------------
edoALL <- enrichDGN(nidALL$ENTREZID)
edoUP <- enrichDGN(nidUP$ENTREZID)
edoDOWN <- enrichDGN(nidDOWN$ENTREZID)
# We can possibly do that if we transform to human orthologous genes.

## BArplots - Dotplots.
# for Jupyter use the following to produce better figures
#options(repr.plot.width=10,repr.plot.height=14))
barplot(ggoALL_MF3, drop = TRUE, showCategory = 200, title = "Barplot GO-MF3 NLS_ALL")
barplot(ggoUP_MF3, drop = TRUE, showCategory = 200, title = "Barplot GO-MF3 NLS_UP")
barplot(ggoDOWN_MF3, drop = TRUE, showCategory = 200, title = "Barplot GO-MF3 NLS_DOWN")
dotplot(egoALL_MF, showCategory = 30) + ggtitle("Dotplot for MF enrichment for NLS_ALL")
dotplot(egoALL_BP, showCategory = 50) + ggtitle("Dotplot for BP enrichment for NLS_ALL")
dotplot(egoUP_MF, showCategory = 50) + ggtitle("Dotplot for MF enrichment for NLS_UP")
dotplot(egoUP_BP, showCategory = 50) + ggtitle("Dotplot for BP enrichment for NLS_UP")
dotplot(egoDOWN_MF, showCategory = 50) + ggtitle("Dotplot for MF enrichment for NLS_DOWN")
dotplot(egoDOWN_BP, showCategory = 50) + ggtitle("Dotplot for BP enrichment for NLS_DOWN")

## Enrichment Map GO plots, these plots help to understand the relationships BETWEEN GOs.
# BF
emapplot(egoALL_BP) + ggtitle("EMAplot BP NLS_ALL")
emapplot(egoUP_BP) + ggtitle("EMAplot BP NLS_UP")
emapplot(egoDOWN_BP) + ggtitle("EMAplot BP NLS_DOWN")
# MF
emapplot(egoALL_MF) + ggtitle("EMAplot MF NLS_ALL")
emapplot(egoUP_MF) + ggtitle("EMAplot MF NLS_UP")
emapplot(egoDOWN_MF) + ggtitle("EMAplot MF NLS_DOWN")


## Category Network (CNET) plots (perhaps the most usefull!)
cnetplot(egoALL_MF, foldChange = geneListNLS, colorEdge = TRUE) + ggtitle("CNETplot MF ALL")
cnetplot(egoUP_MF, foldChange = geneListNLSUP, colorEdge = TRUE) + ggtitle("CNETplot MF UP")
cnetplot(egoDOWN_MF, foldChange = geneListNLSDOWN, colorEdge = TRUE) + ggtitle("CNETplot MF DOWN")
cnetplot(egoALL_BP, foldChange = geneListNLS, colorEdge = TRUE) + ggtitle("CNETplot BP ALL")
cnetplot(egoUP_BP, foldChange = geneListNLSUP, colorEdge = TRUE) + ggtitle("CNETplot BP UP")
cnetplot(egoDOWN_BP, foldChange = geneListNLSDOWN, colorEdge = TRUE) + ggtitle("CNETplot BP DOWN")


## GOplots (extended EMA plots, perhaps confusing)
goplot(egoALL_MF, showCategory = 30, geom = "text", alph = 0.50) + ggtitle("GOplot MF NLS_ALL")
goplot(egoUP_MF, showCategory = 30, geom = "text", alph = 0.50) + ggtitle("GOplot MF NLS_UP")
goplot(egoDOWN_MF, showCategory = 30, geom = "text", alph = 0.50) + ggtitle("GOplot MF NLS_DOWN")
goplot(egoALL_BP, showCategory = 20, geom = "text", alph = 0.50) + ggtitle("GOplot BP NLS_ALL")



## Perform GO gene set enrichment analysis.
# MF
edoALL_MF <- gseGO(geneListNLS, nPerm = 10000, OrgDb = org.Mm.eg.db, ont = "MF", keyType = "SYMBOL")
dotplot(edoALL_MF, title = "GSEA GO MF NLS_ALL")  #NO
edoUP_MF <- gseGO(geneListNLSUP, nPerm = 10000, OrgDb = org.Mm.eg.db, ont = "MF", keyType = "SYMBOL")
dotplot(edoUP_MF, title = "GSEA GO MF NLS_UP")  #NO
edoDOWN_MF <- gseGO(geneListNLSDOWN, nPerm = 10000, OrgDb = org.Mm.eg.db, ont = "MF", keyType = "SYMBOL")
dotplot(edoDOWN_MF, title = "GSEA GO MF NLS_DOWN")  #NO
edoALL_BP <- gseGO(geneListNLS, nPerm = 10000, OrgDb = org.Mm.eg.db, ont = "BP", keyType = "SYMBOL")
dotplot(edoALL_BP, title = "GSEA GO BP NLS_ALL")  #NO
edoUP_BP <- gseGO(geneListNLSUP, nPerm = 10000, OrgDb = org.Mm.eg.db, ont = "BP", keyType = "SYMBOL")
dotplot(edoUP_BP, title = "GSEA GO BP NLS_UP")  #NO
edoDOWN_BP <- gseGO(geneListNLSDOWN, nPerm = 10000, OrgDb = org.Mm.eg.db, ont = "BP", keyType = "SYMBOL")
dotplot(edoDOWN_BP, title = "GSEA GO BP NLS_DOWN")


# Heat plots, grouping and gene expression!
heatplot(edoDOWN_BP, foldChange = geneListNLSDOWN)

# Ridge plot (funcy)
ridgeplot(edoDOWN_BP)

# GSEA plots, USEFULL more difficult to interpret.
gseaplot2(edoDOWN_BP, geneSetID = 1, title = edoDOWN_BP$Description[1])
gseaplot2(edoDOWN_BP, geneSetID = 2, title = edoDOWN_BP$Description[2])
gseaplot2(edoDOWN_BP, geneSetID = 1:5, title = "GSEAs of the top 5 NLS_DOWN")  # This is usefull as we can superimpose many different enrichments!



### Clustering --------------------------------------------
# Collect data.
dd <- read.table("tpm_WT_NLSB.tab", header = TRUE)
ddNLS <- dd[namesNLS, ]
ddNLSUP <- dd[namesNLSUP, ]
ddNLSDOWN <- dd[namesNLSDOWN, ]


#Clean the zeros
dNLS <- ddNLS[apply(ddNLS, 1, function(row) any(row != 0 )),]
dNLSUP <- ddNLSUP[apply(ddNLSUP, 1, function(row) any(row != 0 )),]
dNLSDOWN <- ddNLSDOWN[apply(ddNLSDOWN, 1, function(row) any(row != 0 )),]
summary(dNLSUP)

# Get the logs. (natural) Perhaps we do not need that.
#dCA1 <- log(ddCA1[,c(1,2,3,4,5,6)])
#dKK1 <- log(ddKK1[,c(1,2,3,7,8,9)])
#dCA05 <- log(ddCA05[,c(1,2,3,4,5,6)])

colnames(dNLS) <- c("WT1", "WT2", "WT3", "NLS", "NLS", "NLS")
colnames(dNLSUP) <- c("WT1", "WT2", "WT3", "NLS", "NLS", "NLS")
colnames(dNLSDOWN) <- c("WT1", "WT2", "WT3", "NLS", "NLS", "NLS")


my_palette <- brewer.pal(n = 11, name = "RdYlGn")
### For ALL the NLS DEGs
## Hierarhical clustering of NLS DEGs
# Do the clustering.
clust_NLS <- hclust(dist(dNLS), method = "centroid")
# define clusters (hard thresold)
NLS_Clusts <- cutree(clust_NLS, k = 4)
# Colour vector for clusters side bar.
myClustCols <- rainbow(length(unique(NLS_Clusts)))
myClusters_NLS <- myClustCols[NLS_Clusts]
heatmap.2(as.matrix(dNLS), main = "Hierarchical clustering of NLS DEGs.", Rowv = as.dendrogram(clust_NLS), Colv = FALSE, dendrogram = "row", scale = "row", col = my_palette, cexCol = 1.5, cexRow = 0.4, key.title = NA, keysize = 0.8, key.xlab = NA, ylab = "Genes", RowSideColors = myClusters_NLS)


## Determine number of clusters by plotting sum of squares whitnin groups
wssNLS <- (nrow(dNLS))*sum(apply(dNLS,2,var))
for (i in 2:15) wssNLS[i] <- sum(kmeans(dNLS, centers = i, iter.max = 100, nstart = 250)$withinss)
my_palette <- brewer.pal(n = 11, name = "RdYlGn")
plot(1:15, wssNLS, type = "b", xlab = "Number of Clusters", ylab = "Within groups sum of squares CA")

## Perform k-means clustering
# Kmeans NLS ALL.
set.seed(1)  # Always set the same random seed.
kmNLS <- kmeans(as.matrix(dNLS), 2, iter.max = 200, nstart = 20)
# Append id and cluster
dfclNLS <- cbind(dNLS, id = seq(nrow(dNLS)), cluster = kmNLS$cluster)
# Add idsort, the id number ordered by cluster
dfclNLS$idsort <- dfclNLS$id[order(dfclNLS$cluster)]
dfclNLS$idsort <- order(dfclNLS$idsort)
clusterColsNLS <- as.character(sort(kmNLS$cluster))
## Plot k-means clustering.
heatmap(as.matrix(dNLS)[order(kmNLS$cluster),], Rowv = NA, col = my_palette, Colv = NA, cexCol = 1.5, cexRow = 0.4, RowSideColors = clusterColsNLS, ylab = "Genes", main = "k-means clustering of NLS tmp.")

# Kmeans NLS-UP.
kmNLSUP <- kmeans(as.matrix(dNLSUP), 4, iter.max = 200, nstart = 20)
# Append id and cluster
dfclNLSUP <- cbind(dNLSUP, id = seq(nrow(dNLSUP)), cluster = kmNLSUP$cluster)
# Add idsort, the id number ordered by cluster
dfclNLSUP$idsort <- dfclNLSUP$id[order(dfclNLSUP$cluster)]
dfclNLSUP$idsort <- order(dfclNLSUP$idsort)
clusterColsNLSUP <- as.character(sort(kmNLSUP$cluster))
## Plot k-means clustering.
heatmap(as.matrix(dNLSUP)[order(kmNLSUP$cluster),], Rowv = NA, col = my_palette, Colv = NA, cexCol = 1.5, cexRow = 0.4, RowSideColors = clusterColsNLSUP, ylab = "Genes", main = "k-means clustering of NLS-UP tmp.")


# The PAM method.
pam_NLS <- pam(dNLS, 4)
## Plot PAM clustering.
heatmap(as.matrix(dNLS)[order(kmNLSUP$cluster),], Rowv = NA, col = my_palette, Colv = NA, cexCol = 1.5, cexRow = 0.4, RowSideColors = clusterColsNLSUP, ylab = "Genes", main = "k-means clustering of NLS-UP tmp.")

## Perform some clustering tests and visulaisations.
fviz_cluster(pam_NLS, data = dNLS, ellipse.type = "convex") + theme_minimal()


# Lastly try also the MClust package.
mclust_NLS <- Mclust(dNLS)
summary(mclust_NLS)

fviz_cluster(mclust_NLS, data = dNLS, ellipse.type = "convex") + theme_minimal()

## Clustering raw data of DEGs is not so informative.
