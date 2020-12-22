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
library(RDAVIDWebService)
library(Factoshiny)
library(FactoMineR)
library(corrplot)
library(NetCluster)
library(eulerr)
library(ggcorrplot)


### Load/preprocess data ----------------------------------
# We have preselected for logFC1 and significant genes.
d <- read.table("DESeq2_Genes_NLS_vs_WT.tab", header = TRUE)

# Here are the TPMs of all genes.
dd <- read.table("tpm_WT_NLSB.tab", header = TRUE)

dnls <- d # For LogFC 1
dnls <- subset(d, (d$log2FoldChange >= 1.5 | d$log2FoldChange <= -1.5)) # For LogFC 1.5

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

# Get the universe. (i.e. all the genes that we have identify in our experiment.)
allMm_ProteinCoding <- scan("ensembl_Mouse_ProteinCoding_NAMES.txt", what = "character")
dU <- dd[allMm_ProteinCoding,, ]
# Remove zeros and NAs.
dU <- dU[apply(dU, 1, function(row) all(row !=0 )), ]
dU <- dU[complete.cases(dU),]
universeSYMBOLS <- rownames(dU)
universeENTREZ <- bitr(universeSYMBOLS, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)$ENTREZID



### Enrichments -------------------------------------------

# WikiPaths enrichment.
# Add column with the WikiGene ID.
convt <- read.table("wikiID_geneID.txt", sep = "\t")  # table downloaded from ENSEMBL BioMart
colnames(convt) <- c("WikiGeneID", "GeneID")
namesNLSwiki <- convt[convt$GeneID %in% namesNLS,]$WikiGeneID
namesNLSUPwiki <- convt[convt$GeneID %in% namesNLSUP,]$WikiGeneID
namesNLSDOWNwiki <- convt[convt$GeneID %in% namesNLSDOWN,]$WikiGeneID


# Prepare the reference pathways and the TERM2 objects.
wp2genes <- read.gmt("/home/costas/sysBiol_Diderot/common_data/wikiPatways_Mm/wikipathways-20200610-gmt-Mus_musculus.gmt")
wp2genes <- wp2genes %>% tidyr::separate(term, c("name", "version", "wpid", "org"), "%")
wpid2gene <- wp2genes %>% dplyr::select("wpid", "gene") #TERM2GENE
wpid2name <- wp2genes %>% dplyr::select("wpid", "name") #TERM2NAME


## Perform WikiPaths enrichment analysis.
ewpNLS <- enricher(namesNLSwiki, minGSSize = 10, TERM2GENE = wpid2gene, TERM2NAME = wpid2name)
barplot(ewpNLS, showCategory = 30, title = "WikiPaths enrichment NLS")
ewpNLSUP <- enricher(namesNLSUPwiki, minGSSize = 10, TERM2GENE = wpid2gene, TERM2NAME = wpid2name)
barplot(ewpNLSUP, showCategory = 30, title = "WikiPaths enrichment NLS-UP")
ewpNLSDOWN <- enricher(namesNLSDOWNwiki, minGSSize = 10, TERM2GENE = wpid2gene, TERM2NAME = wpid2name)
barplot(ewpNLS, showCategory = 30, title = "WikiPaths enrichment NLS-DOWN")


# Prepare the geneLists (ordered gene list)
# These are the gene names that we have Wiki IDs.
nn <- unique(convt[convt$WikiGene.ID %in% namesNLSwiki,]$GeneID)

glNLS <- geneListNLS[names(geneListNLS) %in% levels(nn)]

length(unique(convt[convt$GeneID %in% names(glNLS),]$WikiGeneID))
length(names(glNLS))
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


## Perform Gene Set Enrichment of WikiPaths.
ewpgsNLS <- GSEA(geneListNLS, TERM2GENE = wpid2geneID, TERM2NAME = wpid2name, verbose = FALSE)
barplot(ewpgsNLS, title = "WikiPaths GSEA NLS")  # No enrichment.


# MSIGDB enrichment
m_df <- msigdbr(species = "Mus musculus")
#m_df$gs_id <- m_df$gene_symbol # BIG TRICK TO SWAP THE COLUMN NAMES!!!!!!
m_t2g <- m_df %>% dplyr::select(gs_name, gene_symbol)

## Perform GSIGDB gene enrichment.
esigALL <- enricher(namesNLS, minGSSize = 20, TERM2GENE = m_t2g)
dotplot(esigALL, showCategory = 50, title = "MsiGDB enrichment ALL")
esigUP <- enricher(namesNLSUP, minGSSize = 20, TERM2GENE = m_t2g)
dotplot(esigUP, showCategory = 50, title = "MsiGDB enrichment UP")
esigDOWN <- enricher(namesNLSDOWN, minGSSize = 20, TERM2GENE = m_t2g)
dotplot(esigDOWN, showCategory = 50, title = "MsiGDB enrichment DOWN")

## Perform GSIGDB Gene Set Enrichment.
esigsALL <- GSEA(geneListNLS, minGSSize = 20, TERM2GENE = m_t2g)
dotplot(esigsALL, showCategory = 50, title = "MsiGDB GSEA ALL")
esigsUP <- GSEA(geneListNLSUP, minGSSize = 20, TERM2GENE = m_t2g)
dotplot(esigsUP, showCategory = 50, title = "MsiGDB GSEA UP")
esigsDOWN <- GSEA(geneListNLSDOWN, minGSSize = 20, TERM2GENE = m_t2g)
dotplot(esigsDOWN, showCategory = 50, title = "MsiGDB GSEA DOWN")

# Transform the common gene names to ENTREZIDs
nidALL <- bitr(namesNLS, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
nidUP <- bitr(namesNLSUP, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
nidDOWN <- bitr(namesNLSDOWN, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)

## Perform GO group enrichment analysis.
ggoALL_MF3 <- groupGO(gene = nidALL$ENTREZID, OrgDb = org.Mm.eg.db, ont = "MF", level = 3, readable = TRUE)
barplot(ggoALL_MF3, showCategory = 40,  title = "GroupGO NLS-ALL MF3")
ggoUP_MF3 <- groupGO(gene = nidUP$ENTREZID, OrgDb = org.Mm.eg.db, ont = "MF", level = 3, readable = TRUE)
barplot(ggoUP_MF3, showCategory = 40,  title = "GroupGO NLS-UP MF3")
ggoDOWN_MF3 <- groupGO(gene = nidDOWN$ENTREZID, OrgDb = org.Mm.eg.db, ont = "MF", level = 3, readable = TRUE)
barplot(ggoDOWN_MF3, showCategory = 40, title = "GroupGO NLS-DOWN MF3")
# Here one needs to play a lot with the level.


## Perform the enrichment in GO Molecular Functions analysis.
egoALL_MF <- enrichGO(gene = nidALL$ENTREZID, OrgDb = org.Mm.eg.db, minGSSize = 20, maxGSSize = 600, ont = "MF", pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)
dotplot(egoALL_MF, title = "EnrichGO ALL MF", showCategory = 30)
egoUP_MF <- enrichGO(gene = nidUP$ENTREZID, OrgDb = org.Mm.eg.db, minGSSize = 20, maxGSSize = 600, ont = "MF", pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)
dotplot(egoUP_MF, title = "EnrichGO UP MF", showCategory = 30)
egoDOWN_MF <- enrichGO(gene = nidDOWN$ENTREZID, OrgDb = org.Mm.eg.db, minGSSize = 20, maxGSSize = 600, ont = "MF", pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)
dotplot(egoDOWN_MF, title = "EnrichGO DOWN MF", showCategory = 30)  # Only this one yiled results


## Perform enrichment in GO Biological Processes analysis.
egoALL_BP <- enrichGO(gene = nidALL$ENTREZID, OrgDb = org.Mm.eg.db, minGSSize = 20, maxGSSize = 600, ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)
dotplot(egoALL_BP, showCategory = 30, title = "EnrichGO ALL BP")
cnetplot(egoALL_BP)
egoUP_BP <- enrichGO(gene = nidUP$ENTREZID, OrgDb = org.Mm.eg.db, minGSSize = 20, maxGSSize = 600, ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)
dotplot(egoUP_BP, title = "EnrichGO UP BP", showCategory = 30)
egoDOWN_BP <- enrichGO(gene = nidDOWN$ENTREZID, OrgDb = org.Mm.eg.db, minGSSize = 20, maxGSSize = 600, ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)
dotplot(egoDOWN_BP, title = "EnrichGO DOWN BP", showCategory = 30)



## Perform enrichment in GO ALL categories.
egoALL <- enrichGO(gene = nidALL$ENTREZID, OrgDb = org.Mm.eg.db, minGSSize = 20, maxGSSize = 500, ont = "ALL", pAdjustMethod = "BH", pvalueCutoff = 0.02, readable = TRUE)
dotplot(egoALL, showCategory = 40, title = "EnrichGO ALL ALL", color = "p.adjust", x = "GeneRatio")
cnetplot(egoALL, foldChange = geneListNLS, colorEdge = TRUE, showCategory = 25) + ggtitle("Gene-concept network plot NLS ALL")
# CNETplot for the paper.
cnetplot(egoALL, foldChange = geneListNLS, colorEdge = TRUE, showCategory = 16) #+ ggtitle("Gene-concept network plot NLS ALL")

egoUP_BP <- enrichGO(gene = nidUP$ENTREZID, OrgDb = org.Mm.eg.db, minGSSize = 20, maxGSSize = 600, ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)
dotplot(egoUP_BP, title = "EnrichGO UP BP", showCategory = 30)
egoDOWN_BP <- enrichGO(gene = nidDOWN$ENTREZID, OrgDb = org.Mm.eg.db, minGSSize = 20, maxGSSize = 800, ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)
dotplot(egoDOWN_BP, title = "EnrichGO DOWN BP", showCategory = 30)



## Perform gene set enrichment analysis in GO ALL.
egogsALL_ALL <- gseGO(geneList = geneListNLS, ont = "ALL", OrgDb = org.Mm.eg.db, minGSSize = 20, maxGSSize = 1000, pvalueCutoff = 0.1, keyType = "SYMBOL", by = "DOSE", nPerm = 1000)
dotplot(egogsALL_ALL, title = "GSEA GO ALL DEGs all")
egogsUP_ALL <- gseGO(geneList = geneListNLSUP, ont = "ALL", OrgDb = org.Mm.eg.db, minGSSize = 20, maxGSSize = 1000, pvalueCutoff = 0.1, keyType = "SYMBOL", by = "DOSE", nPerm = 1000)
dotplot(egogsUP_ALL, title = "GSEA GO ALL DEGs UP")
egogsDOWN_ALL <- gseGO(geneList = geneListNLSDOWN, ont = "ALL", OrgDb = org.Mm.eg.db, minGSSize = 20, maxGSSize = 1000, pvalueCutoff = 0.1, keyType = "SYMBOL", by = "DOSE", nPerm = 1000)
dotplot(egogsUP_ALL, title = "GSEA GO ALL DEGs DOWN")


## Perform the gene set enrichment in GO biological processes analysis.
egogsALL_BP <- gseGO(geneList = geneListNLS, OrgDb = org.Mm.eg.db, ont = "BP", minGSSize = 20, maxGSSize = 800, pvalueCutoff = 0.1, verbose = FALSE, keyType = "SYMBOL")
dotplot(egogsALL_BP, title = "GSEA GO ALL BP")
egogsUP_BP <- gseGO(geneList = geneListNLSUP, OrgDb = org.Mm.eg.db, ont = "BP", minGSSize = 50, maxGSSize = 500, pvalueCutoff = 0.05, verbose = FALSE, keyType = "SYMBOL")
dotplot(egogsUP_BP, title = "GSEA GO UP BP")
egogsDOWN_BP <- gseGO(geneList = geneListNLSDOWN, OrgDb = org.Mm.eg.db, ont = "BP", minGSSize = 50, maxGSSize = 500, pvalueCutoff = 0.05, verbose = FALSE, keyType = "SYMBOL")
dotplot(egogsDOWN_BP, title = "GSEA GO DOWN BP")


## Perform KEEG pathway enrichment.
ekeALL <- enrichKEGG(gene = nidALL$ENTREZID, organism = "mmu", pvalueCutoff = 0.1, minGSSize = 20, maxGSSize = 1000)
barplot(ekeALL, title = "KEGG enrichment NLS_ALL")
ekeUP <- enrichKEGG(gene = nidUP$ENTREZID, organism = "mmu", pvalueCutoff = 0.05, minGSSize = 20, maxGSSize = 1000)
barplot(ekeUP, title = "KEGG enrichment NLS_UP")
ekeDOWN <- enrichKEGG(gene = nidDOWN$ENTREZID, organism = "mmu", pvalueCutoff = 0.05, minGSSize = 20, maxGSSize = 1000)
barplot(ekeDOWN, showCategory = 60, title = "KEGG enrichment NLS_DOWN")

ekmALL <- enrichMKEGG(gene = nidALL$ENTREZID, organism = "mmu")
barplot(ekmALL, title = "KEGG modules enrichment NLS_ALL")
ekmUP <- enrichMKEGG(gene = nidUP$ENTREZID, organism = "mmu")
barplot(ekmUP, showCategory = 30, title = "KEGG enrichment NLS_UP")
ekmDOWN <- enrichMKEGG(gene = nidDOWN$ENTREZID, organism = "mmu")
barplot(ekmDOWN, title = "KEGG modules enrichment NLS_DOWN")



# Pathways enrichment!
epaALL <- enrichPathway(gene = nidALL$ENTREZID, organism = "mouse", pvalueCutoff = 0.05, universe = universeENTREZ, minGSSize = 20, maxGSSize = 800, readable = TRUE)
barplot(epaALL, showCategory = 30, title = "Pathways enrichment NLS_ALL")
dotplot(epaALL, showCategory = 30, title = "Pathways enrichment NLS_ALL")
epaUP <- enrichPathway(gene = nidUP$ENTREZID, organism = "mouse", pvalueCutoff = 0.05)
barplot(epaUP, title = "Pathways enrichment NLS_UP")
epaDOWN <- enrichPathway(gene = nidDOWN$ENTREZID, organism = "mouse", pvalueCutoff = 0.05)
barplot(epaDOWN, title = "Pathways enrichment NLS-DOWN")



### Visualisations ----------------------------------------
edoALL <- enrichDGN(nidALL$ENTREZID)
edoUP <- enrichDGN(nidUP$ENTREZID)
edoDOWN <- enrichDGN(nidDOWN$ENTREZID)
# We can possibly do that if we transform to human orthologous genes.

## Barplots - Dotplots.
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
emapplot(egoALL_BP) + ggtitle("GO network plot BP NLS-ALL")
emapplot(egoUP_BP) + ggtitle("GO network plot BP NLS-UP")
emapplot(egoDOWN_BP) + ggtitle("GO network plot BP NLS-DOWN")
# MF
emapplot(egoALL_MF) + ggtitle("GO network plot MF NLS-ALL")
emapplot(egoUP_MF) + ggtitle("GO network plot MF NLS-UP")
emapplot(egoDOWN_MF) + ggtitle("GO network plot MF NLS-DOWN")



## Category Network (CNET) plots (perhaps the most usefull!)
cnetplot(egoALL_MF, foldChange = geneListNLS, colorEdge = TRUE, showCategory = 20) + ggtitle("Gene-concept network plot MF NLS-ALL")
dotplot(egoALL_MF, showCategory = 20, color = "p.adjust", title = "Dotplot GO MF NLS-ALL")
# NO MF enrichment for NLS-UP.
cnetplot(egoDOWN_MF, foldChange = geneListNLSDOWN, colorEdge = TRUE, showCategory = 20) + ggtitle("Gene-concept Network plot MF NLS-DOWN")
dotplot(egoDOWN_MF, showCategory = 20, color = "p.adjust", title = "Dotplot GO MF NLS-DOWN")
cnetplot(egoALL_BP, foldChange = geneListNLS, colorEdge = TRUE, showCategory = 20) #+ ggtitle("Gene-concept network plot BP NLS-ALL")
dotplot(egoALL_BP, showCategory = 20, color = "p.adjust", title = "Dotplot GO BP NLS-ALL")
cnetplot(egoUP_BP, foldChange = geneListNLSUP, colorEdge = TRUE, showCategory = 20) + ggtitle("Gene-concept network plot BP NLS-UP")
dotplot(egoUP_BP, showCategory = 20, color = "p.adjust", title = "Dotplot GO BP NLS-UP")
cnetplot(egoDOWN_BP, foldChange = geneListNLSDOWN, colorEdge = TRUE, showCategory = 20) + ggtitle("Gene-concept network plot BP NLS-DOWN")
dotplot(egoDOWN_BP, showCategory = 20, color = "p.adjust", title = "Dotplot GO BP NLS-DOWN")


## GOplots (extended EMA plots, perhaps confusing)
goplot(egoALL_MF, showCategory = 30, geom = "text", alph = 0.50) + ggtitle("GOplot MF NLS_ALL")
goplot(egoUP_MF, showCategory = 30, geom = "text", alph = 0.50) + ggtitle("GOplot MF NLS_UP")
goplot(egoDOWN_MF, showCategory = 30, geom = "text", alph = 0.50) + ggtitle("GOplot MF NLS_DOWN")
goplot(egoALL_BP, showCategory = 20, geom = "text", alph = 0.50) + ggtitle("GOplot BP NLS_ALL")



## Perform GO gene set enrichment analysis.
# MF
edoALL_MF <- gseGO(geneListNLS, OrgDb = org.Mm.eg.db, ont = "MF", keyType = "SYMBOL")
dotplot(edoALL_MF, title = "GSEA GO MF NLS_ALL")  #NO
edoUP_MF <- gseGO(geneListNLSUP, OrgDb = org.Mm.eg.db, ont = "MF", keyType = "SYMBOL")
dotplot(edoUP_MF, title = "GSEA GO MF NLS_UP")  #NO
edoDOWN_MF <- gseGO(geneListNLSDOWN, OrgDb = org.Mm.eg.db, ont = "MF", keyType = "SYMBOL")
dotplot(edoDOWN_MF, title = "GSEA res.pca.CountsGO MF NLS_DOWN")  #NO
# BP
edoALL_BP <- gseGO(geneListNLS, minGSSize = 20, maxGSSize = 1000, OrgDb = org.Mm.eg.db, ont = "BP", keyType = "SYMBOL", pvalueCutoff = 0.1)
dotplot(edoALL_BP, title = "GSEA GO BP NLS_ALL")  #NO
edoUP_BP <- gseGO(geneListNLSUP, nPerm = 10000, OrgDb = org.Mm.eg.db, ont = "BP", keyType = "SYMBOL")
dotplot(edoUP_BP, title = "GSEA GO BP NLS_UP")  #NO
edoDOWN_BP <- gseGO(geneListNLSDOWN, OrgDb = org.Mm.eg.db, ont = "BP", keyType = "SYMBOL", minGSSize = 20, maxGSSize = 1000, pvalueCutoff = 0.1)
dotplot(edoDOWN_BP, title = "GSEA GO BP NLS_DOWN", showCategory = 20)


# Heat plots, grouping and gene expression!
heatplot(edoDOWN_BP, foldChange = geneListNLSDOWN)

# Ridge plot (funcy)
ridgeplot(edoDOWN_BP)

# GSEA plots, USEFULL more difficult to interpret.
gseaplot2(edoDOWN_BP, geneSetID = 1, title = edoDOWN_BP$Description[1])
gseaplot2(edoDOWN_BP, geneSetID = 2, title = edoDOWN_BP$Description[2])
gseaplot2(edoDOWN_BP, geneSetID = 1:5, title = "GSEAs of the top 5 NLS_DOWN")  # This is useful as we can superimpose many different enrichments!


### Quality controls ------------   --------------------------
tpmNLS <- read.table("tpm_WT_NLSB.tab", header = TRUE)
groupsall <- factor(c("WT", "WT", "WT", "NLS", "NLS", "NLS"))
## PCA on TPMs
res.pca.CPMs = PCA(t(tpmNLS), graph = FALSE)
fviz_pca_ind(res.pca.CPMs,
             fill.ind = groupsall, col.var = "black", repel = TRUE,
             col.ind = groupsall, # colored by groups
             palette = c("#00AFBB", "#E7B800"),
             addEllipses = TRUE, # Ellipses
             legend.title = "Groups",
             title = "PCA plot of CPMs WT-NLS.")

# On the detected protein coding genes only!
res.pca.TPMs = PCA(t(dU), graph = FALSE)
fviz_pca_ind(res.pca.TPMs,
             fill.ind = groupsall, col.var = "black", repel = TRUE,
             col.ind = groupsall, # colored by groups
             palette = c("#00AFBB", "#E7B800"),
             addEllipses = TRUE, # Ellipse
             legend.title = "Groups",
             title = "PCA plot of TPMs WT-NLS.")


### Clustering --------------------------------------------
# Collect data.
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

colnames(dNLS) <- c("WT1", "WT2", "WT3", "NLS1", "NLS2", "NLS3")
colnames(dNLSUP) <- c("WT1", "WT2", "WT3", "NLS1", "NLS2", "NLS3")
colnames(dNLSDOWN) <- c("WT1", "WT2", "WT3", "NLS1", "NLS2", "NLS3")


my_palette <- brewer.pal(n = 11, name = "RdYlGn")
### For ALL the NLS DEGs
## Hierarchical clustering of NLS DEGs
# Do the clustering.
clust_NLS <- hclust(dist(dNLS), method = "average")
# define clusters (hard thresold)
NLS_Clusts <- cutree(clust_NLS, k = 4)
# Colour vector for clusters side bar.
myClustCols <- rainbow(length(unique(NLS_Clusts)))
myClusters_NLS <- myClustCols[NLS_Clusts]
heatmap.2(as.matrix(dNLS), main = "Hierarchical clustering of NLS DEGs.", Rowv = as.dendrogram(clust_NLS), Colv = FALSE, dendrogram = "row", scale = "row", col = my_palette, cexCol = 1.5, cexRow = 0.4, key.title = NA, keysize = 0.8, key.xlab = NA, ylab = "Genes", RowSideColors = myClusters_NLS)


## Determine number of clusters by plotting sum of squares within groups
wssNLS <- (nrow(dNLS))*sum(apply(dNLS,2,var))
for (i in 2:15) wssNLS[i] <- sum(kmeans(dNLS, centers = i, iter.max = 100, nstart = 250)$withinss)
my_palette <- brewer.pal(n = 11, name = "RdYlGn")
plot(1:15, wssNLS, type = "b", xlab = "Number of Clusters", ylab = "Within groups sum of squares CA")

## Perform k-means clustering
my_palette <- colorRampPalette(brewer.pal(n = 11, name = "RdYlGn"))
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
heatmap(as.matrix(dNLS)[order(kmNLS$cluster),], Rowv = NA, col = my_palette(16), Colv = NA, cexCol = 1.5, cexRow = 0.3, ylab = "Genes", main = "k-means (2) clustering of NLS TPM.")
# For the paper
heatmap(as.matrix(dNLS)[order(kmNLS$cluster),], Rowv = NA, col = my_palette(16), Colv = NA, cexCol = 1.5, cexRow = 0.8, ylab = "Genes", margins = c(0.5, 5.2), labCol = "")
legend("topleft", legend = round(seq(range(dNLS)[1], range(dNLS)[2], length.out =13)), fill = my_palette(16), border = FALSE, cex = 1, bty = "n", inset = 0.04, title = "TPM")

# Kmeans NLS-UP.
kmNLSUP <- kmeans(as.matrix(dNLSUP), 4, iter.max = 200, nstart = 20)
# Append id and cluster
dfclNLSUP <- cbind(dNLSUP, id = seq(nrow(dNLSUP)), cluster = kmNLSUP$cluster)
# Add idsort, the id number ordered by cluster
dfclNLSUP$idsort <- dfclNLSUP$id[order(dfclNLSUP$cluster)]
dfclNLSUP$idsort <- order(dfclNLSUP$idsort)
clusterColsNLSUP <- as.character(sort(kmNLSUP$cluster))
## Plot k-means clustering.
heatmap(as.matrix(dNLSUP)[order(kmNLSUP$cluster),], Rowv = NA, col = my_palette(16), Colv = NA, cexCol = 1.5, cexRow = 0.4, RowSideColors = clusterColsNLSUP, ylab = "Genes", main = "k-means clustering of NLS-UP tmp.")


# The PAM method.
pam_NLS <- pam(dNLS, 2)
## Plot PAM clustering.
heatmap(as.matrix(dNLS)[order(kmNLSUP$cluster),], Rowv = NA, col = my_palette(16), Colv = NA, cexCol = 1.5, cexRow = 0.4, RowSideColors = clusterColsNLSUP, ylab = "Genes", main = "k-means clustering of NLS-UP tmp.")

## Perform some clustering tests and visulaisations.
fviz_cluster(pam_NLS, data = dNLS, ellipse.type = "convex") + theme_minimal()


# Lastly try also the MClust package.
mclust_NLS <- Mclust(dNLS)
summary(mclust_NLS)

fviz_cluster(mclust_NLS, data = dNLS, ellipse.type = "convex") + theme_minimal()
# Clustering raw data of DEGs is not so informative.



## Correlation - co-expression studies. -----------------
corrDNLS <- cor(t(dNLS))
# Soft threshold, power.
corrDNLS_thres <- corrDNLS**9
corrplot(corrDNLS_thres, diag = FALSE,
         order = "hclust", hclust.method = "ward.D2",
         cl.pos = "n", tl.cex = 0.7,
         main = "GE clustered correlation matrix of DEGs",
         mar = c(0,0,1,0), cex.main = 0.75, sig.level = 0.001, insig = "blank",
         addgrid.col = NA)
# Hard threshold 0.9 correlation.
corrDNLS_thresH <- corrDNLS_thres
corrDNLS_thresH[abs(corrDNLS_thresH) < 0.9 | corrDNLS_thresH == 1.0 ] = 0
idx <- colnames(corrDNLS_thresH[,!!colSums(corrDNLS_thresH)])
corrDNLS_thresHc <- corrDNLS_thresH[idx, idx]
corrplot(corrDNLS_thresHc, diag = FALSE,
         order = "hclust", hclust.method = "ward.D2",
         cl.pos = "n", tl.cex = 0.7,
         main = "GE clustered correlation matrix of DEGs HARD",
         mar = c(0,0,1,0), cex.main = 0.75, sig.level = 0.001, insig = "blank",
         addgrid.col = NA)

# Do the correlation clustering
corrDNLS_clust <- hclust(dist(corrDNLS_thresHc), method = "ward.D2")
clustsDNS <- cutree(corrDNLS_clust, k = 5)
# Put this in a data frame to produce the supplementary table.
clustDNSdf <- data.frame(row.names = names(sort(clustsDNS)), "Cluster" = sort(clustsDNS))
write.table(clustDNSdf, file = "degCorrelationClusters.csv", quote = FALSE)

clustCorrDNS <- clusterCorr(corrDNLS_thresHc, clustsDNS)
clustCorrMatDNLS <- generate_cluster_cor_mat(corrDNLS_thresHc, clustsDNS)
# these two are the SAME!

corrplot(clustCorrDNS, diag = FALSE,
         order = "hclust", hclust.method = "ward.D2",
         cl.pos = "n", tl.cex = 0.7,
         main = "GE clustered correlation matrix of DEGs Clustered",
         mar = c(0,0,1,0), cex.main = 0.75, sig.level = 0.001, insig = "blank",
         addgrid.col = NA)


# By observing the cluster correlation we found that cluster1 does not contain any informative correlation so we will remove it and re-plot the clustered correlation matrix.
clustDNLSclean <- clustsDNS[clustsDNS %in% c(2,3,4,5)]
corrDNLS_clean <- corrDNLS_thresHc[names(clustDNLSclean), names(clustDNLSclean)]

clustCorrDNLS_clean <- clusterCorr(corrDNLS_clean, clustDNLSclean)
corrplot(clustCorrDNLS_clean, diag = FALSE,
         order = "hclust", hclust.method = "ward.D2",
         tl.cex = 0.5, tl.col = "black", tl.pos = "t", tl.srt = 45,
         cl.pos = "r", cl.ratio = 0.1, cl.offset = 0.2, cl.align.text = "l",
         main = "GE clustered correlation matrix of DEGs Clustered_clean",
         mar = c(0,0,0,0), cex.main = 0.5, sig.level = 0.001, insig = "blank",
         addgrid.col = NA, addrect = 4)
# Plot for the paper figure.

# Get the HRNPC targets.
targets_3pUTR <- scan("utr_Analyses/anal15lFC/rbpMAP_UTRs/hrnpc_3pUTR_targets_NAMES.txt", what = "char")
targets_5pUTR <- scan("utr_Analyses/anal15lFC/rbpMAP_UTRs/hrnpc_5pUTR_targets_NAMES.txt", what = "char")

# Get the euler Venn diagram.
targets_UTR <- list("3'UTRs" = targets_3pUTR, "5'UTRs" = targets_5pUTR)
plot(euler(targets_UTR, shape = "ellipse"), quantities = TRUE)

# Retrieve the cluster names.
cluster2 <- names(clustsDNS[clustsDNS == 2])
cluster3 <- names(clustsDNS[clustsDNS == 3])
cluster4 <- names(clustsDNS[clustsDNS == 4])
cluster5 <- names(clustsDNS[clustsDNS == 5])

clusts_3pUTRs <- list("3'UTRs" = targets_3pUTR, "Clust2" = cluster2, "Clust3" = cluster3, "Clust4" = cluster4, "Clust5" = cluster5)
plot(euler(clusts_3pUTRs, shape = "ellipse"), quantities = TRUE)

clusts_5pUTRs <- list("5'UTRs" = targets_5pUTR, "Clust2" = cluster2, "Clust3" = cluster3, "Clust4" = cluster4, "Clust5" = cluster5)
plot(euler(clusts_5pUTRs, shape = "ellipse"), quantities = TRUE)

# GO enrichments of UTR HRNPC targets.
target3UTRsENTREZ <- nidALL[nidALL$SYMBOL %in% targets_3pUTR,]$ENTREZID
egoTarg3UTR <- enrichGO(gene = target3UTRsENTREZ, OrgDb = org.Mm.eg.db, minGSSize = 10, maxGSSize = 500, ont = "ALL", pAdjustMethod = "BH", pvalueCutoff = 0.01, readable = TRUE)
dotplot(egoTarg3UTR, title = "GO ALL enrichment of 3'UTR HRNPC targets", showCategory = 30)



## GO Enrichment of clusters ------------
cluster3ENTREZ <- nidALL[nidALL$SYMBOL %in% cluster3,]$ENTREZID
cluster4ENTREZ <- nidALL[nidALL$SYMBOL %in% cluster4,]$ENTREZID
cluster5ENTREZ <- nidALL[nidALL$SYMBOL %in% cluster5,]$ENTREZID
cluster2ENTREZ <- nidALL[nidALL$SYMBOL %in% cluster2,]$ENTREZID


# Perform the enrichment analysis.
egoClust4 <- enrichGO(gene = cluster4ENTREZ, OrgDb = org.Mm.eg.db, minGSSize = 10, maxGSSize = 500, ont = "ALL", pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)
dotplot(egoClust4, title = "GO ALL enrichment of Clust4", showCategory = 30)
egoClust3 <- enrichGO(gene = cluster3ENTREZ, OrgDb = org.Mm.eg.db, minGSSize = 10, maxGSSize = 500, ont = "ALL", pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)
dotplot(egoClust3, title = "GO ALL enrichment of Clust3", showCategory = 30)
egoClust5 <- enrichGO(gene = cluster5ENTREZ, OrgDb = org.Mm.eg.db, minGSSize = 10, maxGSSize = 500, ont = "ALL", pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)
dotplot(egoClust5, title = "GO ALL enrichment of Clust5", showCategory = 30)
egoClust2 <- enrichGO(gene = cluster2ENTREZ, OrgDb = org.Mm.eg.db, minGSSize = 10, maxGSSize = 200, ont = "ALL", pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)
dotplot(egoClust2, title = "GO ALL enrichment of Clust2", showCategory = 30)
# For the paper.
dotplot(egoClust2, showCategory = 16)



# Get the TPM DEGs table
tpmDEGSs <- tpmNLS[namesNLS,]
write.table(tpmDEGSs, "tpm_NLS_DEGs15lfc.tab", sep = "\t", quote = FALSE)









# TRIM71 analysis from ANOTHER PAPER -------------
namesTRIM71 <- scan("diffExpr_mircroarray_TRIM71KD.txt", what = "character")
trim71all <- bitr(namesTRIM71, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)

egoTrimALL_BP <-  enrichGO(gene = trim71all$ENTREZID, OrgDb = org.Mm.eg.db, minGSSize = 20, ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)
egoTrimALL_MF <-  enrichGO(gene = trim71all$ENTREZID, OrgDb = org.Mm.eg.db, minGSSize = 20, ont = "MF", pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE) # No enrichment.

dotplot(egoTrimALL_BP) #, foldChange = geneListNLS, colorEdge = TRUE, showCategory = 20) + ggtitle("Gene-concept network plot MF NLS-ALL")
