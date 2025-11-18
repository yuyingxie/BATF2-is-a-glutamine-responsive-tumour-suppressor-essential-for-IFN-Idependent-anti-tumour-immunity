library(edgeR)
library(ggplot2)
library(readxl)
library(dplyr)
library(msigdb)
library(msigdbr)
library(gprofiler2)

Project_Path <- "~/Leo/MD_Anderson/10000_WG_Analysis"

Counts <- read.delim(file = file.path(Project_Path, "counts/gene_expected_count.annot.txt"), header = TRUE, sep = "\t",
                     row.names = NULL, stringsAsFactors = FALSE)
Counts$external_gene_name <- ifelse(Counts$external_gene_name==".", Counts$gene_id, Counts$external_gene_name)
rownames(Counts) <- Counts$gene_id
Counts_data <- Counts[,5:34]
colnames(Counts_data) <- gsub("\\.", "_", colnames(Counts_data))

TPM <- read.delim(file = file.path(Project_Path, "counts/gene_TPM.annot.txt"), header = TRUE, sep = "\t",
                  row.names = NULL, stringsAsFactors = FALSE)
TPM$external_gene_name <- ifelse(TPM$external_gene_name==".", TPM$gene_id, TPM$external_gene_name)
rownames(TPM) <- TPM$gene_id
TPM_data <- TPM[,5:34]
colnames(TPM_data) <- gsub("\\.", "_", colnames(TPM_data))


Patient_info <- read_excel(file.path(Project_Path, "HNSCC_patient_info.xlsx"))

clinical_info <- data.frame(Treatment=rep(c("Rep1", "Rep2"), 15),
                            Patient=rep(paste0("P", 1:15),each=2),
                            Sample=paste0("X10000_WG_", 1:30),
                            Response=rep("Pos", 30),
                            Response10=rep("Pos", 30))


# sample 31 is same as sample 10
clinical_info[10,3] <- "X10000_WG_31"
clinical_info[c(7,8,13,14,23,24), 4] <- "Neg"
clinical_info[c(7,8,11,12,13,14,21,22,23,24,27,28), 5] <- "Neg"
clinical_info[] <- lapply(clinical_info, factor)
clinical_info$summary <- paste(clinical_info$Patient, clinical_info$Treatment, clinical_info$Response, sep = "_")

Counts_data <- Counts_data[, clinical_info$Sample]
colnames(Counts_data) <- clinical_info$summary
TPM_data <- TPM_data[, clinical_info$Sample]
colnames(TPM_data) <- clinical_info$summary




####################################################################
# filtering for the protein coding genes only
library(biomaRt) 
mart <- useEnsembl(biomart = "ensembl", mirror = "useast", dataset = "hsapiens_gene_ensembl")
genes <- biomaRt::getBM(attributes = c("external_gene_name", "ensembl_gene_id", "transcript_biotype"),
                        filters = c("transcript_biotype"), values = list("protein_coding"), mart = mart)
rownames(genes) <- genes$ensembl_gene_id
ensembl_protein_coding <- intersect(rownames(Counts_data), genes$ensembl_gene_id)

Counts_protein_coding <- Counts_data[ensembl_protein_coding,]
TPM_protein_coding <- TPM_data[ensembl_protein_coding,]

Counts_protein_coding <- Counts_protein_coding[-which(rowSums(Counts_protein_coding)==0),]
TPM_protein_coding <- TPM_protein_coding[-which(rowSums(TPM_protein_coding)==0),]

####################################################################
####################################################################
# Correlation analysis for BATF2 and Glutamine metabolism markers
# correlation plot using counts data
BATF2_Glutamine <- c("SLC7A5", "MYC", "SLC38A1", "GLS", "GPT2", "GLUD1", "GOT2", "SLC1A5", "SLC38A2", "IFNB1",
                     "IFIT2", "OASL", "ISG15", "IFNG", "CXCL9", "CXCL10", "TNF", "BATF2")

BATF2_Glutamine_EID <- rownames(Counts %>% filter(external_gene_name %in% BATF2_Glutamine))
BATF2_Glutamine_Counts <- Counts[BATF2_Glutamine_EID,]
rownames(BATF2_Glutamine_Counts) <- BATF2_Glutamine_Counts$external_gene_name
BATF2_Glutamine_Counts <- BATF2_Glutamine_Counts[BATF2_Glutamine,]

cor_mat_Counts <- cor(t(BATF2_Glutamine_Counts[,5:34]))
cor_mat_Counts[,"BATF2"]
pdf(file = "~/Documents/10000_WG_Analysis/Glutamine_BATF2_corrplot_counts.pdf", width = 8, height = 8)
corrplot::corrplot(
  cor_mat_Counts, type = 'lower', order = 'original', 
  tl.col = 'black', tl.srt = 90, col=colorRampPalette(c("blue","red"))(100)
)
dev.off()

# correlation plot using TPM data
BATF2_Glutamine <- c("MYC", "SLC7A5", "GLUD1", "SLC1A5", "GOT2", "SLC38A1", "GLS", "SLC38A2", "GPT2", "CXCL9", "CXCL10",
                     "IFIT2", "OASL", "ISG15", "IFNB1", "TNF", "IFNG", "BATF2")
BATF2_Glutamine_TPM <- TPM[BATF2_Glutamine_EID,]
rownames(BATF2_Glutamine_TPM) <- BATF2_Glutamine_TPM$external_gene_name
BATF2_Glutamine_TPM <- BATF2_Glutamine_TPM[BATF2_Glutamine,]

cor_mat_TPM <- cor(t(BATF2_Glutamine_TPM[,5:34]))
cor_mat_TPM[,"BATF2"]
pdf(file = "~/Documents/10000_WG_Analysis/Glutamine_BATF2_corrplot_TPM.pdf", width = 8, height = 8)
corrplot::corrplot(
  cor_mat_TPM, type = 'lower', order = 'original', 
  tl.col = 'black', tl.srt = 90, col=colorRampPalette(c("blue","red"))(100)
)
dev.off()



####################################################################################################################
# What are the significant genes different between Positive and Negative response using protein coding genes only
colSums(Counts_protein_coding)
y <- DGEList(counts = Counts_protein_coding, group = clinical_info$Response)
keep_Counts <- filterByExpr(y)
y <- y[keep_Counts, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
y$samples
dsgn <- model.matrix(~ clinical_info$Response)
rownames(dsgn) <- colnames(y)
fit <- glmQLFit(y, dsgn)
qlf <- glmQLFTest(fit)
Sig_Genes_Counts <- topTags(qlf, n = 19718, adjust.method = "BH")
# BY adjustment is valid under general dependence of multi-testing, produces larger adjusted pvalues
Sig_Genes_Counts
table(Sig_Genes_Counts$table$FDR <= 0.1)
DE_genes_counts <- Counts[rownames(Sig_Genes_Counts$table)[1:72],3]
DE_genes_counts


################################################################
# GSEA
#token <- upload_GMT_file(gmtfile = "~/Downloads/Human_GOBP_AllPathways_no_GO_iea_June_24_2019_symbol_max250gssize.gmt")
token <- "gp__eI5D_xEKK_i8A"
figure_function = function(dat){
  p=ggplot(data=datt, aes(x=forcats::fct_reorder(term_name,log10p),
                          y=log10p,fill=source,label = ifelse(significant ==TRUE, "*",""))) +
    geom_bar(stat="identity", position=position_dodge(),color=NA,width=0.8)+
    #scale_fill_continuous(low='#F0E442', high='red', guide=guide_colorbar(reverse=TRUE)) + 
    scale_fill_manual(values = c("#f6bd60", "#f5cac3", "#bfcc94", "#f28482", "#a11231"))+
    labs(title=paste0(datt$Title, " in 10000-WG HNSC"),x=paste0(""), y = "-log10(p value)")+
    coord_flip()+
    theme_classic()+
    #geom_text(vjust = 0.5) +
    geom_hline(yintercept=-log10(0.1), linetype="dashed", color = "red",size=2)+
    geom_text(vjust = 1, nudge_y = 0.5, size=15)+
    #ylim(0,1)+
    theme(plot.title = element_text(size = 30, color="black"),
          text = element_text(size = 30, color="black"),
          #axis.title = element_text(size = 25, color="black", face="bold"),
          axis.text.x=element_text(size = 30, color="black") ,
          legend.position = "none")# +
  return(p)
}

bg <- Counts$external_gene_name
source = "GO:BP"
gostres <- gost(query = DE_genes_counts, 
                organism = token, ordered_query = FALSE, 
                multi_query = FALSE, significant = FALSE, exclude_iea = TRUE, 
                measure_underrepresentation = FALSE, evcodes = FALSE, 
                user_threshold = 0.1, correction_method = "g_SCS", 
                domain_scope = "custom", custom_bg = bg, 
                numeric_ns = "",  
                as_short_link = FALSE)
gostres$result <- gostres$result[order(gostres$result$p_value),]
datt <- gostres$result[1:10,]
datt$log10p <- -log10(datt$p_value)
head(datt)
datt$Title <- "GSEA"
source <- gsub(":", "_", source)
datt <- data.frame(datt)
print(datt[1,c(3,10,11)])
plt <- figure_function(datt)
pdf(file = "~/Documents/10000_WG_Analysis/GSEA_DE_Genes_edgeR_Counts.pdf", width = 16, height = 8)
plt
dev.off()

plt
ggsave("~/Documents/10000_WG_Analysis/GSEA_DE_Genes_edgeR_Counts.png", width = 16, height = 8, units = "in", dpi = 300)
#############################################################################
# volcano plots
# Assuming 'data' is a data frame with 'log2FoldChange' and 'adj.P.Val' columns
library(EnhancedVolcano)

Sig_genes_names <- genes %>%
  filter(ensembl_gene_id %in% rownames(Sig_Genes_Counts))
Sig_genes_names <- Sig_genes_names[rownames(Sig_Genes_Counts),]
Sig_genes_names$external_gene_name2 <- ifelse(Sig_genes_names$external_gene_name=="",
                                              Sig_genes_names$ensembl_gene_id, Sig_genes_names$external_gene_name)

EnhancedVolcano(Sig_Genes_Counts$table, x="logFC", y="FDR", lab = Sig_genes_names$external_gene_name,
                pCutoff = 0.1, ylim = c(0, 3.5), xlim = c(-4, 6),
                ylab = bquote(~-Log[10] ~ italic(FDR)),
                title = "10000 WG", subtitle = "Volcano Plot",
                legendLabels = c("NS", expression(Log[2] ~ FC), "FDR", expression(FDR ~ and ~ log[2] ~ FC)))
ggsave("~/Documents/10000_WG_Analysis/VolcanoPlot_10000_WG.png", width = 11, height = 8, units = "in", dpi = 300)

#############################################################################
# Saving the logFC and FDR for the significant genes
Sig_Gene_df <- Sig_Genes_Counts$table %>%
  filter(FDR<0.1) %>%
  filter(abs(logFC)>1)

Sig_Gene_df$external_gene_name <- genes[rownames(Sig_Gene_df), 1]
Sig_Gene_df <- Sig_Gene_df %>% arrange(-logFC)
Sig_Gene_df <- Sig_Gene_df[,c(6,1:5)]
writexl::write_xlsx(Sig_Gene_df, path = "~/Documents/10000_WG_Analysis/DE_Genes_10000WG.xlsx")



#############################################################################
# edgeR sums over technical replicates
# Counts_Pooled_protein_coding <- sumTechReps(Counts_protein_coding, ID = clinical_info$Patient)
# 
# y <- DGEList(counts = Counts_Pooled_protein_coding, group = clinical_info$Response[2*(1:15)])
# keep_Counts <- filterByExpr(y)
# y <- y[keep_Counts, , keep.lib.sizes=FALSE]
# y <- calcNormFactors(y)
# y$samples
# dsgn <- model.matrix(~ clinical_info$Response[2*(1:15)])
# rownames(dsgn) <- colnames(y)
# fit <- glmQLFit(y, dsgn)
# qlf <- glmQLFTest(fit)
# Sig_Genes_Counts <- topTags(qlf, n = 19718, adjust.method = "BH")
# # BY adjustment is valid under general dependence of multi-testing, produces larger adjusted pvalues
# Sig_Genes_Counts
# table(Sig_Genes_Counts$table$FDR <= 0.1) # no significant genes
# 
# 
# # 3 negative samples might not be enough, edgeR recommends 6 samples at least
# #############################################################################
# # increasing the CMP cutoff and using the summed technical replicates
# y <- DGEList(counts = Counts_Pooled_protein_coding, group = clinical_info$Response10[2*(1:15)])
# keep_Counts <- filterByExpr(y)
# y <- y[keep_Counts, , keep.lib.sizes=FALSE]
# y <- calcNormFactors(y)
# y$samples
# dsgn <- model.matrix(~ clinical_info$Response10[2*(1:15)])
# rownames(dsgn) <- colnames(y)
# fit <- glmQLFit(y, dsgn)
# qlf <- glmQLFTest(fit)
# Sig_Genes_Counts <- topTags(qlf, n = 19718, adjust.method = "BH")
# # BY adjustment is valid under general dependence of multi-testing, produces larger adjusted pvalues
# Sig_Genes_Counts
# table(Sig_Genes_Counts$table$FDR <= 0.1) # no significant genes

# increasing the negative samples to 6 or even 7 patiensts does not give any significant genes


