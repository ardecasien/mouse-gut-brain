library(limma)
library(edgeR)
library(stringr)
library(umap)
library(variancePartition)
library(BiocParallel)
library(Rtsne)
library(reshape2)
library(ggrepel)
library(egg)
library(VennDiagram)
library(RColorBrewer)
library(corrplot)
library(psych)
library(ViSEAGO)
library(data.table)
library(parallel)
library(BRETIGEA)
library(biomaRt)
library(expss)
library(ggplot2)
library(dplyr)
library(tidyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(cowplot)

`%!in%` = Negate(`%in%`)

species.cols = c('#1b9e77', '#7570b3','#e6ab02')

#################
## load data
#################

fc = read.csv("FC_rawcounts.csv")
rownames(fc) = fc$X
fc = fc[,-1]

meta = read.csv("mouse_meta.csv")
meta.tech = read.csv("mouse_tech_meta.csv")
meta = merge(meta, meta.tech, by = 'ID', all = T)
meta$reads = as.numeric(meta$reads)
meta$scaled_mapping = scale(meta$mapping_rate)
meta$scaled_reads = scale(meta$reads)

meta.qc = read.csv('mouse_qc_meta.csv')
meta.qc <- meta.qc %>% mutate(ID = str_extract(Sample.Name, "(?<=-)[0-9]+(?=-)"))
meta.qc$ID = paste("X",meta.qc$ID,sep="")
meta.qc <- meta.qc %>%
  mutate(
    dups_num = as.numeric(str_remove(X..Dups, "%")),
    gc_num = as.numeric(str_remove(X..GC, "%")))
meta.qc = meta.qc %>% group_by(ID) %>% summarise(dups = mean(dups_num), gc = mean(gc_num))

meta = merge(meta, meta.qc, by = 'ID', all = T)
meta$scaled_dups = scale(meta$dups)
meta$scaled_gc = scale(meta$gc)

# remove outliers

rm.fc = c("X908", "X874", "X923")
fc = fc[,colnames(fc) %!in% rm.fc]

# load Zhu et al. data

hum_mac = readRDS('human_v_macaque_zhu_MFC.rds')

###################
## filter genes
###################

samples.by.species.and.batch = split(meta$ID, meta$Species.Batch, drop = TRUE)

cpm.cutoff = 10

keeps.fc = data.frame(gene = NULL)
for (i in 1:length(samples.by.species.and.batch)){
  keep.now = data.frame(gene = names(which(rowMeans(cpm(fc[,colnames(fc)%in%samples.by.species.and.batch[[i]]])) >= cpm.cutoff)))
  keeps.fc = rbind(keeps.fc, keep.now)}
keeps.fc = unique(keeps.fc)

datanow = fc[rownames(fc) %in% keeps.fc$gene,]
metanow = subset(meta, ID %in% colnames(fc))

mert = useMart(biomart = 'ENSEMBL_MART_ENSEMBL', dataset='mmusculus_gene_ensembl')
sub = getBM(attributes=c('ensembl_gene_id','external_gene_name','gene_biotype'), mart = mert)
colnames(sub)[1] = 'gene'
saveRDS(sub, file = 'sub.rds')

sub = readRDS('sub.rds')
keep.bt = merge(keeps.fc, sub, by = 'gene', all.x = T)
table(keep.bt$gene_biotype)
ribosomal_genes <- grep("^Rp[sl]", keep.bt$external_gene_name, value = TRUE)
ribosomal_ids = subset(keep.bt, external_gene_name %in% ribosomal_genes)$gene
ribosomal_genes

####################
## investigate technical variables (Figure S1)
####################

ribosomal_expr <- datanow[rownames(datanow) %in% ribosomal_ids, ]
total_expr <- colSums(datanow)
ribosomal_sum <- colSums(ribosomal_expr)
pct_ribosomal <- (ribosomal_sum / total_expr) * 100
metanow$ID == names(pct_ribosomal)
metanow$scaled_pct_ribosomal <- scale(pct_ribosomal)

boxplot(metanow$scaled_gc ~ metanow$Species)
boxplot(metanow$scaled_dups ~ metanow$Species)
boxplot(metanow$scaled_mapping ~ metanow$Species)
boxplot(metanow$scaled_reads ~ metanow$Species)
boxplot(metanow$scaled_reads ~ metanow$Species)
boxplot(metanow$scaled_pct_ribosomal ~ metanow$Species)

dgList <- DGEList(counts=datanow, genes=rownames(datanow))
dgList <- calcNormFactors(dgList, method = 'TMM')
y = voomWithQualityWeights(dgList, plot = T)

pca = prcomp(cor(y$E))
plot(pca$sdev^2 / sum(pca$sdev^2))
abline(h=0.01)
top_pcs <- pca$x[, which(pca$sdev^2 / sum(pca$sdev^2)>0.01)]
dim(top_pcs)[2]
all_vars <- c("scaled_mapping", "scaled_reads", "scaled_dups","scaled_gc","Batch", "Species", "scaled_pct_ribosomal")

c = cor(metanow[,c("scaled_mapping", "scaled_reads","scaled_pct_ribosomal","scaled_dups","scaled_gc")])
corrplot::corrplot(c)

results <- data.frame()

for (pc_num in c(1:dim(top_pcs)[2])) {
  pc_name <- paste0("PC", pc_num)
  pc_scores <- top_pcs[, pc_num]
  
  for (var in all_vars) {
    lm_fit <- lm(pc_scores ~ metanow[[var]])
    anova_res <- anova(lm_fit)
    p_val <- anova_res$`Pr(>F)`[1]
    f_stat <- anova_res$`F value`[1]
    
    results <- rbind(results, data.frame(
      PC = pc_name,
      Variable = var,
      Test = "ANOVA",
      Statistic = f_stat,
      P_value = p_val
    ))
  }
}

results$Adj_P_value <- p.adjust(results$P_value, method = "bonferroni")
View(results)

results$sig_label = ifelse(results$Adj_P_value < 0.05, "*", "")

results$PC = factor(results$PC, levels = paste("PC",c(1:10),sep=""))

ggplot(results, aes(x = PC, y = Variable, fill = Statistic)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  geom_text(aes(label = sig_label), color = "black", size = 6) +   # overlay significance symbols
  theme_article() +
  labs(title = "Heatmap of Statistic values by Variable and PC",
       fill = "Statistic") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

metanow = cbind(metanow, pca$x[,c(1,2,4)])

####################
## voom normalization
####################

dgList <- DGEList(counts=datanow, genes=rownames(datanow), group  = factor(metanow$Species.Batch))
dgList <- calcNormFactors(dgList, method = 'TMM')
ribosomal_genes <- rownames(dgList) %in% ribosomal_ids
dgList_filtered <- dgList[!ribosomal_genes, , keep.lib.sizes = FALSE]
y = voomWithQualityWeights(dgList_filtered, plot = T)

#####################
# variance partitioning (Figure 2A)
#####################

n.cores = detectCores() - 4
param = SnowParam(n.cores, "SOCK", progressbar=TRUE)
design = as.formula(paste('~',paste(c('PC1','PC2','PC4','(Species)'),collapse=' + ')))
varPart = fitExtractVarPartModel(y$E, design, metanow, BPPARAM=param)
vp = sortCols(varPart)
colnames(vp) 
colnames(vp) = c('Gut species','PC2','PC1','PC4','Residuals')
cols = c('#1b9e77', '#7570b3','#e6ab02','#e7298a')

mean(vp$`Gut species`) * 100
mean(vp$`PC2`) * 100
mean(vp$`PC1`) * 100
mean(vp$`PC4`) * 100
mean(vp$Residuals) * 100

median(vp$`Gut species`) * 100
median(vp$`PC2`) * 100
median(vp$`PC1`) * 100
median(vp$`PC4`) * 100
median(vp$Residuals) * 100

plotVarPart(vp, col = c(cols,'grey')) +
  geom_boxplot(width = 0.07, alpha=0.8, fill="grey") +
  theme_article() +
  ylab('% variance explained') +
  theme(axis.title = element_text(size=18), 
        axis.text = element_text(size=18),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") 

write.csv(vp, file = 'vp.csv')

######################
# plot expression data (Figures 3A, 3B, S3)
######################

# remove batch and covariate effects

nobatch = removeBatchEffect(y$E, covariates = cbind(metanow$PC1, metanow$PC2, metanow$PC4))

# MDS plot

species_colors <- c("Squirrel.Monkey" = "red","Macaque" = "blue","Human" = "green")
cols <- species_colors[metanow$Species]
plotMDS(nobatch, col = cols)

# PCA

pca = prcomp(cor(nobatch))
plot.pca = as.data.frame(pca$x)
plot.pca$Species = metanow$Species
plot.pca$Species[which(plot.pca$Species == 'Squirrel.Monkey')] = 'Squirrel Monkey'
plot.pca$Batch = factor(metanow$Batch)

centroids <- plot.pca %>%
  group_by(Species) %>%
  summarise(PC1 = mean(PC1), PC2 = mean(PC2))

ggplot(plot.pca,aes(x=PC1,y=PC2,color=Species, shape=Batch)) + 
  scale_color_manual(values=c(species.cols)) +
  geom_point(size=5) + 
  geom_point(data = centroids, aes(x = PC1, y = PC2, color=Species), 
             color = species.cols, shape = 21, alpha=0.4,size = 10,stroke=2) +
  theme_article() + 
  xlab('PC 1') + ylab('PC 2') + 
  theme(axis.ticks=element_blank(),
        axis.text=element_blank(),
        axis.title = element_text(size=18), 
        legend.text = element_text(size=14), 
        legend.title = element_blank())

by_species <- split(plot.pca[, c("PC1", "PC2")], plot.pca$Species)

library(purrr)
pairwise_species_dist <- combn(names(by_species), 2, simplify = FALSE) %>%
  map_df(function(species_pair) {
    sp1 <- species_pair[1]
    sp2 <- species_pair[2]
    dists <- as.matrix(dist(rbind(by_species[[sp1]], by_species[[sp2]])))
    n1 <- nrow(by_species[[sp1]])
    n2 <- nrow(by_species[[sp2]])
    between_distances <- dists[1:n1, (n1+1):(n1+n2)]
    tibble(
      Species1 = sp1,
      Species2 = sp2,
      MeanDistance = mean(between_distances)
    )
  })
View(pairwise_species_dist)

between_dist_df <- combn(names(by_species), 2, simplify = FALSE) %>%
  map_df(function(species_pair) {
    sp1 <- species_pair[1]
    sp2 <- species_pair[2]
    mat <- as.matrix(dist(rbind(by_species[[sp1]], by_species[[sp2]])))
    n1 <- nrow(by_species[[sp1]])
    n2 <- nrow(by_species[[sp2]])
    between_mat <- mat[1:n1, (n1 + 1):(n1 + n2)]
    
    expand.grid(Species1 = rownames(by_species[[sp1]]),
                Species2 = rownames(by_species[[sp2]])) %>%
      mutate(
        Pair = paste0(sp1, " vs ", sp2),
        Distance = as.vector(between_mat),
        SpeciesA = sp1,
        SpeciesB = sp2
      )
  })

ggplot(between_dist_df, aes(x = Pair, y = Distance)) +
  geom_boxplot(outlier.shape = NA, fill = "lightgray") +
  geom_jitter(width = 0.2, alpha = 0.3) +
  theme_minimal() +
  xlab("Species Pair") +
  ylab("Pairwise Distance") +
  ggtitle("Between-Species Pairwise Distances") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# tSNE

a = Rtsne(t(nobatch), dims = 2, perplexity=10, verbose=TRUE, max_iter = 1000)
plot.tsne = as.data.frame(a$Y)
plot.tsne$Species = metanow$Species
plot.tsne$Species[which(plot.tsne$Species == 'Squirrel.Monkey')] = 'Squirrel Monkey'
plot.tsne$Batch = factor(metanow$Batch)

centroids <- plot.tsne %>%
  group_by(Species) %>%
  summarise(V1 = mean(V1), V2 = mean(V2))

ggplot(plot.tsne,aes(x=V1,y=V2,color=Species, shape=Batch)) + 
  scale_color_manual(values=c(species.cols)) +
  geom_point(size=5) + 
  geom_point(data = centroids, aes(x = V1, y = V2, color=Species), 
             color = species.cols, shape = 21, alpha=0.4,size = 10,stroke=2) + theme_article() + 
  xlab('V1') + ylab('V2') + 
  theme(axis.ticks=element_blank(),
        axis.text=element_blank(),
        axis.title = element_text(size=18), 
        legend.text = element_text(size=14), 
        legend.title = element_blank())

# UMAP

set.seed(102)
a = umap(t(nobatch), n_neighbors = 30, min_dist = 0.9, metric = 'manhattan') 

a.umap = data.frame(as.data.frame(a$layout))
a.umap$Species = metanow$Species
a.umap$Species[which(a.umap$Species == 'Squirrel.Monkey')] = 'Squirrel Monkey'
a.umap$Batch = metanow$Batch

centroids <- a.umap %>%
  group_by(Species) %>%
  summarise(V1 = mean(V1), V2 = mean(V2))

ggplot(a.umap,aes(x=V1,y=V2,color=Species,shape =Batch)) + 
  scale_color_manual(values=c(species.cols)) +
  geom_point(size=5) + 
  geom_point(data = centroids, aes(x = V1, y = V2, color=Species), 
             color = species.cols, shape = 21, alpha=0.4,size = 10,stroke=2) + theme_article() + 
  theme_article() +
  xlab('UMAP 1') + ylab('UMAP 2') + 
  theme(axis.ticks=element_blank(),
        axis.text=element_blank(),
        axis.title = element_text(size=18), 
        legend.text = element_text(size=14), 
        legend.title = element_blank())

###########
# run model
###########

metanow$Species.Batch = as.factor(metanow$Species.Batch)
levels(metanow$Species.Batch)
design <- model.matrix(~0 + PC1 + PC2 + PC4 + Species.Batch, data = metanow)
contrasts <- makeContrasts(
  HvM=(Species.BatchHuman.A+Species.BatchHuman.B)/2-(Species.BatchMacaque.E+Species.BatchMacaque.F)/2,
  HvS=(Species.BatchHuman.A+Species.BatchHuman.B)/2-(Species.BatchSquirrel.Monkey.C+Species.BatchSquirrel.Monkey.D)/2,
  SvM=(Species.BatchSquirrel.Monkey.C+Species.BatchSquirrel.Monkey.D)/2-(Species.BatchMacaque.E+Species.BatchMacaque.F)/2,
  HvSM=(Species.BatchHuman.A+Species.BatchHuman.B)/2-(Species.BatchSquirrel.Monkey.C+Species.BatchSquirrel.Monkey.D+Species.BatchMacaque.E+Species.BatchMacaque.F)/4,
  HSvM=(Species.BatchHuman.A+Species.BatchHuman.B+Species.BatchSquirrel.Monkey.C+Species.BatchSquirrel.Monkey.D)/4-(Species.BatchMacaque.E+Species.BatchMacaque.F)/2,
  AvY=(Species.BatchHuman.A+Species.BatchSquirrel.Monkey.C+Species.BatchMacaque.E)/3-(Species.BatchHuman.B+Species.BatchMacaque.F+Species.BatchSquirrel.Monkey.D)/3,
  HAvHY=Species.BatchHuman.A-Species.BatchHuman.B,
  SAvSY=Species.BatchSquirrel.Monkey.C-Species.BatchSquirrel.Monkey.D,
  MAvMY=Species.BatchMacaque.E-Species.BatchMacaque.F,
  HAvALL=Species.BatchHuman.A-(Species.BatchHuman.B+Species.BatchSquirrel.Monkey.C+Species.BatchSquirrel.Monkey.D+Species.BatchMacaque.E+Species.BatchMacaque.F)/5,
  HYvALL=Species.BatchHuman.B-(Species.BatchHuman.A+Species.BatchSquirrel.Monkey.C+Species.BatchSquirrel.Monkey.D+Species.BatchMacaque.E+Species.BatchMacaque.F)/5,
  SAvALL=Species.BatchSquirrel.Monkey.C-(Species.BatchHuman.B+Species.BatchHuman.A+Species.BatchSquirrel.Monkey.D+Species.BatchMacaque.E+Species.BatchMacaque.F)/5,
  SYvALL=Species.BatchSquirrel.Monkey.D-(Species.BatchHuman.B+Species.BatchSquirrel.Monkey.C+Species.BatchHuman.A+Species.BatchMacaque.E+Species.BatchMacaque.F)/5,
  MAvALL=Species.BatchMacaque.E-(Species.BatchHuman.B+Species.BatchSquirrel.Monkey.C+Species.BatchSquirrel.Monkey.D+Species.BatchHuman.A+Species.BatchMacaque.F)/5,
  MYvALL=Species.BatchMacaque.F-(Species.BatchHuman.B+Species.BatchSquirrel.Monkey.C+Species.BatchSquirrel.Monkey.D+Species.BatchMacaque.E+Species.BatchHuman.A)/5,
  levels=colnames(design))

rownames(contrasts) = gsub("Species.Batch", "", rownames(contrasts))
contrasts
round(rowSums(t(contrasts)),3)
colnames(design) = gsub("Species.Batch", "", colnames(design))
head(design)

fit <- lmFit(y, design)
fit <- contrasts.fit(fit, contrasts)
fit <- eBayes(fit)
coef = fit$coefficients
pvals = fit$p.value
stdev = fit$stdev.unscaled
tstat = fit$t

pvals.adj = pvals
for (i in 1:length(colnames(pvals))){
  pvals.adj[,i] = p.adjust(pvals[,i], method = 'BH')
}

coef2 = reshape2::melt(coef)
stdev2 = reshape2::melt(stdev)
colnames(stdev2)[2] = 'Contrasts'
combo = merge(coef2, stdev2, by = c('Var1','Contrasts'))
combo = merge(combo, reshape2::melt(pvals), by = c('Var1','Contrasts'))
combo = merge(combo, reshape2::melt(pvals.adj), by = c('Var1','Contrasts'))
colnames(combo) = c('gene','contrast','beta','stdev','p','padj')
combo = merge(combo, sub, by = 'gene', all.x = T)
saveRDS(combo, 'combo_results.rds')

write.csv(combo, file = 'de-results.csv')

###########
# permutations (Figure S3)
###########

n_perm <- 1000
set.seed(123)  
perm_tstats <- array(NA, dim = c(nrow(y$E), ncol(tstat), n_perm),dimnames = list(rownames(y$E), colnames(tstat), NULL))
covariates <- metanow[, c("PC1", "PC2", "PC4")]
for (perm in 1:n_perm) {
  shuffled_species <- sample(metanow$Species.Batch)
  perm_data <- metanow
  perm_data$Species.Batch <- shuffled_species
  design_perm <- model.matrix(~0 + PC1 + PC2 + PC4 + Species.Batch, data = perm_data)
  colnames(design_perm) = gsub("Species.Batch", "", colnames(design))
  fit_perm <- lmFit(y, design_perm)
  fit_perm <- contrasts.fit(fit_perm, contrasts)
  fit_perm <- eBayes(fit_perm)  
  perm_tstats[,,perm] <- fit_perm$t
}

saveRDS(perm_tstats, 'perm_tstats.rds')

perm_tstats = readRDS('perm_tstats.rds')

empirical_pvals <- matrix(NA, nrow = nrow(y$E), ncol = ncol(tstat),dimnames = list(rownames(y$E), colnames(tstat)))

for (contrast_i in 1:ncol(tstat)) {
  for (gene_i in 1:nrow(y$E)) {
    obs_t <- abs(tstat[gene_i, contrast_i])
    perm_t <- abs(perm_tstats[gene_i, contrast_i, ])
    empirical_pvals[gene_i, contrast_i] <- sum(perm_t >= abs(obs_t)) / length(perm_t)
  }
}
empirical_pvals_adj <- apply(empirical_pvals, 2, p.adjust, method = "BH")

coef_orig <- fit$coefficients
combo_perm <- data.frame(
  gene = rep(rownames(y$E), times = ncol(coef_orig)),
  contrast = rep(colnames(coef_orig), each = nrow(y$E)),
  beta = as.vector(coef_orig),
  empirical_p = as.vector(empirical_pvals),
  empirical_padj = as.vector(empirical_pvals_adj)
)

table(combo_perm$contrast, combo_perm$empirical_padj < 0.2)

saveRDS(combo_perm, 'combo_results_permutation.rds')

combo_perm = readRDS('combo_results_permutation.rds')
write.csv(combo_perm, 'combo_perm.csv')

combo = read.csv('de-results.csv')
combo = subset(combo, contrast %in% c('HvM', 'SvM','HvS'))
combo_perm = read.csv('combo_perm.csv')
combo_perm = subset(combo_perm, contrast %in% c('HvM', 'SvM','HvS'))
pl = merge(combo, combo_perm, by = c('gene','contrast'))

ggplot(pl, aes(x = p, y = empirical_p)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  theme_article() +
  xlab('voom-limma p-values') +
  ylab('permuted p-values')

################################
# volcano plots (Figure 2B-D)
################################

combo = read.csv('de-results.csv')

pl = subset(combo, contrast == 'HvM')
pl$code = ifelse(pl$padj < 0.2 & pl$beta > 0, 'human', ifelse(pl$padj < 0.2 & pl$beta < 0, 'macaque', 'none'))
colsnow = c(species.cols[c(1,2)],'#666666')

pl = subset(combo, contrast == 'HvS')
pl$code = ifelse(pl$padj < 0.2 & pl$beta > 0, 'human', ifelse(pl$padj < 0.2 & pl$beta < 0, 'sm', 'none'))
colsnow = c(species.cols[1],'#666666',species.cols[3])

pl = subset(combo, contrast == 'SvM')
pl$code = ifelse(pl$padj < 0.2 & pl$beta > 0, 'sq', ifelse(pl$padj < 0.2 & pl$beta < 0, 'macaque', 'none'))
colsnow = c(species.cols[2],'#666666',species.cols[3])

table(pl$code)
ggplot(pl, aes(x = beta, y = -log10(padj), color = code)) +
  geom_point(alpha = 0.5) + 
  scale_color_manual(values = colsnow) +
  xlim(c(-2.5,2.5)) +
  ylim(c(0,6)) +
  ylab('-log10(padj)') +
  xlab(expression(~italic(beta))) +
  geom_hline(yintercept = -log10(0.2), linetype = 'dashed') +
  geom_text_repel(data=subset(pl, padj < 0.2),
                  aes(x = beta, y = -log10(padj),label=external_gene_name), 
                  color = 'black', size = 3, max.overlaps = 20) +
  theme_article() +
  theme(legend.position = 'none')

################################
# corr plot (Figure 3D)
################################

comp = coef[,c('HvM','HvS','SvM')]
corr.test(comp)
m = cor(comp)
colnames(m) = rownames(m) = c('HIMs v MIMs','HIMs v SMIMs','SMIMs v MIMs')
COL2(diverging = c("RdBu", "BrBG", "PiYG", "PRGn", "PuOr", "RdYlBu"), n = 200)
corrplot.mixed(m, 
               lower.col = COL2('BrBG'), upper.col = COL2('BrBG'), 
               tl.cex = 1.2, cl.cex = 1, 
               tl.col = 'black')

################
## venn diagrams (Figure 3C)
################

myCol <- brewer.pal(3, "Dark2")

mouse_res = read.csv('de-results.csv')
#mouse_res = read.csv('combo_perm.csv')
#colnames(mouse_res)[6] = 'padj'

hvm = subset(mouse_res, contrast == 'HvM' & padj < 0.2)$gene
hvs = subset(mouse_res, contrast == 'HvS' & padj < 0.2)$gene
svm = subset(mouse_res, contrast == 'SvM' & padj < 0.2)$gene

length(hvm)
table(hvm %in% svm)
table(hvm %in% hvs)

venn.diagram(x = list(hvm, hvs, svm),
             category.names = c("HIMs v MIMs" , "HIMs v SMIMs" , "SMIMs v MIMs"),
             file = 'species_venn_perm.png',
             lwd = 2,
             lty = 'blank',
             fill = myCol,
             fontfamily = "sans",
             cat.cex = 1,
             margin = 0.05,
             cat.fontface = "bold",
             cat.fontfamily = "sans")

################################
# plot effect size distributions (Figure 3C)
################################

coef.plot = reshape2::melt(coef)

ggplot(coef.plot, aes(x = abs(value))) + 
  theme_article() +
  xlab('Absolute Beta') +
  ylab('Count') +
  geom_histogram(data=subset(coef.plot,Contrasts == 'HvS'),fill = myCol[2], alpha = 0.4, bins = 50) +
  geom_histogram(data=subset(coef.plot,Contrasts == 'HvM'),fill = myCol[1], alpha = 0.4, bins = 50) +
  geom_histogram(data=subset(coef.plot,Contrasts == 'SvM'),fill = myCol[3], alpha = 0.4, bins = 50) +
  theme(axis.title = element_text(size = 18)) +
  theme(axis.text = element_text(size = 18)) +
  scale_x_continuous(expand = c(0, 0), limits = c(1e-100, .99)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA))

coef.plot = subset(coef.plot, Contrasts == 'HvS' | Contrasts == 'HvM' | Contrasts == 'SvM')
coef.plot %>% group_by(Contrasts) %>% summarise(mean = mean(abs(value)))
mod = aov(abs(value) ~ Contrasts, data = coef.plot)
summary(mod)
TukeyHSD(mod)

################################
## mimenet (Figure 5)
################################

mouse_res = read.csv('de-results.csv')
mouse_res = subset(mouse_res, contrast == 'HvM')
dim(mouse_res)

de_now = subset(mouse_res, padj < 0.2 & beta > 0)
#de_now = subset(mouse_res, padj < 0.2 & beta < 0)

gene_clusters = read.csv('metabolite_clusters.csv')
dim(gene_clusters)
gene_clusters$X[which(gene_clusters$X %!in% mouse_res$gene)]
gene_clusters = subset(gene_clusters, X %in% mouse_res$gene)
cl = unique(gene_clusters$Cluster)

o = data.frame()
for(i in 1:length(cl)) {
  
  clnow = subset(gene_clusters, Cluster == cl[i])
  dim(clnow)
  
  b = length(subset(clnow, X %in% de_now$gene)$X)
  d = length(subset(de_now, gene %!in% clnow$X)$gene)
  c = length(subset(clnow, X %!in% de_now$gene)$X)
  n = length(subset(mouse_res, gene %!in% c(clnow$X, de_now$gene))$gene)
  m = matrix(c(b,d,c,n), nrow=2, ncol=2)
  m
  f = fisher.test(m, alternative = 'greater')
  o[i,1] = cl[i]
  o[i,2] = f$estimate
  o[i,3] = f$p.value
}
colnames(o) = c('cluster','OR','p')
o$adj = p.adjust(o$p)
View(o)
# OR > 1 and padj < 0.05 for clusters 1,3,0,6

o$cluster = as.character(o$cluster)
ggplot(o, aes(x = cluster, y = OR)) +
  geom_bar(stat = 'identity', fill = '#E6AB02') +
  xlab('Gene cluster') +
  ylab('Enrichments for HIMs > MIMs DE genes\n(odds ratios)') +
  theme_article() 

micro_clusters  = read.csv('microbe_clusters.csv')
table(micro_clusters$Cluster)
sc = read.csv('interaction_score_matrix.csv')
rownames(sc) = sc$X..Pathway
sc = sc[,-1]

gc = unique(gene_clusters$Cluster)
mc = unique(micro_clusters$Cluster)

cl_of_interest = c(0,1,3,6)
u = data.frame()
g = data.frame()
for(i in 1:length(cl_of_interest)) {
  
  genesnow = subset(gene_clusters, Cluster == cl_of_interest[i])$X
  
  for(j in 1:length(mc)){
    
    micronow = c(subset(micro_clusters, Cluster == mc[j])$X..Pathway)
    scnow = sc[micronow, genesnow]
    scmelt = reshape2::melt(scnow)
    me = mean(scmelt$value, na.rm = T)
    su = sum(scmelt$value > 0, na.rm = T)/length(scmelt$value)
    u[i,j] = me
    g[i,j] = su
  }
  
}
u = t(u)
g = t(g)
colnames(u) = colnames(g) = cl_of_interest
rownames(u) = rownames(g) = mc
u = data.frame(u)
g = data.frame(g)
u$mean = rowMeans(u)
g$mean = rowMeans(g)
View(u)
View(g)

# microbe clusters 4, 1 have highest mean interaction score with gene clusters 0 1 3 6

View(subset(micro_clusters, Cluster == 4)) # glucose degredation
View(subset(micro_clusters, Cluster == 1)) # gluconeogensis

pl = u
pl$cluster=rownames(pl)
pl = pl[,c(5,6)]
ggplot(pl, aes(x=cluster, y = 'mean',fill = mean)) +
  geom_tile() +
  scale_fill_gradient2(low = '#666666', high = '#e6ab02', mid = 'white', midpoint = 0) +
  theme_article()

################################
## compare to primate data (Figures 2E 2F S2)
################################

mouse_res = read.csv('de-results.csv')
mouse_res = subset(mouse_res, contrast == 'HvM')
colnames(mouse_res)[2] = 'mmusculus_homolog_ensembl_gene'

mert = useMart(biomart = 'ENSEMBL_MART_ENSEMBL', dataset='hsapiens_gene_ensembl')
conv = getBM(attributes=c('ensembl_gene_id','external_gene_name','mmusculus_homolog_ensembl_gene','mmusculus_homolog_orthology_type'), mart = mert)
conv_121 = subset(conv, mmusculus_homolog_orthology_type == 'ortholog_one2one')
saveRDS(conv_121, file = 'conv_121.rds')

conv_121 = readRDS('conv_121.rds')
mouse_res = merge(mouse_res, conv_121, by = 'mmusculus_homolog_ensembl_gene', all.x = T)

prim_res = hum_mac
prim_res$ensembl_gene_id = rownames(prim_res)

combo = merge(prim_res, mouse_res, by = 'ensembl_gene_id')
combo$code = ifelse(combo$beta.x > 0 & combo$beta.y > 0, 'concordant', 'discordant')
combo$code = ifelse(combo$beta.x < 0 & combo$beta.y < 0, 'concordant', combo$code)
length(combo$ensembl_gene_id)
table(combo$code)
table(combo$code)[1] / sum(table(combo$code))
cor.test(x = combo$beta.x, y = combo$beta.y, method = 'spearman')

cut = data.frame()
p <- round(seq(from = 0.01, to = 1, by = 0.01), digits = 4)
for(i in 1:length(p)) {
  threshold <- round(as.numeric(p[i]), 2)
  #s = subset(combo, padj.x < threshold | padj.y < threshold)
  s = subset(combo, padj.x < threshold & padj.y < threshold)
  if (nrow(s) < 3) {
    cut[i, 1] <- threshold
    cut[i, 2:8] <- NA
  } else {
  sp = cor.test(x = s$beta.x, y = s$beta.y, method = 'spearman')
  pe = cor.test(x = s$beta.x, y = s$beta.y, method = 'pearson')
  cut[i,1] = p[i]
  cut[i,2] = sp$estimate
  cut[i,3] = sp$p.value
  cut[i,4] = pe$estimate
  cut[i,5] = pe$p.value
  upup = length(subset(s, beta.x > 0 & s$beta.y > 0)$external_gene_name.x)
  downdown = length(subset(s, beta.x < 0 & s$beta.y < 0)$external_gene_name.x)
  updown = length(subset(s, beta.x > 0 & s$beta.y < 0)$external_gene_name.x)
  downup = length(subset(s, beta.x < 0 & s$beta.y > 0)$external_gene_name.x)
  m = matrix(ncol=2,nrow=2,c(upup, updown, downup, downdown))
  m
  f = fisher.test(m, alternative = 'greater')
  cut[i,6] = f$estimate
  cut[i,7] = f$p.value
  cut[i,8] = nrow(s)
  }
}
View(cut)
#write.csv(cut, file = 'cutoffs-either.csv')
write.csv(cut, file = 'cutoffs-both.csv')

cut$code = ifelse(cut$V3<0.05, 'sig', 'ns')
ggplot(cut, aes(x = V1, y = V2, color = code)) + 
  geom_point() +
  scale_color_manual(values = c('grey','darkblue')) +
  #geom_text(data = subset(cut, V3 < 0.05), aes(label = "*"), vjust = -1.2, size = 5) +
  labs(x = "FDR cutoff", y = "Spearman correlation") +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  theme_article()
cut$code = ifelse(cut$V7<0.05, 'sig', 'ns')
ggplot(cut, aes(x = V1, y = V6, color = code)) + 
  geom_point() +
  scale_color_manual(values = c('grey','darkblue')) +
  labs(x = "FDR cutoff", y = "Odds ratio") +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  theme_article()

subset(cut, V3<0.05 & V8>10)
pnow = cut$V1[which(cut$V1 == min(subset(cut, V3<0.05 & V8>10)$V1))]
pnow
#combo3 = subset(combo, padj.x < pnow | padj.y < pnow)
combo3 = subset(combo, padj.x < pnow & padj.y < pnow)
dim(combo3)
combo3$code = ifelse(combo3$beta.x > 0 & combo3$beta.y > 0, 'Concordant', 'Discordant')
combo3$code = ifelse(combo3$beta.x < 0 & combo3$beta.y < 0, 'Concordant', combo3$code)
table(combo3$code)
table(combo3$code)[1] / sum(table(combo3$code))

cor.test(x = combo3$beta.x, y = combo3$beta.y, method = 'spearman')
upup = length(subset(combo3, beta.x > 0 & combo3$beta.y > 0)$external_gene_name.x)
downdown = length(subset(combo3, beta.x < 0 & combo3$beta.y < 0)$external_gene_name.x)
updown = length(subset(combo3, beta.x > 0 & combo3$beta.y < 0)$external_gene_name.x)
downup = length(subset(combo3, beta.x < 0 & combo3$beta.y > 0)$external_gene_name.x)
m = matrix(ncol=2,nrow=2,c(upup, updown, downup, downdown))
m
fisher.test(m, alternative = 'greater')

ggplot(combo3, aes(x = beta.x, y = beta.y)) +
  geom_point(size=2.5, aes(x = beta.x, y = beta.y, color = code)) +
  geom_smooth(method = 'lm') +
  theme_article() +
  xlab('Human vs Macaque Brain Tissue') +
  ylab('HIMs vs MIMs Brain Tissue') +
  scale_color_manual(values = c('#BC9CB0','#DDF2EB')) +
  theme(legend.position = 'bottom',
        legend.title = element_blank(),
        legend.text = element_text(size=18),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 18))

h_b = subset(combo3, beta.x > 0)$beta.x
h_b = data.frame(h_b)
colnames(h_b) = 'beta'
h_b$group = 'Human Upregulated'
h_b$group2 = 'Primate Brain'
m_b = subset(combo3, beta.x < 0)$beta.x
m_b = data.frame(m_b)
colnames(m_b) = 'beta'
m_b$group = 'Macaque Upregulated'
m_b$group2 = 'Primate Brain'
h_m = subset(combo3, beta.y > 0)$beta.y
h_m = data.frame(h_m)
colnames(h_m) = 'beta'
h_m$group = 'Human Upregulated'
h_m$group2 = 'Mouse Brain'
m_m = subset(combo3, beta.y < 0)$beta.y
m_m = data.frame(m_m)
colnames(m_m) = 'beta'
m_m$group = 'Macaque Upregulated'
m_m$group2 = 'Mouse Brain'

tt = rbind(h_b, (m_b), (h_m), (m_m))
ggplot(tt, aes(x = group, y=abs(beta), fill = group2)) +
  facet_wrap(~group, scales = 'free_x') +
  geom_violin(alpha = 0.5) +
  theme_article() +
  ylab('Absolute Beta') +
  scale_fill_manual(values = c('#1b9e77', '#7570b3')) +
  theme(legend.position = 'bottom',
        legend.title = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),       
        axis.title.x = element_blank(),
        legend.text = element_text(size=18),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 18))

tt$code = paste(tt$group, tt$group2)
mod = aov(abs(beta) ~ code, data = tt)
summary(mod)
TukeyHSD(mod)

t.test(abs(beta) ~ group2, data = subset(tt, group == 'Human Upregulated'))
t.test(abs(beta) ~ group2, data = subset(tt, group == 'Macaque Upregulated'))

############################
## test for overlap with positive/purifying selection
############################

mouse_res = read.csv('de-results.csv')
mouse_res = subset(mouse_res, contrast == 'HvM')
colnames(mouse_res)[2] = 'mmusculus_homolog_ensembl_gene'
mert = useMart(biomart = 'ENSEMBL_MART_ENSEMBL', dataset='hsapiens_gene_ensembl')
conv = getBM(attributes=c('ensembl_gene_id','external_gene_name','mmusculus_homolog_ensembl_gene','mmusculus_homolog_orthology_type'), mart = mert)
conv_121 = subset(conv, mmusculus_homolog_orthology_type == 'ortholog_one2one')
mouse_res = merge(mouse_res, conv_121, by = 'mmusculus_homolog_ensembl_gene', all.x = T)

dumas = read.csv('dumas.csv')
colnames(dumas)[2] = 'external_gene_name.y'
dumas$mean = rowMeans(dumas[,c(3:5)], na.rm = T)
combo = merge(mouse_res, dumas, by = 'external_gene_name.y')
datanow = combo[complete.cases(combo$ensembl_gene_id),]

ggplot(datanow, aes(x = HS)) + geom_histogram() + theme_article()
ggplot(datanow, aes(x = abs(HS))) + geom_histogram() + theme_article()
ct = mad(abs(datanow$HS))
ct
table(abs(datanow$HS) > ct)
table(abs(datanow$HS) > ct) / length(datanow$HS)
ggplot(datanow, aes(x = HS)) + geom_histogram() + geom_vline(xintercept = c(ct, -ct), linetype = 'dashed') + theme_article()

datanow = subset(datanow, abs(HS) > ct)
datanow = subset(datanow, padj < 0.2)

table(datanow$beta < 0, datanow$HS < 0)
m = matrix(table(datanow$beta < 0, datanow$HS < 0), nrow = 2, ncol = 2)
rownames(m) = c('HIMs','MIMs')
colnames(m) = c('positive','purify')
m
sum(m)
fisher.test(m, alternative='greater')
fisher.test(m)

View(subset(datanow, beta < 0 & HS < 0))
View(subset(datanow, beta > 0 & HS > 0))

############################
## biological processes (Figures 3A S4)
############################

mouse_res = read.csv('de-results.csv')
colnames(mouse_res)[2] = 'mmusculus_homolog_ensembl_gene'

#mert = useMart(biomart = 'ENSEMBL_MART_ENSEMBL', dataset='hsapiens_gene_ensembl')
#conv = getBM(attributes=c('ensembl_gene_id','external_gene_name','mmusculus_homolog_ensembl_gene','mmusculus_homolog_orthology_type'), mart = mert)
#conv_121 = subset(conv, mmusculus_homolog_orthology_type == 'ortholog_one2one')

conv_121 = readRDS('conv_121.rds')
mouse_res = merge(mouse_res, conv_121, by = 'mmusculus_homolog_ensembl_gene', all.x = T)

#contrastnow = 'HvM'
#contrastnow = 'HvS'
contrastnow = 'SvM' 

datanow = subset(mouse_res, contrast == contrastnow)
datanow = datanow[complete.cases(datanow$ensembl_gene_id),]

up_genes <- subset(datanow, padj < 0.2 & beta > 0)$ensembl_gene_id
down_genes <- subset(datanow, padj < 0.2 & beta < 0)$ensembl_gene_id
up_ids <- bitr(up_genes, fromType="ENSEMBL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
down_ids <- bitr(down_genes, fromType="ENSEMBL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
all_ids <- bitr(unique(datanow$ensembl_gene_id), fromType="ENSEMBL", toType="ENTREZID", OrgDb=org.Hs.eg.db)

go_up <- enrichGO(gene = up_ids$ENTREZID,
                  universe = all_ids$ENTREZID,
                  OrgDb = org.Hs.eg.db,
                  ont = "BP",
                  pAdjustMethod = "none",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 1,
                  readable = TRUE)
go_down <- enrichGO(gene = down_ids$ENTREZID,
                    universe = all_ids$ENTREZID,
                    OrgDb = org.Hs.eg.db,
                    ont = "BP",
                    pAdjustMethod = "none",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 1,
                    readable = TRUE)

# plot

convertGeneRatio <- function(gr) {
  sapply(strsplit(gr, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))}

df_up = as.data.frame(go_up)
df_down = as.data.frame(go_down)
View(df_up)
View(df_down)

df_up$GeneRatioNum <- convertGeneRatio(df_up$GeneRatio)
df_down$GeneRatioNum <- convertGeneRatio(df_down$GeneRatio)

df_up <- df_up %>% arrange(desc(GeneRatioNum))
df_down <- df_down %>% arrange(desc(GeneRatioNum))

combined_data <- rbind(df_up[c(1:10),], df_down[c(1:10),])
color_range <- range(-log10(combined_data$pvalue))
color_range[1] = -log10(0.05)
size_range <- range(combined_data$Count)
rat_range = range(combined_data$GeneRatioNum)

p1 <- ggplot(df_up[c(1:10),], aes(x = GeneRatioNum,
                        y = reorder(Description, GeneRatioNum),
                        size = Count,
                        color = -log10(p.adjust))) +
  geom_point() +
  xlim(rat_range) +
  scale_color_continuous(name = "-log10(p.adjust)", limits = color_range, low = "grey", high = "darkred") +
  scale_size(name = "Gene count", limits = size_range, range = c(3, 8)) +
  theme_article() +
  labs(title = "GO BP Enrichment: Upregulated", x = "GeneRatio", y = NULL)

p2 <- ggplot(df_down[c(1:10),], aes(x = GeneRatioNum,
                                  y = reorder(Description, GeneRatioNum),
                                  size = Count,
                                  color = -log10(p.adjust))) +
  geom_point() +
  xlim(rat_range) +
  scale_color_continuous(name = "-log10(p.adjust)", limits = color_range, low = "grey", high = "darkred") +
  scale_size(name = "Gene count", limits = size_range, range = c(3, 8)) +
  theme_article() +
  labs(title = "GO BP Enrichment: Downregulated", x = "GeneRatio", y = NULL)

(p1 / p2) & theme(legend.position = "right")

write.csv(as.data.frame(go_up), file = paste(contrastnow, '_GO_BP_up.csv', sep = ""), row.names = FALSE)
write.csv(as.data.frame(go_down), file = paste(contrastnow, '_GO_BP_down.csv', sep = ""), row.names = FALSE)

hvm = read.csv(paste('HvM', '_GO_BP_up.csv', sep = ""))
svm = read.csv(paste('SvM', '_GO_BP_up.csv', sep = ""))

table(hvm$ID %in% svm$ID)
hvm$Description[hvm$ID %in% svm$ID]

hvm = read.csv(paste('HvM', '_GO_BP_up.csv', sep = ""))
mvh = read.csv(paste('vM', '_GO_BP_down.csv', sep = ""))

table(hvm$ID %in% mvh$ID)
hvm$Description[hvm$ID %in% mvh$ID]

#####################
## disease risk gene enrichments (Figure 4C)
#####################

# Import disease associations from DISEASES dataset (i.e., "Disease Ontology")

do.data = read.table('human_disease_associations.tsv',
                     sep='\t',
                     quote='',
                     col.names=c('protein_id','protein_name','do_id','do_name','z_score','confidence'),
                     stringsAsFactors=FALSE)

do.def = unique(subset(do.data,select=c('do_id','do_name')))

ignore.checkpoints = FALSE
if (ignore.checkpoints || !file.exists('disease_human_orthologs.rds')) {
  library(biomaRt)
  
  hsap = useMart(biomart = 'ENSEMBL_MART_ENSEMBL',
                 dataset='hsapiens_gene_ensembl')
  
  hsap.info = getBM(
    attributes=c('ensembl_gene_id','ensembl_peptide_id','external_gene_name'),
    mart = hsap)
  
  saveRDS(hsap.info,file='disease_human_orthologs.rds')
} else {
  message('Checkpoint found!\nLoading human ortholog annotations from file.')
  
  hsap.info = readRDS('disease_human_orthologs.rds')
}

# For matching ensembl peptides, linking genes is straightforward
do.ensembl = subset(do.data,grepl('^ENSP[0-9]{11}',protein_id))
do.ensembl = merge(do.ensembl,hsap.info,by.x='protein_id','ensembl_peptide_id',all.x=FALSE,all.y=FALSE)
do.ensembl$ensembl_peptide_id = do.ensembl$protein_id

# For everything else
do.proteinname = subset(do.data,!grepl('^ENSP[0-9]{11}',protein_id))
do.proteinname = merge(do.proteinname,hsap.info,by.x='protein_name',by.y='external_gene_name',all.x=FALSE,all.y=FALSE)
do.proteinname$external_gene_name = do.proteinname$protein_name

do.all = rbind(do.ensembl[intersect(names(do.ensembl),names(do.proteinname))],do.proteinname[intersect(names(do.ensembl),names(do.proteinname))])

do.mmus = subset(do.all,select=c('external_gene_name','ensembl_gene_id','ensembl_peptide_id','do_id','do_name','z_score','confidence'))

# run after GO

allgenes = unique(datanow$ensembl_gene_id)
all.res.fet = all.res.kst = numeric(length=length(allgenes))
names(all.res.fet) = names(all.res.kst) = allgenes

up = unique(subset(datanow, padj < 0.2 & beta > 0)$ensembl_gene_id)
length(up)
down = unique(subset(datanow, padj < 0.2 & beta < 0)$ensembl_gene_id)
length(down)

all.res.fet[which(names(all.res.fet) %in% up)] = 1
all.res.fet[which(names(all.res.fet) %in% down)] = -1
table(all.res.fet)
names(all.res.fet) = allgenes

all.res.kst = datanow$beta / datanow$stdev

all.res.join = data.frame(ensembl_gene_id = allgenes, direction = as.integer(all.res.fet), effect = as.numeric(all.res.kst))
all.res.do = merge(all.res.join, do.mmus, by='ensembl_gene_id')

all.res.do.pass = subset(all.res.do,do_id %in% names(which(table(subset(all.res.do,confidence >= 0)$do_id) >= 10)))
all.res.do.gene.pass = unique(all.res.do.pass[c('ensembl_gene_id','direction')])

do.deg.total = as.integer(table(factor(all.res.do.pass$direction != 0,levels=c('TRUE','FALSE'))))
do.inc.total = as.integer(table(factor(all.res.do.pass$direction == 1 | all.res.do.pass$direction == 2,levels=c('TRUE','FALSE'))))
do.dec.total = as.integer(table(factor(all.res.do.pass$direction == -1 | all.res.do.pass$direction == 2,levels=c('TRUE','FALSE'))))

all.res.do.split = split(all.res.do.pass,all.res.do.pass$do_id)

n.cores = detectCores() - 6

all.res.do.test = do.call(rbind,mclapply(names(all.res.do.split),function(i) {
  x = all.res.do.split[[i]]
  
  this.deg.total = as.integer(table(factor(x$direction != 0,levels=c('TRUE','FALSE'))))
  contingency.matrix.deg = matrix(rbind(this.deg.total,do.deg.total - this.deg.total),nrow=2,dimnames=list(c('in group','not in group'),c('in direction','not in direction')))
  
  this.inc.total = as.integer(table(factor(x$direction == 1 | x$direction == 2,levels=c('TRUE','FALSE'))))
  contingency.matrix.inc = matrix(rbind(this.inc.total,do.inc.total - this.inc.total),nrow=2,dimnames=list(c('in group','not in group'),c('in direction','not in direction')))
  
  this.dec.total = as.integer(table(factor(x$direction == -1 | x$direction == 2,levels=c('TRUE','FALSE'))))
  contingency.matrix.dec = matrix(rbind(this.dec.total,do.dec.total - this.dec.total),nrow=2,dimnames=list(c('in group','not in group'),c('in direction','not in direction')))
  
  deg.fet.test = fisher.test(contingency.matrix.deg,alternative='greater')
  inc.fet.test = fisher.test(contingency.matrix.inc,alternative='greater')
  dec.fet.test = fisher.test(contingency.matrix.dec,alternative='greater')
  inc.kst.test = ks.test(x$effect,subset(all.res.do.pass,do_id != i)$effect,alternative='less')
  dec.kst.test = ks.test(x$effect,subset(all.res.do.pass,do_id != i)$effect,alternative='greater')
  
  data.frame(
    do_id = unique(x$do_id),
    do.size = sum(this.deg.total),
    deg.n = this.deg.total[1],
    inc.n = this.inc.total[1],
    dec.n = this.dec.total[1],
    deg.fet.score = deg.fet.test$estimate,
    inc.fet.score = inc.fet.test$estimate,
    dec.fet.score = dec.fet.test$estimate,
    inc.kst.score = inc.kst.test$statistic,
    dec.kst.score = dec.kst.test$statistic,
    deg.fet.pval = deg.fet.test$p.value,
    inc.fet.pval = inc.fet.test$p.value,
    dec.fet.pval = dec.fet.test$p.value,
    inc.kst.pval = inc.kst.test$p.value,
    dec.kst.pval = dec.kst.test$p.value
  )
},mc.cores=n.cores))

all.res.do.test = within(all.res.do.test,{
  dec.kst.qval = p.adjust(dec.kst.pval,'BH')
  inc.kst.qval = p.adjust(inc.kst.pval,'BH')
  dec.fet.qval = p.adjust(dec.fet.pval,'BH')
  inc.fet.qval = p.adjust(inc.fet.pval,'BH')
  deg.fet.qval = p.adjust(deg.fet.pval,'BH')
})

all.res.do.results = merge(all.res.do.test,do.def,by='do_id')
all.res.do.results$dataset = 'DISEASES'
View(all.res.do.results)

write.csv(all.res.do.results, file = paste(contrastnow,"DO.csv",sep=""))

# plot

pl = read.csv(file = "HvMDO.csv")
pl = subset(pl, do_name %in% c('Autistic disorder',
                                'Bipolar disorder',
                                'Schizophrenia',
                                'Intellectual disability',
                                'Attention deficit hyperactivity disorder'))
pl$do_name = as.factor(pl$do_name)
levels(pl$do_name)
levels(pl$do_name) = c('ADHD','ASD','BPD','ID','SCZ')
pl$do_name <- factor(pl$do_name, levels=c('ID','ASD','BPD','SCZ','ADHD'))

ggplot(pl, aes(x=do_name, y = -log10(dec.kst.qval))) +
  geom_bar(stat = 'identity', fill = '#7570b3') +
  coord_flip() +
  ylab("Adjusted P-value (-log10)") +
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed', color = 'grey') +
  theme_classic() +
  theme(axis.text = element_text(size=18),
        legend.position = "none", 
        legend.title = element_blank(),
        axis.title.x = element_text(size=18),
        legend.text = element_text(size=18),
        axis.title.y = element_blank()) 

############
## estimate cell type proportions (Figure S5)
###########

conv_121 = readRDS('conv_121.rds')
markers = markers_df_mouse_brain
markers$markers = str_to_upper(markers$markers)
colnames(markers)[1] = 'external_gene_name'
markers = merge(markers, conv_121, by = 'external_gene_name')
markers = markers[,c(4,2)]
colnames(markers) = c('markers','cell')

#spv = findCells(inputMat=y$E, markers=markers, nMarker = 50, method = "SVD", scale = TRUE)
spv = findCells(inputMat=nobatch, markers=markers, nMarker = 50, method = "SVD", scale = TRUE)
spv = reshape2::melt(spv)
colnames(spv)[1] = 'ID'
spv = merge(spv, metanow, by = 'ID')
levels(spv$Var2) = c('Astrocytes','Endothelial','Oligodendrocytes','Neurons','Microglia','OPCs')
spv$Species = factor(spv$Species, levels = c('Human','Squirrel.Monkey','Macaque'))
levels(spv$Species) = c('HIMs','SMIMs','MIMs')

ggplot(spv, aes(x=Species, y = value, fill = Species)) +
  geom_boxplot() +
  scale_fill_manual(values = species.cols[c(1,3,2)]) +
  facet_wrap(~Var2) +
  ylab('Relative cell type proportions (SPVs)') +
  geom_hline(yintercept = 0, color = 'grey', linetype = 'dashed') +
  theme_article() +
  theme(legend.position = 'none',
        axis.text.y = element_text(size = 16),
        strip.text.x = element_text(size = 16),
        axis.text.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_blank())

anova_results <- data.frame()
tukey_results <- data.frame()
mean_values <- data.frame()

for(i in seq_along(levels(spv$Var2))){
  
  var_name <- levels(spv$Var2)[i]
  now <- subset(spv, Var2 == var_name)
  
  group_means <- now %>%
    group_by(Species) %>%
    summarise(Mean = mean(value), .groups = 'drop') %>%
    mutate(Variable = var_name)
  
  mean_values <- rbind(mean_values, group_means)
  
  # Run ANOVA
  mod <- aov(value ~ Species, data = now)
  anova_summary <- summary(mod)[[1]]
  anova_row <- data.frame(
    Variable = var_name,
    Df = anova_summary["Species", "Df"],
    F_value = anova_summary["Species", "F value"],
    P_value = anova_summary["Species", "Pr(>F)"]
  )
  anova_results <- rbind(anova_results, anova_row)
  
  # Run Tukey HSD
  tu <- TukeyHSD(mod)
  tu_df <- as.data.frame(tu$Species)
  tu_df$Comparison <- rownames(tu_df)
  tu_df$Variable <- var_name
  tukey_results <- rbind(tukey_results, tu_df)
}

tukey_results <- tukey_results[, c("Variable", "Comparison", "diff", "lwr", "upr", "p adj")]

colnames(tukey_results) <- c("Variable", "Comparison", "Difference", "Lower_CI", "Upper_CI", "Adj_P")
colnames(anova_results) <- c("Variable", "Df", "F_value", "P_value")
anova_results$padj = p.adjust(anova_results$P_value, method = 'BH')
colnames(mean_values) <- c("Species", "Mean", "Variable")

write.csv(tukey_results, file = 'tukey_results.csv')
write.csv(anova_results, file = 'anova_results.csv')
write.csv(mean_values, file = 'mean_values.csv')

#######################
## plot expression of specific genes (Figure 4B)
#######################

## DLG4 ENSMUSG00000020886	ENSG00000132535
## SYP ENSMUSG00000031144	ENSG00000102003

genenow = 'ENSMUSG00000020886'
genenow = 'ENSMUSG00000031144'

expnow = nobatch[genenow,]
expnow = data.frame(expnow)
expnow$ID = rownames(expnow)
expnow = merge(expnow, metanow[,c('ID','Species')], by = 'ID')
expnow$Species = factor(expnow$Species, levels = c('Human','Squirrel.Monkey','Macaque'))
levels(expnow$Species) = c('HIMs','SMIMs','MIMs')

ggplot(expnow, aes(x = Species, y = expnow, fill = Species)) +
  geom_boxplot() +
  scale_fill_manual(values = species.cols[c(1,3,2)]) +
  theme_article() +
  ylab('Normalized Expression') +
  theme(legend.position = 'none',
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 18),
        axis.text = element_text(size = 18))
