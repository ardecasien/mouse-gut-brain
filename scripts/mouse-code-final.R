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

`%!in%` = Negate(`%in%`)

#################
## load data
#################

fc = read.csv("counts.csv")
rownames(fc) = fc$X
fc = fc[,-1]

meta = read.csv("mouse_meta.csv")
meta.tech = read.csv("mouse_tech_meta.csv")
meta = merge(meta, meta.tech, by = 'ID', all = T)
meta$reads = as.numeric(meta$reads)
meta$scaled_mapping = scale(meta$mapping_rate)
meta$scaled_reads = scale(meta$reads)

# remove outliers (identified using plotMDS below)

rm.fc = c("X908", "X874", "X923")
fc = fc[,colnames(fc) %!in% rm.fc]

# Zhu et al. data

hum_mac1a = readRDS('human_v_macaque_zhu_MFC.rds')

###################
## filter genes
###################

samples.by.species.and.age = split(meta$ID, list(meta$Species,meta$Age), drop = TRUE)

cpm.cutoff = 10

keeps.fc = data.frame(gene = NULL)
for (i in 1:length(samples.by.species.and.age)){
  keep.now = data.frame(gene = names(which(rowMeans(cpm(fc[,colnames(fc)%in%samples.by.species.and.age[[i]]])) >= cpm.cutoff)))
  keeps.fc = rbind(keeps.fc, keep.now)}
keeps.fc = unique(keeps.fc)

####################
## voom 
####################

datanow = fc[rownames(fc) %in% keeps.fc$gene,]
metanow = subset(meta, ID %in% colnames(fc))
metanow = metanow[order(match(metanow$ID, colnames(datanow))),]

design <- model.matrix(~0 + scaled_mapping + scaled_reads + Species.Batch, data = metanow)
metanow$Species.Batch = as.factor(metanow$Species.Batch)
levels(metanow$Species.Batch)
contrasts <- makeContrasts(
  HvM=(Species.BatchHuman.A+Species.BatchHuman.B)/2-(Species.BatchMacaque.A+Species.BatchMacaque.B)/2,
  HvS=(Species.BatchHuman.A+Species.BatchHuman.B)/2-(Species.BatchSquirrel.Monkey.A+Species.BatchSquirrel.Monkey.B)/2,
  SvM=(Species.BatchSquirrel.Monkey.A+Species.BatchSquirrel.Monkey.B)/2-(Species.BatchMacaque.A+Species.BatchMacaque.B)/2,
  HvSM=(Species.BatchHuman.A+Species.BatchHuman.B)/2-(Species.BatchSquirrel.Monkey.A+Species.BatchSquirrel.Monkey.B+Species.BatchMacaque.A+Species.BatchMacaque.B)/4,
  HSvM=(Species.BatchHuman.A+Species.BatchHuman.B+Species.BatchSquirrel.Monkey.A+Species.BatchSquirrel.Monkey.B)/4-(Species.BatchMacaque.A+Species.BatchMacaque.B)/2,
  AvY=(Species.BatchHuman.A+Species.BatchSquirrel.Monkey.A+Species.BatchMacaque.A)/3-(Species.BatchHuman.B+Species.BatchMacaque.B+Species.BatchSquirrel.Monkey.B)/3,
  HAvHY=Species.BatchHuman.A-Species.BatchHuman.B,
  SAvSY=Species.BatchSquirrel.Monkey.A-Species.BatchSquirrel.Monkey.B,
  MAvMY=Species.BatchMacaque.A-Species.BatchMacaque.B,
  HAvALL=Species.BatchHuman.A-(Species.BatchHuman.B+Species.BatchSquirrel.Monkey.A+Species.BatchSquirrel.Monkey.B+Species.BatchMacaque.A+Species.BatchMacaque.B)/5,
  HYvALL=Species.BatchHuman.B-(Species.BatchHuman.A+Species.BatchSquirrel.Monkey.A+Species.BatchSquirrel.Monkey.B+Species.BatchMacaque.A+Species.BatchMacaque.B)/5,
  SAvALL=Species.BatchSquirrel.Monkey.A-(Species.BatchHuman.B+Species.BatchHuman.A+Species.BatchSquirrel.Monkey.B+Species.BatchMacaque.A+Species.BatchMacaque.B)/5,
  SYvALL=Species.BatchSquirrel.Monkey.B-(Species.BatchHuman.B+Species.BatchSquirrel.Monkey.A+Species.BatchHuman.A+Species.BatchMacaque.A+Species.BatchMacaque.B)/5,
  MAvALL=Species.BatchMacaque.A-(Species.BatchHuman.B+Species.BatchSquirrel.Monkey.A+Species.BatchSquirrel.Monkey.B+Species.BatchHuman.A+Species.BatchMacaque.B)/5,
  MYvALL=Species.BatchMacaque.B-(Species.BatchHuman.B+Species.BatchSquirrel.Monkey.A+Species.BatchSquirrel.Monkey.B+Species.BatchMacaque.A+Species.BatchHuman.A)/5,
  levels=colnames(design))

rownames(contrasts) = gsub("Species.Batch", "", rownames(contrasts))
contrasts
round(rowSums(t(contrasts)),3)

dgList <- DGEList(counts=datanow, genes=rownames(datanow), 
                  group  = factor(metanow$Species.Batch)
)
dgList <- calcNormFactors(dgList)

y = voomWithQualityWeights(dgList, design, plot = F)
saveRDS(y, 'voom_y.rds')

#####################
# variance partitioning (Figure 2A)
#####################

y = readRDS('voom_y.rds')

n.cores = detectCores() - 4
param = SnowParam(n.cores, "SOCK", progressbar=TRUE)
design = as.formula(paste('~',paste(c('scaled_mapping','scaled_reads','(Species)','(Batch)'),collapse=' + ')))
varPart = fitExtractVarPartModel(y$E, design, metanow, BPPARAM=param)
vp = sortCols(varPart)
rownames(vp) = vp$X
vp = vp[,-1]
cols = c('#1b9e77', '#7570b3', '#e6ab02','#e7298a')
colnames(vp) = c('Gut species','Mapping rate','Read count','Batch','Residuals')

mean(vp$`Gut species`) * 100
mean(vp$`Mapping rate`) * 100
mean(vp$`Read count`) * 100
mean(vp$Batch) * 100
mean(vp$Residuals) * 100

plotVarPart(vp, col = c(cols,'grey')) +
  geom_boxplot(width = 0.07, alpha=0.8, fill="grey") +
  theme_article() +
  ylab('% variance explained') +
  theme(axis.title = element_text(size=18), 
        axis.text = element_text(size=18),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") 

######################
# plot expression data (Figure 3A)
######################

# remove batch and covariate effects

nobatch = removeBatchEffect(y$E, batch = metanow$Batch, covariates = cbind(metanow$scaled_mapping, metanow$scaled_reads))

# MDS plot

plotMDS(nobatch, col = as.numeric(metanow$Species))

# PCA

pca = prcomp(cor(nobatch))
plot.pca = as.data.frame(pca$x)
plot.pca$Species = metanow$Species
plot.pca$Species[which(plot.pca$Species == 'Squirrel.Monkey')] = 'Squirrel Monkey'
plot.pca$Batch = factor(metanow$Batch)

ggplot(plot.pca,aes(x=PC1,y=PC2,color=Species, shape=Batch)) + 
  scale_color_manual(values=c(cols)) +
  geom_point(size=5) + 
  theme_classic() + 
  xlab('PC 1') + ylab('PC 2') + 
  theme(axis.ticks=element_blank(),
        axis.text=element_blank(),
        axis.title = element_text(size=18), 
        legend.text = element_text(size=14), 
        legend.title = element_blank())

# tSNE

a = Rtsne(t(nobatch), dims = 2, perplexity=10, verbose=TRUE, max_iter = 1000)
plot.tsne = as.data.frame(a$Y)
plot.tsne$Species = metanow$Species
plot.tsne$Species[which(plot.tsne$Species == 'Squirrel.Monkey')] = 'Squirrel Monkey'
plot.tsne$Batch = factor(metanow$Batch)

ggplot(plot.tsne,aes(x=V1,y=V2,color=Species, shape=Batch)) + 
  scale_color_manual(values=c(cols)) +
  geom_point(size=5) + 
  theme_classic() + 
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

ggplot(a.umap,aes(x=V1,y=V2,color=Species,shape =Batch)) + 
  scale_color_manual(values=c(cols)) +
  geom_point(size=5) + 
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

colnames(design) = gsub("Species.Batch", "", colnames(design))
head(design)

fit <- lmFit(y, design)
fit <- contrasts.fit(fit, contrasts)
fit <- eBayes(fit)
coef = fit$coefficients
pvals = fit$p.value
stdev = fit$stdev.unscaled

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
saveRDS(combo, 'combo_results.rds')

combo = readRDS('combo_results.rds')
colnames(combo)[1] = 'mmusculus_homolog_ensembl_gene'
mert = useMart(biomart = 'ENSEMBL_MART_ENSEMBL', dataset='hsapiens_gene_ensembl')
conv = getBM(attributes=c('ensembl_gene_id','external_gene_name','mmusculus_homolog_ensembl_gene','mmusculus_homolog_orthology_type'), mart = mert)
conv = getBM(attributes=c('ensembl_gene_id','external_gene_name','gene_biotype'), mart = mert)
conv_121 = subset(conv, mmusculus_homolog_orthology_type == 'ortholog_one2one')
combo2 = merge(combo, conv_121, by = 'mmusculus_homolog_ensembl_gene', all.x = T)

write.csv(combo2, file = 'model_results_qualweights_tech.csv')

################################
# volcano plots (Figure 2B-D)
################################

combo = read.csv('model_results_qualweights_tech.csv')

pl = subset(combo, contrast == 'HvM')
pl$code = ifelse(pl$padj < 0.1 & pl$beta > 0, 'human', ifelse(pl$padj < 0.1 & pl$beta < 0, 'macaque', 'none'))
colsnow = c(cols[c(1,2)],'#666666')
pl$ribo = ifelse(str_sub(pl$external_gene_name,1,3) %in% c('RPL','RPS'),'ribo','gene')
pl$ribo = ifelse(str_sub(pl$external_gene_name,1,4) %in% c('MRPL','MRPS'),'ribo',pl$ribo)
table(pl$ribo, pl$code)

pl = subset(combo, contrast == 'HvS')
pl$code = ifelse(pl$padj < 0.1 & pl$beta > 0, 'human', ifelse(pl$padj < 0.1 & pl$beta < 0, 'sm', 'none'))
colsnow = c(cols[1],'#666666',cols[3])
pl$ribo = ifelse(str_sub(pl$external_gene_name,1,3) %in% c('RPL','RPS'),'ribo','gene')
pl$ribo = ifelse(str_sub(pl$external_gene_name,1,4) %in% c('MRPL','MRPS'),'ribo',pl$ribo)
table(pl$ribo, pl$code)

pl = subset(combo, contrast == 'SvM')
pl$code = ifelse(pl$padj < 0.1 & pl$beta > 0, 'sq', ifelse(pl$padj < 0.1 & pl$beta < 0, 'macaque', 'none'))
colsnow = c(cols[2],'#666666',cols[3])
pl$ribo = ifelse(str_sub(pl$external_gene_name,1,3) %in% c('RPL','RPS'),'ribo','gene')
pl$ribo = ifelse(str_sub(pl$external_gene_name,1,4) %in% c('MRPL','MRPS'),'ribo',pl$ribo)
table(pl$ribo, pl$code)

table(pl$code)

ggplot(pl, aes(x = beta, y = -log10(padj), color = code)) +
  geom_point(alpha = 0.5) + 
  scale_color_manual(values = colsnow) +
  xlim(c(-2,2)) +
  ylim(c(0,8)) +
  ylab('-log10(padj)') +
  xlab(expression(~italic(beta))) +
  geom_hline(yintercept = -log10(0.1), linetype = 'dashed') +
  geom_text_repel(data=subset(pl, padj < 0.1),
                  aes(x = beta, y = -log10(padj),label=external_gene_name), 
                  color = 'black', size = 3, max.overlaps = 10) +
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
## venn diagrams (Figure 3B)
################

myCol <- brewer.pal(3, "Dark2")

mouse_res = read.csv('model_results_qualweights_tech.csv')
hvm = subset(mouse_res, contrast == 'HvM' & padj < 0.1)$gene
hvs = subset(mouse_res, contrast == 'HvS' & padj < 0.1)$gene
svm = subset(mouse_res, contrast == 'SvM' & padj < 0.1)$gene

venn.diagram(x = list(hvm, hvs, svm),
             category.names = c("HIMs v MIMs" , "HIMs v SMIMs" , "SMIMs v MIMs"),
             file = 'species_venn.png',
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

coef.plot = melt(coef)

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

# run using outputs from microbiome.sh

mouse_res = read.csv('model_results_qualweights_tech.csv')
mouse_res = subset(mouse_res, contrast == 'HvM')
dim(mouse_res)

de_hum = subset(mouse_res, padj < 0.1 & beta > 0)
de_mac = subset(mouse_res, padj < 0.1 & beta < 0)

gene_clusters = read.csv('metabolite_clusters.csv')
dim(gene_clusters)
table(gene_clusters$X %in% mouse_res$mmusculus_homolog_ensembl_gene)
cl = unique(gene_clusters$Cluster)

o = data.frame()
for(i in 1:length(cl)) {
  
  clnow = subset(gene_clusters, Cluster == cl[i])
  dim(clnow)
  
  b = length(subset(clnow, X %in% de_hum$mmusculus_homolog_ensembl_gene)$X)
  d = length(subset(de_hum, mmusculus_homolog_ensembl_gene %!in% clnow$X)$mmusculus_homolog_ensembl_gene)
  c = length(subset(clnow, X %!in% de_hum$mmusculus_homolog_ensembl_gene)$X)
  n = length(subset(mouse_res, mmusculus_homolog_ensembl_gene %!in% c(clnow$X, de_hum$mmusculus_homolog_ensembl_gene))$mmusculus_homolog_ensembl_gene)
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
# OR > 1 and padj < 0.05 for clusters 1,3,0

o$cluster = as.character(o$cluster)
ggplot(o, aes(x = cluster, y = OR)) +
  geom_bar(stat = 'identity', fill = '#E6AB02') +
  xlab('Gene cluster') +
  ylab('Enrichments for HIMs > MIMs DE genes\n(odds ratios)') +
  theme_article() 

################################
## compare to primate data (Figures 2E 2F)
################################

mouse_res = read.csv('model_results_qualweights_tech.csv')
mouse_res = subset(mouse_res, contrast == 'HvM')

prim_res = hum_mac1a
prim_res$ensembl_gene_id = rownames(prim_res)

combo = merge(prim_res, mouse_res, by = 'ensembl_gene_id')
write.csv(combo, file = 'mouse_human_combo_tech.csv')
combo$code = ifelse(combo$beta.x > 0 & combo$beta.y > 0, 'concordant', 'discordant')
combo$code = ifelse(combo$beta.x < 0 & combo$beta.y < 0, 'concordant', combo$code)

combo3 = subset(combo, padj.x < 0.2 & padj.y < 0.2)
dim(combo3)
combo3$code = ifelse(combo3$beta.x > 0 & combo3$beta.y > 0, 'Concordant', 'Discordant')
combo3$code = ifelse(combo3$beta.x < 0 & combo3$beta.y < 0, 'Concordant', combo3$code)
table(combo3$code)
table(combo3$code)[1] / sum(table(combo3$code))

cor.test(x = combo3$beta.x, y = combo3$beta.y, method = 'spearman')
upup = length(subset(combo3, beta.x > 0 & combo3$beta.y > 0)$external_gene_name)
downdown = length(subset(combo3, beta.x < 0 & combo3$beta.y < 0)$external_gene_name)
updown = length(subset(combo3, beta.x > 0 & combo3$beta.y < 0)$external_gene_name)
downup = length(subset(combo3, beta.x < 0 & combo3$beta.y > 0)$external_gene_name)
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

mouse_res = read.csv('model_results_qualweights_tech.csv')
dumas = read.csv('dumas.csv')
colnames(dumas)[2] = 'external_gene_name'
dumas$mean = rowMeans(dumas[,c(3:5)], na.rm = T)

combo = merge(mouse_res, dumas, by = 'external_gene_name')

contrastnow = 'HvM'
datanow = subset(combo, contrast == contrastnow)
datanow = datanow[complete.cases(datanow$ensembl_gene_id),]

datanow = subset(datanow, abs(HS) > 0.5)
datanow = subset(datanow, padj < 0.1)

table(datanow$beta < 0, datanow$HS < 0)
m = matrix(table(datanow$beta < 0, datanow$HS < 0), nrow = 2, ncol = 2)
rownames(m) = c('HIMs','MIMs')
colnames(m) = c('positive','purify')
m
sum(m)
fisher.test(m, alternative = 'greater')

############################
## biological processes (Figure 4A)
############################

mouse_res = read.csv('model_results_qualweights_tech.csv')
prim_res = hum_mac1a
prim_res$ensembl_gene_id = rownames(prim_res)
combo = merge(prim_res, mouse_res, by = 'ensembl_gene_id', all = T)

Ensembl<-listEnsembl()
match.arg("genes",Ensembl$biomart)
mart<-useEnsembl("genes")
Ensembl = new("genomic_ressource",
              db="Ensembl",
              stamp=paste("www.ensembl.org",Ensembl$version[Ensembl$biomart=="genes"]),
              data=data.table(),
              mart=list(mart),
              organisms=data.table(listDatasets(mart)))
myGENE2GO<-ViSEAGO::annotate(
  "hsapiens_gene_ensembl",
  Ensembl)
saveRDS(myGENE2GO, file = 'myGENE2GO.rds')

myGENE2GO = readRDS('myGENE2GO.rds')

# select contrast

contrastnow = 'HvM'
#contrastnow = 'HvS'
#contrastnow = 'SvM'

datanow = subset(combo, contrast == contrastnow)
datanow = datanow[complete.cases(datanow$ensembl_gene_id),]

allgenes = unique(datanow$ensembl_gene_id)
all.region = numeric(length=length(allgenes))
names(all.region) = allgenes

up = unique(subset(datanow, contrast == contrastnow & padj.y < 0.1 & beta.y > 0)$ensembl_gene_id)
length(up)
down = unique(subset(datanow, contrast == contrastnow & padj.y < 0.1 & beta.y < 0)$ensembl_gene_id)
length(down)

all.region[which(names(all.region) %in% up)] = 1
all.region[which(names(all.region) %in% down)] = -1
table(all.region)
all.region

# run models

BP.m<-ViSEAGO::create_topGOdata(
  geneSel=names(all.region[which(all.region == 1)]),
  allGenes=names(all.region),
  gene2GO=myGENE2GO, 
  ont="BP",
  nodeSize=5)

BP.f<-ViSEAGO::create_topGOdata(
  geneSel=names(all.region[which(all.region == -1)]),
  allGenes=names(all.region),
  gene2GO=myGENE2GO, 
  ont="BP",
  nodeSize=5)

test.m<-topGO::runTest(
  BP.m,
  algorithm ="parentchild",
  statistic = "fisher")

test.f<-topGO::runTest(
  BP.f,
  algorithm ="parentchild",
  statistic = "fisher")

BP_sResults<-merge_enrich_terms(
  Input=list(up=c("BP.m","test.m"),down=c("BP.f","test.f")),
  cutoff = 0.05)

BP_sResults2 = BP_sResults@data

myGOs<-ViSEAGO::build_GO_SS(
  gene2GO=myGENE2GO,
  enrich_GO_terms=BP_sResults)

myGOs<-ViSEAGO::compute_SS_distances(
  myGOs,
  distance="Wang")

Wang_clusters_wardD2<-ViSEAGO::GOterms_heatmap(
  myGOs,
  showIC=TRUE,
  showGOlabels=TRUE,
  GO.tree=list(
    tree=list(
      distance="Wang",
      aggreg.method="ward.D2"
    ),
    cut=list(
      dynamic=list(
        pamStage=TRUE,
        pamRespectsDendro=TRUE,
        deepSplit=2,
        minClusterSize =2))),
  samples.tree=NULL)

ViSEAGO::show_heatmap(
  Wang_clusters_wardD2,
  "GOterms")

Wang = Wang_clusters_wardD2@enrich_GOs@data

saveRDS(Wang_clusters_wardD2, file = paste(contrastnow, '_Wang_clusters_wardD2_tech.rds', sep = ""))
write.csv(Wang, file = paste(contrastnow, '_GO_tech.csv', sep = ""))

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

do.mmul = subset(do.all,
                 select=c('external_gene_name','ensembl_gene_id','ensembl_peptide_id','do_id','do_name','z_score','confidence'))

# select contrast

contrastnow = 'HvM'
#contrastnow = 'HvS'
#contrastnow = 'SvM'

datanow = subset(combo, contrast == contrastnow)
datanow = datanow[complete.cases(datanow$ensembl_gene_id),]

allgenes = unique(datanow$ensembl_gene_id)
all.region.fet = all.region.kst = numeric(length=length(allgenes))
names(all.region.fet) = names(all.region.kst) = allgenes

up = unique(subset(datanow, padj.y < 0.1 & beta.y > 0)$ensembl_gene_id)
length(up)
down = unique(subset(datanow, padj.y < 0.1 & beta.y < 0)$ensembl_gene_id)
length(down)

all.region.fet[which(names(all.region.fet) %in% up)] = 1
all.region.fet[which(names(all.region.fet) %in% down)] = -1
table(all.region.fet)
names(all.region.fet) = allgenes

all.region.kst = datanow$beta.y / datanow$stdev.y

all.region.join = data.frame(ensembl_gene_id = allgenes, direction = as.integer(all.region.fet), effect = as.numeric(all.region.kst))
all.region.do = merge(all.region.join, do.mmul, by='ensembl_gene_id')

all.region.do.pass = subset(all.region.do,do_id %in% names(which(table(subset(all.region.do,confidence >= 0)$do_id) >= 10)))
all.region.do.gene.pass = unique(all.region.do.pass[c('ensembl_gene_id','direction')])

do.deg.total = as.integer(table(factor(all.region.do.pass$direction != 0,levels=c('TRUE','FALSE'))))
do.inc.total = as.integer(table(factor(all.region.do.pass$direction == 1 | all.region.do.pass$direction == 2,levels=c('TRUE','FALSE'))))
do.dec.total = as.integer(table(factor(all.region.do.pass$direction == -1 | all.region.do.pass$direction == 2,levels=c('TRUE','FALSE'))))

all.region.do.split = split(all.region.do.pass,all.region.do.pass$do_id)

n.cores = detectCores() - 6

all.region.do.test = do.call(rbind,mclapply(names(all.region.do.split),function(i) {
  x = all.region.do.split[[i]]
  
  this.deg.total = as.integer(table(factor(x$direction != 0,levels=c('TRUE','FALSE'))))
  contingency.matrix.deg = matrix(rbind(this.deg.total,do.deg.total - this.deg.total),nrow=2,dimnames=list(c('in group','not in group'),c('in direction','not in direction')))
  
  this.inc.total = as.integer(table(factor(x$direction == 1 | x$direction == 2,levels=c('TRUE','FALSE'))))
  contingency.matrix.inc = matrix(rbind(this.inc.total,do.inc.total - this.inc.total),nrow=2,dimnames=list(c('in group','not in group'),c('in direction','not in direction')))
  
  this.dec.total = as.integer(table(factor(x$direction == -1 | x$direction == 2,levels=c('TRUE','FALSE'))))
  contingency.matrix.dec = matrix(rbind(this.dec.total,do.dec.total - this.dec.total),nrow=2,dimnames=list(c('in group','not in group'),c('in direction','not in direction')))
  
  deg.fet.test = fisher.test(contingency.matrix.deg,alternative='greater')
  inc.fet.test = fisher.test(contingency.matrix.inc,alternative='greater')
  dec.fet.test = fisher.test(contingency.matrix.dec,alternative='greater')
  inc.kst.test = ks.test(x$effect,subset(all.region.do.pass,do_id != i)$effect,alternative='less')
  dec.kst.test = ks.test(x$effect,subset(all.region.do.pass,do_id != i)$effect,alternative='greater')
  
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

all.region.do.test = within(all.region.do.test,{
  dec.kst.qval = p.adjust(dec.kst.pval,'BH')
  inc.kst.qval = p.adjust(inc.kst.pval,'BH')
  dec.fet.qval = p.adjust(dec.fet.pval,'BH')
  inc.fet.qval = p.adjust(inc.fet.pval,'BH')
  deg.fet.qval = p.adjust(deg.fet.pval,'BH')
})

all.region.do.results = merge(all.region.do.test,do.def,by='do_id')
all.region.do.results$dataset = 'DISEASES'
all.region.do.results$set = 'union'
all.region.do.results$region = 'all'

write.csv(all.region.do.results, file = paste(contrastnow,"DO.csv",sep=""))

############
## estimate cell type proportions (Figure 4D)
###########

markers = markers_df_mouse_brain
markers$markers = str_to_upper(markers$markers)
colnames(markers)[1] = 'external_gene_name'
markers = merge(markers, conv, by = 'external_gene_name')
markers = markers[,c(4,2)]
colnames(markers) = c('markers','cell')

spv = findCells(inputMat=y$E, markers=markers, nMarker = 50, method = "SVD", scale = TRUE)
spv = reshape2::melt(spv)
colnames(spv)[1] = 'ID'
spv = merge(spv, metanow, by = 'ID')
levels(spv$Var2) = c('Astrocytes','Endothelial','Oligodendrocytes','Neurons','Microglia','OPCs')
spv$Species = factor(spv$Species, levels = c('Human','Squirrel.Monkey','Macaque'))
levels(spv$Species) = c('HIMs','SMIMs','MIMs')

ggplot(spv, aes(x=Species, y = value, fill = Species)) +
  geom_boxplot() +
  scale_fill_manual(values = cols) +
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

for(i in 1:length(levels(spv$Var2))){
  print(levels(spv$Var2)[i])
  now = subset(spv, Var2 == levels(spv$Var2)[i])
  bb = now %>% group_by(Species) %>% summarise(mean = mean(value))
  print(bb)
  mod = aov(value ~ Species, data = now)
  print(summary(mod))
  tu = TukeyHSD(mod)
  print(tu)
}

#######################
## plot residual expression of specific genes (Figure 4B)
#######################

## DLG4 ENSMUSG00000020886	ENSG00000132535
## MOG ENSMUSG00000076439	ENSG00000204655
## MOBP ENSMUSG00000032517	ENSG00000168314
## SLC2A1 ENSMUSG00000028645	ENSG00000117394
## APOE ENSMUSG00000002985	ENSG00000130203
## SYP ENSMUSG00000031144	ENSG00000102003

genenow = 'ENSMUSG00000020886'
genenow = 'ENSMUSG00000076439'
genenow = 'ENSMUSG00000032517'
genenow = 'ENSMUSG00000028645'
genenow = 'ENSMUSG00000002985'
genenow = 'ENSMUSG00000031144'

expnow = nobatch[genenow,]
expnow = data.frame(expnow)
expnow$ID = rownames(expnow)
expnow = merge(expnow, metanow[,c('ID','Species')], by = 'ID')
expnow$Species = factor(expnow$Species, levels = c('Human','Squirrel.Monkey','Macaque'))
levels(expnow$Species) = c('HIMs','SMIMs','MIMs')

ggplot(expnow, aes(x = Species, y = expnow, fill = Species)) +
  geom_boxplot() +
  scale_fill_manual(values = cols[c(1,3,2)]) +
  theme_article() +
  ylab('Normalized Expression') +
  theme(legend.position = 'none',
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 18),
        axis.text = element_text(size = 18))

