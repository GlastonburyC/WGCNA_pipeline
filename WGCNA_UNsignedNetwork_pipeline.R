library(WGCNA)
library(data.table)
allowWGCNAThreads()

rpkm<-fread("~/../DTR/Expression/EUROBATS/Counts/Residuals_RPKM_F_All_family.txt",head=F,sep="\t")
rpkm<-as.data.frame(rpkm)
rpkm$V772<-NULL
info=rpkm[,1:5]
rpkm=rpkm[,-1:-5]

datExpr0=rpkm

datExpr0=t(datExpr0)

datExpr0[1:5,1:5]
covs<-read.table('~/../DTR/Expression/EUROBATS/Counts/qc_F_freezev2.txt',head=T,sep="\t")

row.names(datExpr0)= covs$EurobatsID
covs<-read.table('~/../DTR/Expression/EUROBATS/Counts/qc_F_freezev1.txt',head=T,sep="\t")

datExpr0=datExpr0[row.names(datExpr0) %in% covs$SampleID,]


sampleTree = hclust(dist(datExpr0), method = "average");
pdf(file = "sampleClustering.pdf", width = 20, height = 15);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
    cex.axis = 1.5, cex.main = 2)
dev.off()

abline(h = 200, col = "red");
clust = cutreeStatic(sampleTree, cutHeight = 200, minSize = 10)
table(clust)
keepSamples = (clust==1)
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)


# Samples removed, recluster.
sampleTree = hclust(dist(datExpr), method = "average");
pdf(file = "sampleClustering_outliersRM.pdf", width = 20, height = 15);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering without outliers", sub="", xlab="", cex.lab = 1.5,
    cex.axis = 1.5, cex.main = 2)
dev.off()


options(stringsAsFactors = FALSE);


powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
pdf(file = "pickSoftThresholdingPower.pdf", width = 9, height = 9);
par(mfrow = c(1,2));
cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
    main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
    labels=powers,cex=cex1,col="red");

plot(sft$fitIndices[,1], sft$fitIndices[,5],
    xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
    main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

dev.off()

# SET SOFTPOWER THRESHOLD BASED ON RESULTS OF SCALE FREE TOPOLOGY PLOTS
softPower = 10 ; # Using a unsigned network
adjacency = adjacency(datExpr, power = softPower, type="unsigned");


# Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency, TOMType="unsigned");
dissTOM = 1-TOM

# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)
pdf(file = "GeneClustering.pdf", width = 12, height = 9);
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
    labels = FALSE, hang = 0.04);

dev.off()
# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 30;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
              deepSplit = 2, pamRespectsDendro = FALSE,
              minClusterSize = minModuleSize);
table(dynamicMods)

# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
pdf(file = "GeneClustering_Cut_Modules.pdf", width = 12, height = 9);
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                  dendroLabels = FALSE, hang = 0.03,
                  addGuide = TRUE, guideHang = 0.05,
                  main = "Gene dendrogram and module colors")
dev.off()

# MERGE SIMILAR MODULES BASED ON EIGENGENE CORRELATION > 0.75.
# Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
pdf(file = "Clustering_of_ME.pdf", width = 12, height = 9);
plot(METree, main = "Clustering of module eigengenes",
xlab = "", sub = "")

dev.off()

MEDissThres = 0.25
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
mergedColors = merge$colors;
mergedMEs = merge$newMEs;


pdf(file = "gene_dendrogram_merged_modules.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
c("Dynamic Tree Cut", "Merged dynamic"),
dendroLabels = FALSE, hang = 0.03,
addGuide = TRUE, guideHang = 0.05)
dev.off()

moduleColors = mergedColors


# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;
# Save module colors and labels for use in subsequent parts
save(MEs, moduleLabels, moduleColors, geneTree, file="WGCNA_module_construction.RData")

covs=covs[covs[,1] %in% row.names(datExpr),]
 datTraits=covs$BMI

# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

pdf('MEs_cor_BMI.pdf',height=9,width=12)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
xLabels = names(datTraits),
yLabels = names(MEs),
ySymbols = names(MEs),
colorLabels = FALSE,
colors = greenWhiteRed(50),
textMatrix = textMatrix,
setStdMargins = FALSE,
cex.text = 0.6,
zlim = c(-1,1),
main = paste("Module-trait relationships"))
dev.off()

colnames(datExpr)=info$V2
# Define variable weight containing the weight column of datTrait
BMI = as.data.frame(datTraits);
names(BMI) = "BMI"
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datExpr, BMI, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(BMI), sep="");
names(GSPvalue) = paste("p.GS.", names(BMI), sep="");

module = "lightcyan"
column = match(module, modNames);
moduleGenes = moduleColors==module;
pdf('BMI_posmodule_MM_GS.pdf',height=7,width=7)
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
abs(geneTraitSignificance[moduleGenes, 1]),
xlab = paste("Module Membership in", module, "module"),
ylab = "Gene significance for BMI",
main = paste("Module membership vs. gene significance\n"),
cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
dev.off()

colnames(datExpr)[moduleColors=="lightcyan"]

geneInfo0 = data.frame(geneSymbol = colnames(datExpr),
moduleColor = moduleColors,
geneTraitSignificance,
GSPvalue)
# Order modules by their significance for BMI
modOrder = order(-abs(cor(MEs, BMI, use = "p")));
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
oldNames = names(geneInfo0)
geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
MMPvalue[, modOrder[mod]]);
names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
paste("p.MM.", modNames[modOrder[mod]], sep=""))
}

# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.BMI));
geneInfo = geneInfo0[geneOrder, ]

write.table(geneInfo,"geneBMIsignificance-power_modules.txt",col.names=T,row.names=F,sep="\t",quote=F)

#pdf("BMI_vs_darkgrey.pdf")
#plot(covs$BMI,MEs$MEdarkgrey,main="BMI vs Dark grey ME\n cor= 0.54, p < 2.2x10-16",xlab="BMI",ylab="Dark Grey ME")
#abline(lm(MEs$MEdarkgrey ~ covs$BMI))
#dev.off()
#
#pdf("BMI_vs_midnightblue.pdf")
#plot(covs$BMI,MEs$MEmidnightblue,main="BMI vs Midnight blue ME\n cor= -0.62, p=7.87x10-76",xlab="BMI",ylab="Midnight Blue ME")
#abline(lm(MEs$MEmidnightblue ~ covs$BMI))
#dev.off()