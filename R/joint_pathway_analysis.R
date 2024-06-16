download.file("https://www.metaboanalyst.ca/MetaboAnalyst/resources/data/integ_genes.txt", "integ_genes.txt", "curl")
download.file("https://www.metaboanalyst.ca/MetaboAnalyst/resources/data/integ_cmpds.txt", "integ_cmpds.txt", "curl")
library(MetaboAnalystR)
# Initiate MetaboAnalyst
mSet<-InitDataObjects("conc", "pathinteg", FALSE)
# Set organism library
mSet<-SetOrganism(mSet, "hsa")
# Set the name of your file containing your gene list
geneListFile<-"integ_genes.txt"
# Read in your gene list file
geneList<-readChar(geneListFile, file.info(geneListFile)$size)
# Perform gene mapping of your file
mSet<-PerformIntegGeneMapping(mSet, geneList, "hsa", "symbol");
# Set the name of your file containing your compound list
cmpdListFile<-"integ_cmpds.txt"
# Read in your compound list file
cmpdList<-readChar(cmpdListFile, file.info(cmpdListFile)$size)
# Perform compound mapping of your file
mSet<-PerformIntegCmpdMapping(mSet, cmpdList, "hsa", "kegg");
# Create a mapping result table
mSet<-CreateMappingResultTable(mSet)
# Prepare data for joint pathway analysis
mSet<-PrepareIntegData(mSet)
#### OPTION 1 ####
# Perform integrated pathway analysis, using hypergeometric test, degree centrality, and the gene-metabolite pathways
# Use the current integrated metabolic pathways
# Saves the output as MetaboAnalyst_result_pathway.csv
mSet<-PerformIntegPathwayAnalysis(mSet, "dc", "hyper", "current", "integ", "query");
# View the output of the pathway analysis
mSet$dataSet$path.mat
#### OPTION 2 ####
# Perform integrated pathway analysis using the current library of all pathways (not just metabolic)
# Fisher's exact test for enrichment and betweenness centrality for topology measure
mSet<-PerformIntegPathwayAnalysis(mSet, "bc", "fisher", "current", "all", "query")
# Perform pathway analysis
mSet<-PlotPathSummary(mSet, "path_view_0_", "png", 72, width=NA)
PreparePDFReport(mSet, "My Name")
