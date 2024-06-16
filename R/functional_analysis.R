# This example data set is a peak list from a COVID-19 study
# Two groups (COVID-19 vs. Healthy Controls) are included for this example

download.file("https://raw.githubusercontent.com/Zhiqiang-PANG/MetaboRaw/master/examples/peaks_ms1.txt", 
              destfile = "peaks_ms1.txt", 
              mode = "auto")

# Load MetaboAnalystR
library(MetaboAnalystR)

# Clean global environment
rm(list = ls())

# Create objects for storing processed data
mSet <- InitDataObjects("mass_all", "mummichog", FALSE)

# Set Peak format as "mpt" here.
# This option can be "mp", "mt", "mpt", "mprt", "mrt", "mpr", "rmp", "rmt".
# "rmp" and "rmt" refers to peaks ranked by p values or t scores;
# For other options, "m" means "mz"; "p" means "p value"; "t" means "t score"; "r" means "retention time";
mSet <- SetPeakFormat(mSet, "mpt")

# Set parameters for analysis, in this case the mass accuracy is set to 15 ppm, 
# the mode of the MS instrument is "mixed" (contains both positive and negative), 
# and the p-value cut-off is 0.02
mSet <- UpdateInstrumentParameters(mSet, 15.0, "mixed", "yes", 0.02)

# Read in peak-list data
mSet <- Read.PeakListData(mSet, "peaks_ms1.txt");

# set retention time included, unit is "seconds"
mSet <- SetRTincluded(mSet, "seconds")

# Sanity check of the uploaded data
mSet <- SanityCheckMummichogData(mSet)

# set peak enrichment method. Method can be one of the "mum", "gsea" or "integ";
# Method, "integ" means integration of both mummichog and GSEA algorithm;
# version can be "v1" or "v2" ("v1" will use m/z only; "v2" will use both "m/z" and "retention time")
mSet <- SetPeakEnrichMethod(mSet, "mum", "v2")

# Here we use the top 10% peaks as the p value cutoff
pval <- sort(mSet[["dataSet"]][["mummi.proc"]][["p.value"]])[ceiling(length(mSet[["dataSet"]][["mummi.proc"]][["p.value"]])*0.1)]
mSet <- SetMummichogPval(mSet, pval)

# Perform the mummichog algorith, in this case the model is the human MFN model, using "current" version by default
# This function may take sometime for processing, and will output the pathway-results and the compound matching tables in your working directory
mSet <- PerformPSEA(mSet, "hsa_mfn", "current", 3 , 100)

# To view the results of the pathway analysis
mSet <- PlotPeaks2Paths(mSet, "peaks_to_paths_ms1_", "png", 72, width=8)

## In this case study, we are downloading a peak table as an example dataset
## This peak table is from a Malaria study, which includes two immune status (Naive vs. Semi-immune)
## Six samples are included in each group
download.file("https://raw.githubusercontent.com/Zhiqiang-PANG/MetaboRaw/master/examples/malaria_feature_table.csv", 
              destfile = "malaria_feature_table.csv", mode = "auto")

# Load MetaboAnalystR
library(MetaboAnalystR)

# Clean global environment
rm(list = ls())

# Create objects for storing processed data
mSet <- InitDataObjects("mass_table", "mummichog", FALSE)

# Set parameters, ppm is 4.3 here
# Only positive mode (ESI+) included
mSet <- SetPeakFormat(mSet, "colu")
mSet <- UpdateInstrumentParameters(mSet, 4.3, "positive", "yes", 0.02);
mSet <- SetRTincluded(mSet, "seconds")

# Read Peak table
mSet <- Read.TextData(mSet, "malaria_feature_table.csv", "mpt", "disc");
mSet <- SanityCheckMummichogData(mSet)

# Replace minimum value and data filtration
mSet <- ReplaceMin(mSet);
mSet <- SanityCheckMummichogData(mSet)

mSet <- FilterVariable(mSet, "none", -1, "F", 25, F)

# Perform data normalization
# The normailization distribution will be displaced in figure "norm_0_dpi72.png" and "snorm_0_dpi72.png"
# in your working directory
mSet <- PreparePrenormData(mSet)
mSet <- Normalization(mSet, "NULL", "LogNorm", "NULL", ratio=FALSE, ratioNum=20)

mSet <- PlotNormSummary(mSet, "norm_0_", "png", 72, width=NA)
mSet <- PlotSampleNormSummary(mSet, "snorm_0_", "png", 72, width=NA)

# Perform functional analysis with mummichog algorithm
mSet <- SetPeakEnrichMethod(mSet, "mum", "v2")
mSet <- PreparePeakTable4PSEA(mSet)

mSet <- SetMummichogPval(mSet, 0.003)
mSet <- PerformPSEA(mSet, "hsa_mfn", "current", 3 , 100)

mSet <- PlotPeaks2Paths(mSet, "peaks_to_paths_0_", "png", 72, width=NA)

## In this case study, we are downloading a peak list and the corresponding compound list as an example dataset
## This dataset is from a COVID-19 study, which includes two groups (COVID vs. Healthy Controls)
## Six samples are included in each group
## This example is same as "Case Study 1" above
download.file("https://raw.githubusercontent.com/Zhiqiang-PANG/MetaboRaw/master/examples/peaks_ms1.txt", 
              destfile = "peaks_ms1.txt", mode = "auto")

download.file("https://raw.githubusercontent.com/Zhiqiang-PANG/MetaboRaw/master/examples/compound_ms2.txt", 
              destfile = "compound_ms2.txt", mode = "auto")

# Load MetaboAnalystR
library(MetaboAnalystR)

# Clean global environment
rm(list = ls())

# Create objects for storing processed data
mSet <- InitDataObjects("mass_all", "mummichog", FALSE)

# Set parameters, ppm is 4.3 here
# Only positive mode (ESI+) included
mSet <- SetPeakFormat(mSet, "mpt")
mSet <- SetMS2IDType(mSet, "inchikeys")
mSet <- UpdateInstrumentParameters(mSet, 15.0, "mixed", "yes", 0.02)

# Read in MS peaks and MS/MS based compounds 
mSet <- Read.PeakMS2ListData(mSet, msfile = "peaks_ms1.txt", 
                             msmsfile = "compound_ms2.txt")

# Set parameters
# Set Retention time unit as "seconds"
# set algorithm as "mummichog" and version 2
# Here we use the top 10% peaks as the p value cutoff
mSet <- SetRTincluded(mSet, "seconds")
mSet <- SanityCheckMummichogData(mSet)
mSet<-SetPeakEnrichMethod(mSet, "mum", "v2")
pval <- sort(mSet[["dataSet"]][["mummi.proc"]][["p.value"]])[ceiling(length(mSet[["dataSet"]][["mummi.proc"]][["p.value"]])*0.1)]
mSet<-SetMummichogPval(mSet, pval)

# Perform functional analysis
mSet<-PerformPSEA(mSet, "hsa_mfn", "current", 3 , 100)

mSet<-PlotPeaks2Paths(mSet, "peaks_to_paths_msms_", "png", 72, width=8)

# Load MetaboAnalystR
library(MetaboAnalystR)

# Clean global environment
rm(list = ls())

# Create objects for storing processed data
mSet <- InitDataObjects("mass_all", "mummichog", FALSE)

# Set peak formart - contains m/z features, p-values and t-scores
mSet <- SetPeakFormat(mSet, "mpt")
mSet <- UpdateInstrumentParameters(mSet, 5.0, "mixed");
mSet <- Read.PeakListData(mSet, "https://www.metaboanalyst.ca/MetaboAnalyst/resources/data/mummichog_mixed.txt");
mSet <- SanityCheckMummichogData(mSet)

# Customize currency
curr.vec <- c("Water (C00001)", "Proton (C00080)", "Oxygen (C00007)", "NADPH (C00005)",
              "NADP (C00006)", "NADH (C00004)", "NAD (C00003)", "Adenosine triphosphate (C00002)",
              "Pyrophosphate (C00013)","Phosphate (C00009)","Carbon dioxide (C00011)",
              "Hydrogen (C00282)","Hydrogen peroxide (C00027)","Sodium (C01330)");

# Map selected currency to user's data
mSet <- Setup.MapData(mSet, curr.vec)
mSet <- PerformCurrencyMapping(mSet)

# Now customize adducts
add.vec <- c("M [1+]","M+H [1+]","M+2H [2+]","M+3H [3+]","M+Na [1+]",
             "M+H+Na [2+]","M+K [1+]","M+H2O+H [1+]","M-H2O+H [1+]",
             "M-H4O2+H [1+]","M(C13)+H [1+]","M(C13)+2H [2+]","M(C13)+3H [3+]",
             "M(Cl37)+H [1+]","M-NH3+H [1+]","M-CO+H [1+]","M-CO2+H [1+]",
             "M-HCOOH+H [1+]","M+HCOONa [1+]","M-HCOONa+H [1+]","M+NaCl [1+]",
             "M-C3H4O2+H [1+]","M+HCOOK [1+]","M-HCOOK+H [1+]","M-H [1-]",
             "M-2H [2-]","M-H2O-H [1-]","M-H+O [1-]","M+K-2H [1-]","M+Na-2H [1- ]",
             "M+Cl [1-]","M+Cl37 [1-]","M+HCOO [1-]","M+CH3COO [1-]")

# Set up the selected adducts
mSet <- Setup.AdductData(mSet, add.vec)
mSet <- PerformAdductMapping(mSet, "mixed")

# Perform mummichog algorithm using selected currency and adducts, using Version1 of the mummichog algorithm
mSet <- SetPeakEnrichMethod(mSet, "mum", "v1")
mSet <- SetMummichogPval(mSet, 1.0E-5)
mSet <- PerformPSEA(mSet, "hsa_mfn", "current", 100)

# Image is not shown below to avoid to be overwhelmed.
mSet <- PlotPeaks2Paths(mSet, "peaks_to_paths_2_", "png", 72, width=NA)

# Disable force primary ion
mSet <- UpdateInstrumentParameters(mSet, instrumentOpt, msModeOpt, force_primary_ion ="no", rt_frac =0.02,rt_tol =NA)
# Change retention time fraction when calculating the retention time window
mSet <- UpdateInstrumentParameters(mSet, instrumentOpt, msModeOpt,force_primary_ion ="yes", rt_frac =0.025, rt_tol =NA)
# Set the retention time window (in seconds)
mSet <- UpdateInstrumentParameters(mSet, instrumentOpt, msModeOpt, force_primary_ion ="yes", rt_frac =0.02, rt_tol =25)

# Load MetaboAnalystR
library(MetaboAnalystR)

# Clean global environment
rm(list = ls())

# Create objects for storing processed data

mSet<-InitDataObjects("mass_all", "mummichog", FALSE)

mSet<-SetPeakFormat(mSet, "mpt")
mSet<-UpdateInstrumentParameters(mSet, 5.0, "negative", "yes", 0.02);
mSet<-Read.PeakListData(mSet, "https://www.metaboanalyst.ca/MetaboAnalyst/resources/data/mummichog_ibd.txt");
mSet<-SetRTincluded(mSet, "no")
mSet<-SanityCheckMummichogData(mSet)

## Set method as "GSEA"
mSet<-SetPeakEnrichMethod(mSet, "gsea", "v2")
mSet<-PerformPSEA(mSet, "hsa_mfn", "current", 3 , 100)

mSet<-PlotPeaks2Paths(mSet, "peaks_to_paths_3_", "png", 72, width=NA)
