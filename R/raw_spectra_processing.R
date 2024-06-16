# Load MetaboAnalystR
library(MetaboAnalystR)
# Load OptiLCMS
library(OptiLCMS)

download.file("https://rest.xialab.ca/api/download/metaboanalyst/malaria_r_example.zip",
              destfile = "malaria_raw.zip",
              method = "curl")
unzip("malaria_raw.zip", exdir = "upload")

# Load MetaboAnalystR and OptiLCMS
library(MetaboAnalystR)
library(OptiLCMS)

# Here, we extract ROIs from 3 QC samples.
DataFiles <- list.files("upload/QC/", full.names = TRUE)
mSet <- PerformROIExtraction(datapath = DataFiles, rt.idx = 0.9, rmConts = TRUE)

# Here we use PerformParamsOptimization to optimize parameters based on 
# the extracted ROI (stored in 'mSet') before process the entire dataset

best_params <- PerformParamsOptimization(mSet, param = SetPeakParam(platform = "UPLC-Q/E"), ncore = 4)

# "path" is used to specify the path to the folder containing the raw MS spectra to be processed.
# BPI and TIC plotting can be enabled with parameter, 
# "plotSettings = SetPlotParam(Plot = T)", or disabled by changing "T" into "F";

mSet <- ImportRawMSData(path = c("upload"), plotSettings = SetPlotParam(Plot = T))

# "mSet" include complete raw MS spectra to be processed.
# "params" is using the "best_params" generated above
# Plotting functions can be enabled with parameter, 
# "plotSettings = SetPlotParam(Plot = T)", or disabled by changing "T" into "F";
mSet <- PerformPeakProfiling(mSet, Params = best_params, plotSettings = SetPlotParam(Plot=TRUE))

# We firstly define the parameters for feature annotation

# 'polarity' is required, can be either 'negative' or 'positive';
# 'perc_fwhm' is used to set the percentage of the width of the FWHM for peak grouping. 
#              Default is set to 0.6;
# 'mz_abs_iso' is used to set the allowed variance for the search (for isotope annotation). 
#              The default is set to 0.005;
# 'max_charge' is set the maximum number of the isotope charge. 
#              For example, the default is 2, therefore the max isotope charge is 2+/-;
# 'max_iso' is used to set the maximum number of isotope peaks.
#              For example, the default is 2, therefore the max number of isotopes per peaks is 2;
# 'corr_eic_th' is used to set the threshold for intensity correlations across samples. 
#              Default is set to 0.85.
# 'mz_abs_add' is used to set the allowed variance for the search (for adduct annotation). 
#              Default is set to 0.001.
# 'adducts' is used to specify the adducts based on your instrument settings.

annParams <- SetAnnotationParam(polarity = 'positive', mz_abs_add = 0.015);

# "mSet" include processed raw MS spectra to be processed.
# "annParams" is the parameters used for annotation

mSet <- PerformPeakAnnotation(mSet, annParams)

# Here we format and filter the peak list for following analysis with MetaboAnalystR

# Parameters are explained as below,
# annParams, is the object created using the SetAnnotationParam function above;
# filtIso, is used to decide to filter out all isotopes (TRUE) or not (FALSE);
# filtAdducts, is used to decide to filter out all adducts (TRUE) or not (FALSE);
# missPercent, specify the threshold to remove features missing in a certain percentage
#              of all samples in a group.

mSet <- FormatPeakList(mSet, annParams, filtIso = FALSE, filtAdducts = FALSE, missPercent = 1)

# Export annotation results, the annotation will be save as "annotated_peaklist.csv";
Export.Annotation(mSet)

# Export complete feature table. It will be saved as "metaboanalyst_input.csv";
# This table can be used for statistic analysis, functional analysis, biomarker analysis module directly.
Export.PeakTable(mSet)

# Export a summary table (peak_result_summary.txt) to summarize the information of all peaks in
# different samples. The columns are sample names, groups, retention time range, m/z range of all peaks,
# number of all peaks and percentage of missing features.
Export.PeakSummary(mSet)

# There are 6 raw spectral files (.mzML) and 1 MS feature table included in the downloaded folder
# The MS feature table was generated with the MS1 data processing pipeline (as described in section 2)

download.file("https://rest.xialab.ca/api/download/metaboanalyst/ms2_dda_example.zip",
              destfile = "ms2_dda_example.zip",
              method = "curl")
unzip("ms2_dda_example.zip", exdir = "ms2_dda")

## Here we are downloading biology MS/MS reference database.

download.file("https://rest.xialab.ca/api/download/metaboanalyst/MS2ID_Bio_v09102023.zip",
              destfile = "MS2ID_Bio.zip",
              method = "curl")
unzip("MS2ID_Bio.zip", exdir = "ms2_db")

## We load MetaboAnalystR and OptiLCMS
library(MetaboAnalystR)
library(OptiLCMS)

## Clean the environment first
rm(list = ls())

## Read the MS1 feature table as the target features for MS/MS data processing
## This table include four columns (mzmin, mzmax, rtmin, rtmax)
## mzmin and mzmax is the minimum and the maximum value of m/z for the feature;
## rtmin and rtmax is the minimum and the maximum value of retention time for the feature;
ft_dt <- qs::qread("ms2_dda/ms2_dda_example/ft_dt.qs")

## Here we use function PerformMSnImport to read MS/MS data
## This step may take seconds to minutes depending on the size of your dataset
mSet <- PerformMSnImport(filesPath = c(list.files("ms2_dda/ms2_dda_example/",
                                                  pattern = ".mzML",
                                                  full.names = T, recursive = T)),
                         targetFeatures = ft_dt,
                         acquisitionMode = "DDA")

# In this step, we are processing the DDA spectra deconvolution;

# The deconvolution is based on the Biology MS/MS spectral database;
# Parallel computing is supported, users are encouraged to use multiple cores to speed up;
# This step may take minutes to hours to finish, depending on the size of datasets

system.time(mSet <- PerformDDADeconvolution(mSet,
                                            ppm1 = 5,
                                            ppm2 = 10,
                                            sn = 12,
                                            filtering = 0,
                                            window_size = 1.5,
                                            intensity_thresh = 1.6e5,
                                            database_path = "ms2_db/MS2ID_Bio_v09102023.sqlite",
                                            ncores = 6L))

mSet <- PerformSpectrumConsenus (mSet,
                                 ppm2 = 15,
                                 concensus_fraction = 0.2,
                                 database_path = "",
                                 use_rt = FALSE,
                                 user_dbCorrection = FALSE)

# PerformDBSearchingBatch is used to seatching MS/MS reference library
# Results will be scored based on the similarity rules above
# Parallel computing is allowed. CPU cores used are controlled by argument "ncores".

mSet <- PerformDBSearchingBatch (mSet,
                                 ppm1 = 10,
                                 ppm2 = 25,
                                 rt_tol = 5,
                                 database_path = "ms2_db/MS2ID_Bio_v09102023.sqlite",
                                 use_rt = FALSE,
                                 enableNL = FALSE,
                                 ncores = 6L)

mSet <- PerformResultsExport (mSet,
                              type = 0L,
                              topN = 10L,
                              ncores = 3L)

# 2nd argument, TopN, is the argument used to specifiy the number of compounds to export
# If the dataset is a lipidomics dataset, please set 3rd argument as "TRUE"
# to extract the lipid classification information
# All identified compounds have been summarized as a data.frame (dtx) here

dtx <- FormatMSnAnnotation(mSet, 5L, F)

## Here we are downloading biology MS/MS reference database.

download.file("https://rest.xialab.ca/api/download/metaboanalyst/FragsAnnotateDB_v02042024.zip",
              destfile = "FragsAnnotateDB.zip",
              method = "curl")
unzip("FragsAnnotateDB.zip", exdir = "ms2_db")

save(mSet, file = "mSet_raw_dda.rda") # you are suggested to save this, in case of lossing data

# Here, user need to convert raw spec mSet object into regular analysis object;
# Several parameters need to be specified
# 1. peak_idx: Index of the peak. For example, 11 refers to the 11th peak in the target peak list;
# 2. sub_idx: Index of the identification result. For example, 1 refers to the first identified compound (the one with the highest matching score.)
# 3. interative: can be TRUE or FALSE. If TRUE, will plot an interactive figure.
# 4. ppm, the parameter used to determine the mactching results.
mSet <- Convert2AnalObject(mSet, "raw", "spec", F)

mSet <- PerformMirrorPlotting(mSetObj = mSet, 
                              fragDB_path = "ms2_db/FragsAnnotateDB_v02042024.sqlite", 
                              peak_idx = 11, sub_idx = 1, 
                              interactive = T, ppm = 25, 
                              dpi = 72, format = "png", width = 8, height = 6)

# There are 6 raw spectral files (.mzML) and 1 text file, which include the design of SWATH windows

download.file("https://rest.xialab.ca/api/download/metaboanalyst/ms2_dia_example.zip",
              destfile = "ms2_dia_example.zip",
              method = "curl")

unzip("ms2_dia_example.zip", exdir = "ms2_dia")

## Load OptiLCMS
library(OptiLCMS)
rm(list = ls())

## Construct meta data
meta_dt <- data.frame(samples = c("210210-SWATH-NEG-Covid-Cov-16.mzML", 
                                  "210210-SWATH-NEG-Covid-Cov-17.mzML", 
                                  "210210-SWATH-NEG-Covid-Cov-18.mzML",
                                  "210210-SWATH-NEG-Covid-Ct-1.mzML",
                                  "210210-SWATH-NEG-Covid-Ct-2.mzML",
                                  "210210-SWATH-NEG-Covid-Ct-3.mzML"),
                      groups = c(rep("COVID",3), rep("Control",3)))

## Import raw MS1 data
mSet1 <- ImportRawMSData(path = "ms2_dia/ms2_dia_example/", metadata = meta_dt)

## Process MS1 peaks by using customized parameters
mSet1 <- PerformPeakProfiling(mSet1, 
                              Params = SetPeakParam(ppm = 25,
                                                    bw = 3,
                                                    mzdiff = 0.001,
                                                    max_peakwidth = 35,
                                                    min_peakwidth = 5,
                                                    noise = 200, 
                                                    minFraction = 0.2),
                              ncore = 6,
                              plotSettings = SetPlotParam(Plot = F))

## Annotation
annParams <- SetAnnotationParam(polarity = 'negative',
                                mz_abs_add = 0.035);

mSet1 <- PerformPeakAnnotation(mSet = mSet1,
                               annotaParam = annParams,
                               ncore =1)

## Format MS1 feature table
mSet1 <- FormatPeakList(mSet = mSet1,
                        annParams,
                        filtIso =FALSE,
                        filtAdducts = FALSE,
                        missPercent = 1)

mSet <- PerformMSnImport(mSet = mSet1,
                         filesPath = c(list.files("ms2_dia/ms2_dia_example/", pattern = ".mzML", full.names = T, recursive = T)),
                         acquisitionMode = "DIA",
                         SWATH_file = "ms2_dia/ms2_dia_example/DIA_SWATH_MS_experiment_file_neg.txt")

mSet <- PerformDIADeconvolution(mSet,
                                min_width = 5,
                                span = 0.3,
                                ppm2 = 30,
                                sn = 12,
                                filtering = 0,
                                ncores = 6L)

mSet <- PerformSpectrumConsenus (mSet,
                                 ppm2 = 30,
                                 concensus_fraction = 0.25,
                                 database_path = "",
                                 use_rt = FALSE,
                                 user_dbCorrection = FALSE)

mSet <- PerformDBSearchingBatch (mSet,
                                 ppm1 = 15,
                                 ppm2 = 30,
                                 rt_tol = 5,
                                 database_path = "ms2_db/MS2ID_Bio_v09102023.sqlite",
                                 use_rt = FALSE,
                                 enableNL = FALSE,
                                 ncores = 4L)

mSet <- PerformResultsExport (mSet,
                              type = 0L,
                              topN = 10L,
                              ncores = 4L)

dtx2 <- FormatMSnAnnotation(mSet, 5L, F)

save(mSet, file = "mSet_raw_dia.rda")
mSet <- Convert2AnalObject(mSet, "raw", "spec", F)

mSet <- PerformMirrorPlotting(mSetObj = mSet, 
                              fragDB_path = "ms2_db/FragsAnnotateDB_v02042024.sqlite", 
                              peak_idx = 5, sub_idx = 1, 
                              interactive = T, ppm = 25, 
                              dpi = 72, format = "png", width = 8, height = 6)
