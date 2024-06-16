# Use Google API for data downloading peak feature. 
# Please "install.packages('googledrive')" and "install.packages('httpuv')"first.
library(googledrive);
temp <- tempfile(fileext = ".csv")
# Please authorize your google account to access the data
dl3 <- drive_download(
  as_id("1cPs3vZhB1gVCOV3ER9BquMAFSjSvmDYe"), path = temp, overwrite = TRUE)
mSet <- InitDataObjects("pktable", "utils", FALSE)
## we set samples in "row" according to the table format. If your samples are in column, set it as "col".
mSet <- Read.SignalDriftData(mSet, dl3$local_path, "row")
mSet <- PerformSignalDriftCorrection(mSet)
