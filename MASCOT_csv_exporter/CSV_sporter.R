######################################################################
#CSV_SPORTER.R: 
#This script create a xml file with filted data from Mascot csv export 
######################################################################


list.of.packages <- c("readr", "readxl","WriteXLS")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library('readr'); library('readxl'); library('WriteXLS')

#Suppress warnings globally
options(warn = -1)
print("Ejecutando... Puede tardar un poco. ")


#Read files from directory. Take the control from 
files_glob <- (Sys.glob("*.csv")) 
control <- grep('Control', files_glob, value=TRUE)
files_glob_no_control <- files_glob[lapply(files_glob,function(x) length(grep("Control",x,value=FALSE))) == 0]



F001437 <- read_csv("D:/CAI-Proteomics/MASCOT_csv_exporter/F001437.csv", skip = 67, sep=",")
