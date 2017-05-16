######################################################################
#TEST_ACCESION: 
#This script create a xml file with merged data and Venn Diagram plot. 
#IMPORTANT: CHANGE NAME OF THE FILE EACH TIME THAT YOU WANT TO USE IT. 
######################################################################

####################################################################################
#REPLACE NAME OF THE FILE: // REEMPALZAR LOS NOMBRES DE LOS xlsx QUE ESTAN ENTRE COMILLAS
# File1 <- file.path(getwd(), "17-17_Resultados ProteinPilot SPHuman Muestra C.xlsx")
# File2 <- file.path(getwd(), "17-17_Resultados ProteinPilot SPHuman Muestra X.xlsx")
# File3 <- file.path(getwd(), "17-17_Resultados ProteinPilot SPHuman Muestra Tg.xlsx") 
# File4 <- file.path(getwd(),"17-17_Resultados ProteinPilot SPHuman Muestra Fl.xlsx")

####################################################################################

list.of.packages <- c("VennDiagram", "readxl","WriteXLS")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)


# if (!require('VennDiagram')) install.packages('VennDiagram')
 library('VennDiagram')
# if (!require('readxl')) install.packages('readxl')
 library('readxl')
# if (!require('WriteXLS')) install.packages('WriteXLS')
 library('WriteXLS')

#Suppress warnings globally
options(warn = -1)
cat("Ejecutando... Puede tardar un poco... \n\n")


#Read files from directory. 
files_glob <- (Sys.glob("*.xlsx")) 
control <- grep('Control', files_glob, value=TRUE)
files_glob_no_control <- files_glob[files_glob != control]


repeat{
  cat("�Quieres hacer la comparativa con alg�n control?\n ")
  cat("1. Yes | 2. No\n")
  x <- readLines(file("stdin"),1)#enter "yes"
  if (x=="yes" | x==1 | x=="Yes")  {
    #Bucle if primero:
    

    
    
    for(i in 1:length(files_glob_no_control)){
      
      df1 <- read_excel(file.path(getwd(), grep('Control', files_glob, value=TRUE)))
      df2 <- read_excel(file.path(getwd(), files_glob_no_control[i]))
      
      #Combine data frame using reduce function
      df_final <- Reduce(function(x, y) merge (x, y, by = c("Name", "Accession"), all = TRUE), list(df1, df2))
      
      
      #Export data frame to table.
      WriteXLS(df_final, ExcelFileName = paste(gsub("\\.xlsx*","",files_glob_no_control[i]), "_Multiconsenso_Control.xlsx", sep = "", na=""), SheetNames = NULL, perl = "perl",
               verbose = FALSE, Encoding = c("UTF-8", "latin1", "cp1252"),
               row.names = FALSE, col.names = TRUE,
               AdjWidth = FALSE, AutoFilter = FALSE, BoldHeaderRow = FALSE,
               na = "",
               FreezeRow = 0, FreezeCol = 0,
               envir = parent.frame())
      
      
      #Create a Venn Diagram: 
      #Dissable .log files
      futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
      
      #Function ven.diagram and grid. 
      tiff( width=5, height=5, units="in",
            pointsize=8, compression="lzw", bg="white", res=600,
            restoreConsole=TRUE,paste(gsub("\\.xlsx*","",files_glob_no_control[i]), "_VenDiagram_Multiconsenso_Control.tiff", sep = "", na=""))
      
      
      v <- venn.diagram(list(Control=df1$Accession,Muestra=df2$Accession),
                        fill = c("red", "blue"),
                        cat.cex = 1.5, cex=1.5,cat.pos=0,
                        filename=NULL)
      
      
      
      
      # have a look at the default plot
      grid.newpage()
      grid.draw(v)
      garbage <- dev.off() 
    }
    
    
    break

  }else if(x=="no" | x==2 | x=="No"){
    
    if (length(files_glob_no_control) == 1){ 
      print("There is only 1 xlsx file")
      
    } else if (length(files_glob_no_control) == 2){
      #Load excel file
      
      df1 <- read_excel(file.path(getwd(), files_glob_no_control[1]))
      df2 <- read_excel(file.path(getwd(), files_glob_no_control[2]))
      
      
      
      #Combine data frame using reduce function
      df_final <- Reduce(function(x, y) merge (x, y, by = c("Name", "Accession"), all = TRUE), list(df1, df2))
      
      #Export data frame to table.
      WriteXLS(df_final, ExcelFileName = "Multiconsenso_2.xls", SheetNames = NULL, perl = "perl",
               verbose = FALSE, Encoding = c("UTF-8", "latin1", "cp1252"),
               row.names = FALSE, col.names = TRUE,
               AdjWidth = FALSE, AutoFilter = FALSE, BoldHeaderRow = FALSE,
               na = "",
               FreezeRow = 0, FreezeCol = 0,
               envir = parent.frame())
      

      
      
      
      #Create a Venn Diagram: 
      #Dissable .log files
      futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
      
      #Function ven.diagram and grid. 
      pdf("Venn_Diagram_2.pdf")
      
      
      v <- venn.diagram(list(Muestra1=df1$Accession,Muestra2=df2$Accession),
                        fill = c("red", "blue"),
                        cat.cex = 1.5, cex=1.5,
                        filename=NULL)
      # have a look at the default plot
      grid.newpage()
      grid.draw(v)
      garbage <- dev.off()
      
      
      
      
      
    } else if (length(files_glob_no_control) == 3){
      #Load excel file
      df1 <- read_excel(file.path(getwd(), files_glob_no_control[1]))
      df2 <- read_excel(file.path(getwd(), files_glob_no_control[2]))
      df3 <- read_excel(file.path(getwd(), files_glob_no_control[3]))
      
      #Combine data frame using reduce function
      df_final <- Reduce(function(x, y) merge (x, y, by = c("Name", "Accession"), all = TRUE), list(df1, df2, df3))
      
      #Export data frame to table.
      WriteXLS(df_final, ExcelFileName = "Multiconsenso_3.xls", SheetNames = NULL, perl = "perl", 
               verbose = FALSE, Encoding = c("UTF-8", "latin1", "cp1252"),
               row.names = FALSE, col.names = TRUE,
               AdjWidth = FALSE, AutoFilter = FALSE, BoldHeaderRow = FALSE,
               na = "",
               FreezeRow = 0, FreezeCol = 0,
               envir = parent.frame())
      
      
      #Create a Venn Diagram: 
      #Dissable .log files
      futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
      
      #Function ven.diagram and grid. 
      pdf("Venn_Diagram_3.pdf")
      v <- venn.diagram(list(Muestra1=df1$Accession,Muestra2=df2$Accession, Muestra3=df3$Accession),
                        fill = c("red", "blue", "green"),
                        cat.cex = 1.5, cex=1.5,
                        filename=NULL)
      # have a look at the default plot
      grid.newpage()
      grid.draw(v)
      garbage <- dev.off()
      
    } else if (length(files_glob_no_control) == 4){
      df1 <- read_excel(file.path(getwd(), files_glob_no_control[1]))
      df2 <- read_excel(file.path(getwd(), files_glob_no_control[2]))
      df3 <- read_excel(file.path(getwd(), files_glob_no_control[3]))
      df4 <- read_excel(file.path(getwd(), files_glob_no_control[4]))
      
      
      #Combine data frame using reduce function
      df_final <- Reduce(function(x, y) merge (x, y, by = c("Name", "Accession"), all = TRUE), list(df1, df2, df3, df4))
      
      #Export data frame to table.
      WriteXLS(df_final, ExcelFileName = "Multiconsenso_4.xls", SheetNames = NULL, perl = "perl",
               verbose = FALSE, Encoding = c("UTF-8", "latin1", "cp1252"),
               row.names = FALSE, col.names = TRUE,
               AdjWidth = FALSE, AutoFilter = FALSE, BoldHeaderRow = FALSE,
               na = "",
               FreezeRow = 0, FreezeCol = 0,
               envir = parent.frame())
      
      
      #Create a Venn Diagram: 
      #Dissable .log files
      futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
      
      #Function ven.diagram and grid. 
      pdf("Venn_Diagram_4.pdf")
      v <- venn.diagram(list(Muestra1=df1$Accession,Muestra2=df2$Accession,Muestra3=df3$Accession,Muestra4=df4$Accession),
                        fill = c("red", "blue", "green", "yellow"),
                        cat.cex = 1.5, cex=1.5,
                        filename=NULL)
      # have a look at the default plot
      grid.newpage()
      grid.draw(v)
      garbage <- dev.off()
      
      
      
    }else{
      df1 <- read_excel(file.path(getwd(), files_glob_no_control[1]))
      df2 <- read_excel(file.path(getwd(), files_glob_no_control[2]))
      df3 <- read_excel(file.path(getwd(), files_glob_no_control[3]))
      df4 <- read_excel(file.path(getwd(), files_glob_no_control[4]))
      df5 <- read_excel(file.path(getwd(), files_glob_no_control[5]))
      
      #Load excel file
      df1 <- read_excel(file.path(File1))
      df2 <- read_excel(file.path(File2))
      df3 <- read_excel(file.path(File3))
      df4 <- read_excel(file.path(File4))
      df5 <- read_excel(file.path(File5))
      
      #Combine data frame using reduce function
      df_final <- Reduce(function(x, y) merge (x, y, by = c("Name", "Accession"), all = TRUE), list(df1, df2, df3, df4, df5))
      
      #Export data frame to table.
      WriteXLS(df_final, ExcelFileName = "Multiconsenso_5.xls", SheetNames = NULL, perl = "perl",
               verbose = FALSE, Encoding = c("UTF-8", "latin1", "cp1252"),
               row.names = FALSE, col.names = TRUE,
               AdjWidth = FALSE, AutoFilter = FALSE, BoldHeaderRow = FALSE,
               na = "",
               FreezeRow = 0, FreezeCol = 0,
               envir = parent.frame())
      
      #Create a Venn Diagram: 
      #Dissable .log files
      futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
      
      #Function ven.diagram and grid. 
      pdf("Venn_Diagram_5.pdf")
      v <- venn.diagram(list(Muestra1=df1$Accession,Muestra2=df2$Accession,Muestra3=df3$Accession,Muestra4=df4$Accession, Muestra5=df5$Accession),
                        fill = c("red", "blue", "green", "yellow", "black"),
                        cat.cex = 1.5, cex=1.5,
                        filename=NULL)
      # have a look at the default plot
      grid.newpage()
      grid.draw(v)
      garbage <- dev.off()
      
      
    }
    
    
    
    break
  }else{
    cat("No han incluido nada\n")
  }
}



