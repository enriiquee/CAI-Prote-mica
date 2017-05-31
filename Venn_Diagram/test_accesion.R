######################################################################
#TEST_ACCESION: 
#This script create a xml file with merged data, comparate and create
#Venn Diagram plot. 
#IMPORTANT: CHANGE NAME OF THE FILE EACH TIME THAT YOU WANT TO USE IT. 
######################################################################

#Packages
list.of.packages <- c("VennDiagram", "readxl","WriteXLS", "dplyr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

#Load packages
# if (!require('VennDiagram')) install.packages('VennDiagram')
 library('VennDiagram')
# if (!require('readxl')) install.packages('readxl')
 library('readxl')
# if (!require('WriteXLS')) install.packages('WriteXLS')
 library('WriteXLS'); library("dplyr")

#Suppress warnings globally. Avoid error printed
options(warn = -1)
cat("Ejecutando... Puede tardar un poco... \n\n")

#Read files from directory. Take the control from 
files_glob <- (Sys.glob("*.xlsx")) 
control <- grep('Control', files_glob, value=TRUE)
files_glob_no_control <- files_glob[lapply(files_glob,function(x) length(grep("Control",x,value=FALSE))) == 0]

#control2 <- files_glob[lapply(files_glob,function(x) length(grep("Control",x,value=FALSE))) != 0]


#Check how many controls there are. 
if(length(control2)>1){
  print("Hay 2 o más archivos con la palabra Control")
  break
  
}else{
  
  repeat{
    cat("¿Quieres hacer la comparativa con algún control?\n ")
    cat("1. Yes | 2. No\n")
    x <- readLines(file("stdin"),1)#enter "yes"
    x <- "yes"
    if (x=="yes" | x==1 | x=="Yes")  {
      
      #Bucle if primero:
      dfList <- list()
      
      for(i in 1:length(files_glob_no_control)){
        
        df1 <- read_excel(file.path(getwd(), grep('Control', files_glob, value=TRUE)))
        df2 <- read_excel(file.path(getwd(), files_glob_no_control[i]))
        
        #Combine data frame using reduce function
        df_final <- Reduce(function(x, y) merge (x, y, by = c("Name", "Accession"), all = TRUE), list(df1, df2))
        
        #Clasification. Se filtran los resultados
        #Filter with the same peptides
        test1 <- filter(df_final,  !is.na(df_final$`Peptides(95%).x` & df_final$`Peptides(95%).y`))
        test1 <- with(test1,  test1[order(-test1$`Peptides(95%).x`) , ])
        #Different peptides between 1º and the 2º. 
        test2 <- filter(df_final,  !is.na(df_final$`Peptides(95%).x`) & is.na(df_final$`Peptides(95%).y`))
        test2 <- with(test2,  test2[order(-test2$`Peptides(95%).x`) , ])
        #Different peptides between 2º and the 1º.
        test3 <- filter(df_final,  is.na(df_final$`Peptides(95%).x`) & !is.na(df_final$`Peptides(95%).y`))
        test3 <- with(test3,  test3[order(-test3$`Peptides(95%).y`) , ])
        
        #Merge data. 
        test_final<-rbind(test1, test2, test3)
        
        #Añadimos columna N 
        test_final2 <- data.frame(cbind(N = 1:nrow(test_final), test_final))
        test_final2$N.x <- NULL
        test_final2$N.y <- NULL; 
        
        #Eliminales la columna X__1
        if("X__1"%in% colnames(test_final2)){
          test_final2 <- subset(test_final2, select = -c(X__1) )
          colnames(test_final2) <- c("N","Name","Accesion","Score","%Cov(95)","Peptides(95%)","Species","Score","%Cov(95)","Peptides(95%)","Species")
          
        }
        
        #Cambiamos nombres
        colnames(test_final2) <- c("N","Name","Accesion","Score","%Cov(95)","Peptides(95%)","Species","Score","%Cov(95)","Peptides(95%)","Species")
        
        #Round: 
        round_df <- function(df, digits) {
          nums <- vapply(df, is.numeric, FUN.VALUE = logical(1))
          
          df[,nums] <- round(df[,nums], digits = digits)
          
          (df)
        }
        #Round
        test_final3 <- round_df(test_final2, digits=2)
        
        #Save in a list
        dfList[[i]] = data.frame(test_final3)
        
        names(dfList)<-sprintf(paste(gsub(".*","",files_glob_no_control[1:length(dfList)]), "Sumary_vs_Control",1:length(dfList), sep = "", na=""),1:length(dfList))
      

        #Export data frame to table.
        WriteXLS(dfList, ExcelFileName = "Summaries_vs_Control.xlsx", names(dfList))
      
   

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
        
         #Clasification. Se filtran los resultados
        #Filter with the same peptides
        test1 <- filter(df_final,  !is.na(df_final$`Peptides(95%).x` & df_final$`Peptides(95%).y`))
        test1 <- with(test1,  test1[order(-test1$`Peptides(95%).x`) , ])
        #Different peptides between 1º and the 2º. 
        test2 <- filter(df_final,  !is.na(df_final$`Peptides(95%).x`) & is.na(df_final$`Peptides(95%).y`))
        test2 <- with(test2,  test2[order(-test2$`Peptides(95%).x`) , ])
        #Different peptides between 2º and the 1º.
        test3 <- filter(df_final,  is.na(df_final$`Peptides(95%).x`) & !is.na(df_final$`Peptides(95%).y`))
        test3 <- with(test3,  test3[order(-test3$`Peptides(95%).y`) , ])
        
        
        #Merge data. 

        Comunes2a2 <- test1
        Comunes1a1 <- rbind(test2,test3)
        
        mylist_comunes <- mget(ls(pattern = "Comunes*"))
        dfList <- list()
        
        for(i in 1:length(mylist_comunes)){
          
          
          #Añadimos columna N 
          test_final <- mylist_comunes[[i]]

        #Añadimos columna N 
          test_final2 <- data.frame(cbind(N = 1:nrow(test_final), test_final))
          test_final2$N.x <- NULL
          test_final2$N.y <- NULL 
        
        #Eliminales la columna X__1
          if("X__1"%in% colnames(test_final2)){
            test_final2 <- subset(test_final2, select = -c(X__1) )
            colnames(test_final2) <- c("N","Name","Accesion","Score","%Cov(95)","Peptides(95%)","Species","Score","%Cov(95)","Peptides(95%)","Species")
            
          }
        
        #Cambiamos nombres
          colnames(test_final2) <- c("N","Name","Accesion","Score","%Cov(95)","Peptides(95%)","Species","Score","%Cov(95)","Peptides(95%)","Species")
        
        #Round: 
          round_df <- function(df, digits) {
            nums <- vapply(df, is.numeric, FUN.VALUE = logical(1))
            
            df[,nums] <- round(df[,nums], digits = digits)
            
            (df)
          }
        #Round
          test_final3 <- round_df(test_final2, digits=2)
          
          
          #Save in a list
          dfList[[1]] = data.frame(test_final3)
          
          names(dfList)<-sprintf(paste(gsub(".*","",mylist_comunes[1:length(dfList)]), "Sumary_vs_Control",1:length(dfList), sep = "", na=""),1:length(dfList))
          
          
          #Export data frame to table.
          WriteXLS(dfList, ExcelFileName = "Summaries_vs_Control.xlsx", names(dfList))
        

        #Export data frame to table.
        WriteXLS(test_final3, ExcelFileName = "Multiconsenso_2.xls", SheetNames = NULL, perl = "perl",
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
              restoreConsole=TRUE,"Venn_Diagram_2.tiff")
        
        
        v <- venn.diagram(list(Muestra1=df1$Accession,Muestra2=df2$Accession),
                          fill = c("red", "blue"),
                          cat.cex = 1.5, cex=1.5,cat.pos=0,
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
        
        #Clasification. Se filtran los resultados
        #Filter with the same peptides
        test1 <- filter(df_final,  !is.na(df_final$`Peptides(95%).x`) & !is.na(df_final$`Peptides(95%).y`)& !is.na(df_final$`Peptides(95%)`))
        test1 <- with(test1,  test1[order(-test1$`Peptides(95%).x`) , ])
        #Different peptides between 1º and the 2º. 
        test2 <- filter(df_final,  !is.na(df_final$`Peptides(95%).x`) & !is.na(df_final$`Peptides(95%).y`)& is.na(df_final$`Peptides(95%)`))
        test2 <- with(test2,  test2[order(-test2$`Peptides(95%).x`) , ])
        #Different peptides between 2º and the 1º.
        test3 <- filter(df_final,  is.na(df_final$`Peptides(95%).x`) & !is.na(df_final$`Peptides(95%).y`)& !is.na(df_final$`Peptides(95%)`))
        test3 <- with(test3,  test3[order(-test3$`Peptides(95%).y`) , ])
        #Different between 3º and 2º,1º
        test4 <- filter(df_final,  !is.na(df_final$`Peptides(95%).x`) & is.na(df_final$`Peptides(95%).y`)& !is.na(df_final$`Peptides(95%)`))
        test4 <- with(test4,  test4[order(-test4$`Peptides(95%)`) , ])
        
        test5 <- filter(df_final,  !is.na(df_final$`Peptides(95%).x`) & is.na(df_final$`Peptides(95%).y`)& is.na(df_final$`Peptides(95%)`))
        test5 <- with(test5,  test5[order(-test5$`Peptides(95%).x`) , ])
        
        test6 <- filter(df_final,  is.na(df_final$`Peptides(95%).x`) & !is.na(df_final$`Peptides(95%).y`)& is.na(df_final$`Peptides(95%)`))
        test6 <- with(test6,  test6[order(-test6$`Peptides(95%).y`) , ])
        
        test7 <- filter(df_final,  is.na(df_final$`Peptides(95%).x`) &  is.na(df_final$`Peptides(95%).y`)& !is.na(df_final$`Peptides(95%)`))
        test7 <- with(test7,  test7[order(-test7$`Peptides(95%)`) , ])
      
        
        #Merge data. 
        test_final<-rbind(test1, test2, test3, test4, test5, test6, test7)
        

        #Añadimos columna N 
        test_final2 <- data.frame(cbind(N = 1:nrow(test_final), test_final))
        test_final2$N.x <- NULL
        test_final2$N.y <- NULL 
        test_final2$N.1 <- NULL
        
        #Eliminales la columna X__1
        if("X__1"%in% colnames(test_final2)){
          test_final2 <- subset(test_final2, select = -c(X__1) )
          colnames(test_final2) <- c("N","Name","Accesion","Score","%Cov(95)","Peptides(95%)","Species","Score","%Cov(95)","Peptides(95%)","Species")
          
        }
        
        #Cambiamos nombres
        colnames(test_final2) <- c("N","Name","Accesion","Score","%Cov(95)","Peptides(95%)","Species","Score","%Cov(95)","Peptides(95%)","Species","Score","%Cov(95)","Peptides(95%)","Species")
        
        #Round: 
        round_df <- function(df, digits) {
          nums <- vapply(df, is.numeric, FUN.VALUE = logical(1))
          
          df[,nums] <- round(df[,nums], digits = digits)
          
          (df)
        }
        #Round
        test_final3 <- round_df(test_final2, digits=2)
        
        
        
        #Export data frame to table.
        WriteXLS(test_final3, ExcelFileName = "Multiconsenso_3.xls", SheetNames = NULL, perl = "perl", 
                 verbose = FALSE, Encoding = c("UTF-8", "latin1", "cp1252"),
                 row.names = FALSE, col.names = TRUE,
                 AdjWidth = FALSE, AutoFilter = FALSE, BoldHeaderRow = FALSE,
                 na = "",
                 FreezeRow = 0, FreezeCol = 0,
                 envir = parent.frame())
        
        
        #Function ven.diagram and grid.
        tiff( width=5, height=5, units="in",
              pointsize=8, compression="lzw", bg="white", res=600,
              restoreConsole=TRUE,"Venn_Diagram_3.tiff")
        
        futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
        
        v <- venn.diagram(list(Muestra1=df1$Accession,Muestra2=df2$Accession, Muestra3=df3$Accession),
                          fill = c("red", "blue", "green"),
                          cat.cex = 1.5, cex=1.5,cat.pos=0,
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
        df_final <- Reduce(function(x, y) merge (x, y, by = c("Name", "Accession"), all = TRUE), list(df1, df2, df3))
        
        #Clasification. Se filtran los resultados
        #Filter with the same peptides
        test1 <- filter(df_final,  !is.na(df_final$`Peptides(95%).x`) & !is.na(df_final$`Peptides(95%).y`)& !is.na(df_final$`Peptides(95%)`))
        test1 <- with(test1,  test1[order(-test1$`Peptides(95%).x`) , ])
        #Different peptides between 1º and the 2º. 
        test2 <- filter(df_final,  !is.na(df_final$`Peptides(95%).x`) & !is.na(df_final$`Peptides(95%).y`)& is.na(df_final$`Peptides(95%)`))
        test2 <- with(test2,  test2[order(-test2$`Peptides(95%).x`) , ])
        #Different peptides between 2º and the 1º.
        test3 <- filter(df_final,  is.na(df_final$`Peptides(95%).x`) & !is.na(df_final$`Peptides(95%).y`)& !is.na(df_final$`Peptides(95%)`))
        test3 <- with(test3,  test3[order(-test3$`Peptides(95%).y`) , ])
        #Different between 3º and 2º,1º
        test4 <- filter(df_final,  !is.na(df_final$`Peptides(95%).x`) & is.na(df_final$`Peptides(95%).y`)& !is.na(df_final$`Peptides(95%)`))
        test4 <- with(test4,  test4[order(-test4$`Peptides(95%)`) , ])
        
        test5 <- filter(df_final,  !is.na(df_final$`Peptides(95%).x`) & is.na(df_final$`Peptides(95%).y`)& is.na(df_final$`Peptides(95%)`))
        test5 <- with(test5,  test5[order(-test5$`Peptides(95%).x`) , ])
        
        test6 <- filter(df_final,  is.na(df_final$`Peptides(95%).x`) & !is.na(df_final$`Peptides(95%).y`)& is.na(df_final$`Peptides(95%)`))
        test6 <- with(test6,  test6[order(-test6$`Peptides(95%).y`) , ])
        
        test7 <- filter(df_final,  is.na(df_final$`Peptides(95%).x`) &  is.na(df_final$`Peptides(95%).y`)& !is.na(df_final$`Peptides(95%)`))
        test7 <- with(test7,  test7[order(-test7$`Peptides(95%)`) , ])
        
        test8 <- filter(df_final,  is.na(df_final$`Peptides(95%).x`) &  is.na(df_final$`Peptides(95%).y`)& !is.na(df_final$`Peptides(95%)`))
        test8 <- with(test8,  test8[order(-test8$`Peptides(95%)`) , ])
        
        test9 <- filter(df_final,  is.na(df_final$`Peptides(95%).x`) &  is.na(df_final$`Peptides(95%).y`)& !is.na(df_final$`Peptides(95%)`))
        test9 <- with(test9,  test9[order(-test9$`Peptides(95%)`) , ])
        
        test10 <- filter(df_final,  is.na(df_final$`Peptides(95%).x`) &  is.na(df_final$`Peptides(95%).y`)& !is.na(df_final$`Peptides(95%)`))
        test10 <- with(test10,  test10[order(-test10$`Peptides(95%)`) , ])
        
        test11 <- filter(df_final,  is.na(df_final$`Peptides(95%).x`) &  is.na(df_final$`Peptides(95%).y`)& !is.na(df_final$`Peptides(95%)`))
        test11 <- with(test11,  test11[order(-test11$`Peptides(95%)`) , ])
        
        test12 <- filter(df_final,  is.na(df_final$`Peptides(95%).x`) &  is.na(df_final$`Peptides(95%).y`)& !is.na(df_final$`Peptides(95%)`))
        test12 <- with(test12,  test12[order(-test12$`Peptides(95%)`) , ])
        
        test13 <- filter(df_final,  is.na(df_final$`Peptides(95%).x`) &  is.na(df_final$`Peptides(95%).y`)& !is.na(df_final$`Peptides(95%)`))
        test13 <- with(test13,  test13[order(-test13$`Peptides(95%)`) , ])
        
        test14 <- filter(df_final,  is.na(df_final$`Peptides(95%).x`) &  is.na(df_final$`Peptides(95%).y`)& !is.na(df_final$`Peptides(95%)`))
        test14 <- with(test14,  test14[order(-test14$`Peptides(95%)`) , ])
        
        test15 <- filter(df_final,  is.na(df_final$`Peptides(95%).x`) &  is.na(df_final$`Peptides(95%).y`)& !is.na(df_final$`Peptides(95%)`))
        test15 <- with(test15,  test15[order(-test15$`Peptides(95%)`) , ])
        
        
        #Merge data. 
        test_final<-rbind(test1, test2, test3, test4, test5, test6, test7)
        
        
        #Añadimos columna N 
        test_final2 <- data.frame(cbind(N = 1:nrow(test_final), test_final))
        test_final2$N.x <- NULL
        test_final2$N.y <- NULL 
        test_final2$N.1 <- NULL
        
        #Eliminales la columna X__1
        if("X__1"%in% colnames(test_final2)){
          test_final2 <- subset(test_final2, select = -c(X__1) )
          colnames(test_final2) <- c("N","Name","Accesion","Score","%Cov(95)","Peptides(95%)","Species","Score","%Cov(95)","Peptides(95%)","Species")
          
        }
        
        #Cambiamos nombres
        colnames(test_final2) <- c("N","Name","Accesion","Score","%Cov(95)","Peptides(95%)","Species","Score","%Cov(95)","Peptides(95%)","Species","Score","%Cov(95)","Peptides(95%)","Species")
        
        #Round: 
        round_df <- function(df, digits) {
          nums <- vapply(df, is.numeric, FUN.VALUE = logical(1))
          
          df[,nums] <- round(df[,nums], digits = digits)
          
          (df)
        }
        #Round
        test_final3 <- round_df(test_final2, digits=2)
        
        
        
        #Export data frame to table.
        WriteXLS(test_final3, ExcelFileName = "Multiconsenso_4.xls", SheetNames = NULL, perl = "perl", 
                 verbose = FALSE, Encoding = c("UTF-8", "latin1", "cp1252"),
                 row.names = FALSE, col.names = TRUE,
                 AdjWidth = FALSE, AutoFilter = FALSE, BoldHeaderRow = FALSE,
                 na = "",
                 FreezeRow = 0, FreezeCol = 0,
                 envir = parent.frame())
        
        
        #Function ven.diagram and grid.
        tiff( width=5, height=5, units="in",
              pointsize=8, compression="lzw", bg="white", res=600,
              restoreConsole=TRUE,"Venn_Diagram_4.tiff")
        
        futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
        
        v <- venn.diagram(list(Muestra1=df1$Accession,Muestra2=df2$Accession, Muestra3=df3$Accession, Muestra4=df4$Accession),
                          fill = c("red", "blue", "green", "yellow"),
                          cat.cex = 1.5, cex=1.5,cat.pos=0,
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
        df5 <- read_excel(file.path(getwd(), files_glob_no_control[4]))
        

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
}



