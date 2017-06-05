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



repeat{
  cat("¿Quieres hacer la comparativa con un control?\n ")
  cat("1. Yes | 2. No\n")
  x <- readLines(file("stdin"),1)#enter "yes"

  if (x=="yes" | x==1 | x=="Yes")  {

    #Check how many controls there are. 
    if (length(control)<1){
      print("No hay un archivo con t�tulo Control")
      break 
    }

    #Bucle if primero:
    dfList <- list()

    for(i in 1:length(files_glob_no_control)){
      df1 <- read_excel(file.path(getwd(), grep('Control', files_glob, value=TRUE)))
      df2 <- read_excel(file.path(getwd(), files_glob_no_control[i]))

      if("X__1"%in% colnames(df1)){
        df1 <- subset(df1, select = -c(X__1) )
      } 
      if("X__1"%in% colnames(df2)){
        df2 <- subset(df2, select = -c(X__1) )
      }

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

      names(dfList)<-sprintf(paste(gsub(".*","",files_glob_no_control[1:length(dfList)]), "Control_vs_Summary",1:length(dfList), sep = "", na=""),1:length(dfList))

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

    #Export data frame to table.
    WriteXLS(dfList, ExcelFileName = "Summaries_vs_Control.xlsx", names(dfList))

    break

  }else if(x=="no" | x==2 | x=="No"){

    if (length(files_glob_no_control) == 1){ 
      print("There is only 1 xlsx file")

      break

    }else if (length(files_glob_no_control) == 2){

      #Load excel file
      df1 <- read_excel(file.path(getwd(), files_glob_no_control[1]))
      df2 <- read_excel(file.path(getwd(), files_glob_no_control[2]))

      #Remove X__1
      if("X__1"%in% colnames(df1)){
        df1 <- subset(df1, select = -c(X__1) )
      } 
      if("X__1"%in% colnames(df2)){
        df2 <- subset(df2, select = -c(X__1) )
      }


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
      Comunes_ALL<- rbind(test1, test2,test3)
      Comunes2a2 <- test1
      Comunes1a1 <- rbind(test2,test3)

      #Check variables that have comunes
      mylist_comunes <- mget(ls(pattern = "Comunes*"))
      #Create a empty list
      dfList <- list()

      for(i in 1:length(mylist_comunes)){
        #Load dataframe in test_final
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
        colnames(test_final3) <- c("N","Name","Accesion","Score","%Cov(95)","Peptides(95%)","Species","Score","%Cov(95)","Peptides(95%)","Species")

        #Save in a list
        dfList[[i]] = data.frame(test_final3)

        names(dfList)<-sprintf(paste(gsub(".*","",mylist_comunes[1:length(dfList)]), "",names(mylist_comunes)[1:length(dfList)], sep = "", na=""),1:length(dfList))

      }

      #Export data frame to table.

      WriteXLS(dfList, ExcelFileName = "Multiconsenso_2.xls", names(dfList))

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

      break


    } else if (length(files_glob_no_control) == 3){
      #Load excel file
      df1 <- read_excel(file.path(getwd(), files_glob_no_control[1]))
      df2 <- read_excel(file.path(getwd(), files_glob_no_control[2]))
      df3 <- read_excel(file.path(getwd(), files_glob_no_control[3]))

      #Remove X__1
      if("X__1"%in% colnames(df1)){
        df1 <- subset(df1, select = -c(X__1) )
      } 
      if("X__1"%in% colnames(df2)){
        df2 <- subset(df2, select = -c(X__1) )
      } 
      if("X__1"%in% colnames(df3)){
        df3 <- subset(df3, select = -c(X__1) )
      }


      #Combine data frame using reduce function
      df_final <- Reduce(function(x, y) merge (x, y, by = c("Name", "Accession"), all = TRUE), list(df1, df2, df3))
      #colnames(df_final) <- c("Name","Accesion","N","Score","%Cov(95)","Peptides(95%)","Species","N.y","Score.y","%Cov(95).y","Peptides(95%).y","Species.y","N.z","Score.z","%Cov(95).z","Peptides(95%).z","Species.z")


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


      #Lista de comunes
      Comunes_ALL <- rbind(test1,test2,test3,test4, test5, test6,test7)
      Comunes1a1 <- rbind(test5, test6,test7)
      Comunes2a2 <- rbind(test2,test3,test4)
      Comunes3a3 <- test1

      #Check variables that have comunes
      mylist_comunes <- mget(ls(pattern = "Comunes*"))
      #Create a empty list
      dfList <- list()

      for(i in 1:length(mylist_comunes)){
        #Load dataframe in test_final
        test_final <- mylist_comunes[[i]]

        #Añadimos columna N 
        test_final2 <- data.frame(cbind(N = 1:nrow(test_final), test_final))
        test_final2$N.x <- NULL
        test_final2$N.y <- NULL 

        #Cambiamos nombres
        #colnames(test_final2) <- c("N","Name","Accesion","Score","%Cov(95)","Peptides(95%)","Species","Score","%Cov(95)","Peptides(95%)","Species","Score","%Cov(95)","Peptides(95%)","Species")

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

        names(dfList)<-sprintf(paste(gsub(".*","",mylist_comunes[1:length(dfList)]), "",names(mylist_comunes)[1:length(dfList)], sep = "", na=""),1:length(dfList))
      }

      #Export data frame to table.

      WriteXLS(dfList, ExcelFileName = "Multiconsenso_3.xls", names(dfList))

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

      break


    } else if (length(files_glob_no_control) == 4){
      df1 <- read_excel(file.path(getwd(), files_glob_no_control[1]))
      df2 <- read_excel(file.path(getwd(), files_glob_no_control[2]))
      df3 <- read_excel(file.path(getwd(), files_glob_no_control[3]))
      df4 <- read_excel(file.path(getwd(), files_glob_no_control[4]))

      if("X__1"%in% colnames(df1)){
        df1 <- subset(df1, select = -c(X__1) )
      } 
      if("X__1"%in% colnames(df2)){
        df2 <- subset(df2, select = -c(X__1) )
      } 
      if("X__1"%in% colnames(df3)){
        df3 <- subset(df3, select = -c(X__1) )
      }
      if("X__1"%in% colnames(df4)){
        df4 <- subset(df4, select = -c(X__1) )
      }


      #Combine data frame using reduce function
      df_final <- Reduce(function(x, y) merge (x, y, by = c("Name", "Accession"), all = TRUE), list(df1, df2, df3, df4))

      colnames(df_final) <- c("Name","Accesion","N.x","Score.x","%Cov(95).x","Peptides(95%).x","Species.x", "N.y","Score.y","%Cov(95).y","Peptides(95%).y","Species.y","N.z","Score.z","%Cov(95).z","Peptides(95%).z","Species.z","N.w","Score.w","%Cov(95).w","Peptides(95%).w","Species.w")


      #Clasification. Se filtran los resultados
      #ALL
      test1 <- filter(df_final, !is.na(df_final$`Peptides(95%).x`)& !is.na(df_final$`Peptides(95%).y`)& !is.na(df_final$`Peptides(95%).z`)& !is.na(df_final$`Peptides(95%).w`))
      test1 <- with(test1,  test1[order(-test1$`Peptides(95%).x`) , ])
      #Different peptides between 1º and the 2º. 
      test2 <- filter(df_final, !is.na(df_final$`Peptides(95%).x`)& !is.na(df_final$`Peptides(95%).y`)& !is.na(df_final$`Peptides(95%).z`)& is.na(df_final$`Peptides(95%).w`))
      test2 <- with(test2,  test2[order(-test2$`Peptides(95%).x`) , ])
      #Different peptides between 2º and the 1º.
      test3 <- filter(df_final, !is.na(df_final$`Peptides(95%).x`)& !is.na(df_final$`Peptides(95%).y`)& is.na(df_final$`Peptides(95%).z`)& !is.na(df_final$`Peptides(95%).w`))
      test3 <- with(test3,  test3[order(-test3$`Peptides(95%).x`) , ])
      #Different between 3º and 2º,1º
      test4 <- filter(df_final, !is.na(df_final$`Peptides(95%).x`)& is.na(df_final$`Peptides(95%).y`)& !is.na(df_final$`Peptides(95%).z`)& !is.na(df_final$`Peptides(95%).w`))
      test4 <- with(test4,  test4[order(-test4$`Peptides(95%).y`) , ])

      test5 <- filter(df_final, is.na(df_final$`Peptides(95%).x`)& !is.na(df_final$`Peptides(95%).y`)& !is.na(df_final$`Peptides(95%).z`)& !is.na(df_final$`Peptides(95%).w`))        
      test5 <- with(test5,  test5[order(-test5$`Peptides(95%).y`) , ])

      test6 <- filter(df_final, !is.na(df_final$`Peptides(95%).x`)& !is.na(df_final$`Peptides(95%).y`)& is.na(df_final$`Peptides(95%).z`)& is.na(df_final$`Peptides(95%).w`))        
      test6 <- with(test6,  test6[order(-test6$`Peptides(95%).x`) , ])

      test7 <- filter(df_final, !is.na(df_final$`Peptides(95%).x`)& is.na(df_final$`Peptides(95%).y`)& !is.na(df_final$`Peptides(95%).z`)& is.na(df_final$`Peptides(95%).w`))
      test7 <- with(test7,  test7[order(-test7$`Peptides(95%).x`) , ])

      test8 <- filter(df_final, !is.na(df_final$`Peptides(95%).x`)& is.na(df_final$`Peptides(95%).y`)& is.na(df_final$`Peptides(95%).z`)& !is.na(df_final$`Peptides(95%).w`))
      test8 <- with(test8,  test8[order(-test8$`Peptides(95%).x`) , ])

      test9 <- filter(df_final, is.na(df_final$`Peptides(95%).x`)& !is.na(df_final$`Peptides(95%).y`)& !is.na(df_final$`Peptides(95%).z`)& is.na(df_final$`Peptides(95%).w`))
      test9 <- with(test9,  test9[order(-test9$`Peptides(95%).y`) , ])

      test10 <- filter(df_final, is.na(df_final$`Peptides(95%).x`)& !is.na(df_final$`Peptides(95%).y`)& is.na(df_final$`Peptides(95%).z`)& !is.na(df_final$`Peptides(95%).w`))
      test10 <- with(test10,  test10[order(-test10$`Peptides(95%).y`) , ])

      test11 <- filter(df_final, is.na(df_final$`Peptides(95%).x`)& is.na(df_final$`Peptides(95%).y`)& !is.na(df_final$`Peptides(95%).z`)& !is.na(df_final$`Peptides(95%).w`))
      test11 <- with(test11,  test11[order(-test11$`Peptides(95%).z`) , ])

      test12 <- filter(df_final, !is.na(df_final$`Peptides(95%).x`)& is.na(df_final$`Peptides(95%).y`)& is.na(df_final$`Peptides(95%).z`)& is.na(df_final$`Peptides(95%).w`))
      test12 <- with(test12,  test12[order(-test12$`Peptides(95%).x`) , ])

      test13 <- filter(df_final, is.na(df_final$`Peptides(95%).x`)& !is.na(df_final$`Peptides(95%).y`)& is.na(df_final$`Peptides(95%).z`)& is.na(df_final$`Peptides(95%).w`))
      test13 <- with(test13,  test13[order(-test13$`Peptides(95%).y`) , ])

      test14 <- filter(df_final, is.na(df_final$`Peptides(95%).x`)& is.na(df_final$`Peptides(95%).y`)& !is.na(df_final$`Peptides(95%).z`)& is.na(df_final$`Peptides(95%).w`))
      test14 <- with(test14,  test14[order(-test14$`Peptides(95%).z`) , ])

      test15 <- filter(df_final, is.na(df_final$`Peptides(95%).x`)& is.na(df_final$`Peptides(95%).y`)& is.na(df_final$`Peptides(95%).z`)& !is.na(df_final$`Peptides(95%).w`))
      test15 <- with(test15,  test15[order(-test15$`Peptides(95%).w`) , ])


      #Merge data. 
      Comunes_ALL<-rbind(test1, test2, test3, test4, test5, test6, test7, test8, test9, test10, test11, test12, test13, test14, test15)

      Comunes4a4 <- test1
      Comunes3a3 <- rbind(test2,test3,test4,test5)
      Comunes2a2 <- rbind(test6, test7, test8, test9, test10, test11)
      Comunes1a1 <- rbind(test12, test13, test14, test15)


      #Check variables that have comunes
      mylist_comunes <- mget(ls(pattern = "Comunes*"))
      #Create a empty list
      dfList <- list()

      for(i in 1:length(mylist_comunes)){
        #Load dataframe in test_final
        test_final <- mylist_comunes[[i]]

        #Añadimos columna N 
        test_final2 <- data.frame(cbind(N = 1:nrow(test_final), test_final))
        test_final2$N.x <- NULL
        test_final2$N.y <- NULL 
        test_final2$N.z <- NULL
        test_final2$N.w <- NULL

        #Cambiamos nombres
        #colnames(test_final2) <- c("N","Name","Accesion","Score","%Cov(95)","Peptides(95%)","Species","Score","%Cov(95)","Peptides(95%)","Species","Score","%Cov(95)","Peptides(95%)","Species")

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

        names(dfList)<-sprintf(paste(gsub(".*","",mylist_comunes[1:length(dfList)]), "",names(mylist_comunes)[1:length(dfList)], sep = "", na=""),1:length(dfList))

      }

      #Export data frame to table.

      WriteXLS(dfList, ExcelFileName = "Multiconsenso_4.xls", names(dfList))

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

      break

    }else if (length(files_glob_no_control) == 5){
      df1 <- read_excel(file.path(getwd(), files_glob_no_control[1]))
      df2 <- read_excel(file.path(getwd(), files_glob_no_control[2]))
      df3 <- read_excel(file.path(getwd(), files_glob_no_control[3]))
      df4 <- read_excel(file.path(getwd(), files_glob_no_control[4]))
      df5 <- read_excel(file.path(getwd(), files_glob_no_control[5]))

      if("X__1"%in% colnames(df1)){
        df1 <- subset(df1, select = -c(X__1) )
      } 
      if("X__1"%in% colnames(df2)){
        df2 <- subset(df2, select = -c(X__1) )
      } 
      if("X__1"%in% colnames(df3)){
        df3 <- subset(df3, select = -c(X__1) )
      }
      if("X__1"%in% colnames(df4)){
        df4 <- subset(df4, select = -c(X__1) )
      }
      if("X__1"%in% colnames(df5)){
        df5 <- subset(df5, select = -c(X__1) )
      }


      #Combine data frame using reduce function
      df_final <- Reduce(function(x, y) merge (x, y, by = c("Name", "Accession"), all = TRUE), list(df1, df2, df3, df4,df5))

      colnames(df_final) <- c("Name","Accesion","N.x","Score.x","%Cov(95).x","Peptides(95%).x","Species.x", "N.y","Score.y","%Cov(95).y","Peptides(95%).y","Species.y","N.z","Score.z","%Cov(95).z","Peptides(95%).z","Species.z","N.w","Score.w","%Cov(95).w","Peptides(95%).w","Species.w","N.u","Score.u","%Cov(95).u","Peptides(95%).u","Species.u")


      #Clasification. Se filtran los resultados
      #ALL
      test1 <- filter(df_final, !is.na(df_final$`Peptides(95%).x`)& !is.na(df_final$`Peptides(95%).y`)& !is.na(df_final$`Peptides(95%).z`)& !is.na(df_final$`Peptides(95%).w`)& !is.na(df_final$`Peptides(95%).u`))
      test1 <- with(test1,  test1[order(-test1$`Peptides(95%).x`) , ])
      #Different peptides between 1º and the 2º. 
      test2 <- filter(df_final, !is.na(df_final$`Peptides(95%).x`)& !is.na(df_final$`Peptides(95%).y`)& !is.na(df_final$`Peptides(95%).z`)& !is.na(df_final$`Peptides(95%).w`)& is.na(df_final$`Peptides(95%).u`))
      test2 <- with(test2,  test2[order(-test2$`Peptides(95%).x`) , ])
      #Different peptides between 2º and the 1º.
      test3 <- filter(df_final, !is.na(df_final$`Peptides(95%).x`)& !is.na(df_final$`Peptides(95%).y`)& !is.na(df_final$`Peptides(95%).z`)& is.na(df_final$`Peptides(95%).w`)& !is.na(df_final$`Peptides(95%).u`))
      test3 <- with(test3,  test3[order(-test3$`Peptides(95%).x`) , ])
      #Different between 3º and 2º,1º
      test4 <- filter(df_final, !is.na(df_final$`Peptides(95%).x`)& !is.na(df_final$`Peptides(95%).y`)& is.na(df_final$`Peptides(95%).z`)& !is.na(df_final$`Peptides(95%).w`)& !is.na(df_final$`Peptides(95%).u`))
      test4 <- with(test4,  test4[order(-test4$`Peptides(95%).y`) , ])

      test5 <- filter(df_final, !is.na(df_final$`Peptides(95%).x`)& is.na(df_final$`Peptides(95%).y`)& !is.na(df_final$`Peptides(95%).z`)& !is.na(df_final$`Peptides(95%).w`)& !is.na(df_final$`Peptides(95%).u`))        
      test5 <- with(test5,  test5[order(-test5$`Peptides(95%).y`) , ])

      test6 <- filter(df_final, is.na(df_final$`Peptides(95%).x`)& !is.na(df_final$`Peptides(95%).y`)& !is.na(df_final$`Peptides(95%).z`)& !is.na(df_final$`Peptides(95%).w`)& !is.na(df_final$`Peptides(95%).u`))        
      test6 <- with(test6,  test6[order(-test6$`Peptides(95%).x`) , ])

      test7 <- filter(df_final, !is.na(df_final$`Peptides(95%).x`)& !is.na(df_final$`Peptides(95%).y`)& !is.na(df_final$`Peptides(95%).z`)& is.na(df_final$`Peptides(95%).w`)& is.na(df_final$`Peptides(95%).u`))
      test7 <- with(test7,  test7[order(-test7$`Peptides(95%).x`) , ])

      test8 <- filter(df_final, !is.na(df_final$`Peptides(95%).x`)& !is.na(df_final$`Peptides(95%).y`)& is.na(df_final$`Peptides(95%).z`)& !is.na(df_final$`Peptides(95%).w`)& is.na(df_final$`Peptides(95%).u`))
      test8 <- with(test8,  test8[order(-test8$`Peptides(95%).x`) , ])

      test9 <- filter(df_final, !is.na(df_final$`Peptides(95%).x`)& !is.na(df_final$`Peptides(95%).y`)& is.na(df_final$`Peptides(95%).z`)& is.na(df_final$`Peptides(95%).w`)& !is.na(df_final$`Peptides(95%).u`))
      test9 <- with(test9,  test9[order(-test9$`Peptides(95%).y`) , ])

      test10 <- filter(df_final, !is.na(df_final$`Peptides(95%).x`)& is.na(df_final$`Peptides(95%).y`)& !is.na(df_final$`Peptides(95%).z`)& !is.na(df_final$`Peptides(95%).w`)& is.na(df_final$`Peptides(95%).u`))
      test10 <- with(test10,  test10[order(-test10$`Peptides(95%).y`) , ])

      test11 <- filter(df_final, !is.na(df_final$`Peptides(95%).x`)& is.na(df_final$`Peptides(95%).y`)& !is.na(df_final$`Peptides(95%).z`)& is.na(df_final$`Peptides(95%).w`)& !is.na(df_final$`Peptides(95%).u`))
      test11 <- with(test11,  test11[order(-test11$`Peptides(95%).z`) , ])

      test12 <- filter(df_final, !is.na(df_final$`Peptides(95%).x`)& is.na(df_final$`Peptides(95%).y`)& is.na(df_final$`Peptides(95%).z`)& !is.na(df_final$`Peptides(95%).w`)& !is.na(df_final$`Peptides(95%).u`))
      test12 <- with(test12,  test12[order(-test12$`Peptides(95%).x`) , ])

      test13 <- filter(df_final, is.na(df_final$`Peptides(95%).x`)& !is.na(df_final$`Peptides(95%).y`)& !is.na(df_final$`Peptides(95%).z`)& !is.na(df_final$`Peptides(95%).w`)& is.na(df_final$`Peptides(95%).u`))
      test13 <- with(test13,  test13[order(-test13$`Peptides(95%).y`) , ])

      test14 <- filter(df_final, is.na(df_final$`Peptides(95%).x`)& !is.na(df_final$`Peptides(95%).y`)& is.na(df_final$`Peptides(95%).z`)& !is.na(df_final$`Peptides(95%).w`)& !is.na(df_final$`Peptides(95%).u`))
      test14 <- with(test14,  test14[order(-test14$`Peptides(95%).z`) , ])

      test15 <- filter(df_final, is.na(df_final$`Peptides(95%).x`)& !is.na(df_final$`Peptides(95%).y`)& !is.na(df_final$`Peptides(95%).z`)& is.na(df_final$`Peptides(95%).w`)& !is.na(df_final$`Peptides(95%).u`))
      test15 <- with(test15,  test15[order(-test15$`Peptides(95%).w`) , ])

      test16 <- filter(df_final, is.na(df_final$`Peptides(95%).x`)& is.na(df_final$`Peptides(95%).y`)& !is.na(df_final$`Peptides(95%).z`)& !is.na(df_final$`Peptides(95%).w`)& !is.na(df_final$`Peptides(95%).u`))
      test16 <- with(test16,  test16[order(-test16$`Peptides(95%).w`) , ])

      test17 <- filter(df_final, !is.na(df_final$`Peptides(95%).x`)& is.na(df_final$`Peptides(95%).y`)& is.na(df_final$`Peptides(95%).z`)& !is.na(df_final$`Peptides(95%).w`)& is.na(df_final$`Peptides(95%).u`))
      test17 <- with(test17,  test17[order(-test17$`Peptides(95%).w`) , ])

      test18 <- filter(df_final, !is.na(df_final$`Peptides(95%).x`)& is.na(df_final$`Peptides(95%).y`)& is.na(df_final$`Peptides(95%).z`)& is.na(df_final$`Peptides(95%).w`)& !is.na(df_final$`Peptides(95%).u`))
      test18 <- with(test18,  test18[order(-test18$`Peptides(95%).w`) , ])

      test19 <- filter(df_final, !is.na(df_final$`Peptides(95%).x`)& is.na(df_final$`Peptides(95%).y`)& !is.na(df_final$`Peptides(95%).z`)& is.na(df_final$`Peptides(95%).w`)& is.na(df_final$`Peptides(95%).u`))
      test19 <- with(test19,  test19[order(-test19$`Peptides(95%).w`) , ])

      test20 <- filter(df_final, !is.na(df_final$`Peptides(95%).x`)& !is.na(df_final$`Peptides(95%).y`)& is.na(df_final$`Peptides(95%).z`)& is.na(df_final$`Peptides(95%).w`)& is.na(df_final$`Peptides(95%).u`))
      test20 <- with(test20,  test20[order(-test20$`Peptides(95%).w`) , ])

      test21 <- filter(df_final, is.na(df_final$`Peptides(95%).x`)& !is.na(df_final$`Peptides(95%).y`)& is.na(df_final$`Peptides(95%).z`)& !is.na(df_final$`Peptides(95%).w`)& is.na(df_final$`Peptides(95%).u`))
      test21 <- with(test21,  test21[order(-test21$`Peptides(95%).w`) , ])

      test22 <- filter(df_final, is.na(df_final$`Peptides(95%).x`)& !is.na(df_final$`Peptides(95%).y`)& !is.na(df_final$`Peptides(95%).z`)& is.na(df_final$`Peptides(95%).w`)& is.na(df_final$`Peptides(95%).u`))
      test22 <- with(test22,  test22[order(-test22$`Peptides(95%).w`) , ])

      test23 <- filter(df_final, is.na(df_final$`Peptides(95%).x`)& !is.na(df_final$`Peptides(95%).y`)& is.na(df_final$`Peptides(95%).z`)& is.na(df_final$`Peptides(95%).w`)& !is.na(df_final$`Peptides(95%).u`))
      test23 <- with(test23,  test23[order(-test23$`Peptides(95%).w`) , ])

      test24 <- filter(df_final, is.na(df_final$`Peptides(95%).x`)& is.na(df_final$`Peptides(95%).y`)& !is.na(df_final$`Peptides(95%).z`)& is.na(df_final$`Peptides(95%).w`)& !is.na(df_final$`Peptides(95%).u`))
      test24 <- with(test24,  test24[order(-test24$`Peptides(95%).w`) , ])

      test25 <- filter(df_final, is.na(df_final$`Peptides(95%).x`)& is.na(df_final$`Peptides(95%).y`)& !is.na(df_final$`Peptides(95%).z`)& !is.na(df_final$`Peptides(95%).w`)& is.na(df_final$`Peptides(95%).u`))
      test25 <- with(test25,  test25[order(-test25$`Peptides(95%).w`) , ])

      test26 <- filter(df_final, is.na(df_final$`Peptides(95%).x`)& is.na(df_final$`Peptides(95%).y`)& is.na(df_final$`Peptides(95%).z`)& !is.na(df_final$`Peptides(95%).w`)& !is.na(df_final$`Peptides(95%).u`))
      test26 <- with(test26,  test26[order(-test26$`Peptides(95%).w`) , ])

      test27 <- filter(df_final, !is.na(df_final$`Peptides(95%).x`)& is.na(df_final$`Peptides(95%).y`)& is.na(df_final$`Peptides(95%).z`)& is.na(df_final$`Peptides(95%).w`)& is.na(df_final$`Peptides(95%).u`))
      test27 <- with(test27,  test27[order(-test27$`Peptides(95%).w`) , ])

      test28 <- filter(df_final, is.na(df_final$`Peptides(95%).x`)& !is.na(df_final$`Peptides(95%).y`)& is.na(df_final$`Peptides(95%).z`)& is.na(df_final$`Peptides(95%).w`)& is.na(df_final$`Peptides(95%).u`))
      test28 <- with(test28,  test28[order(-test28$`Peptides(95%).w`) , ])

      test29 <- filter(df_final, is.na(df_final$`Peptides(95%).x`)& is.na(df_final$`Peptides(95%).y`)& !is.na(df_final$`Peptides(95%).z`)& is.na(df_final$`Peptides(95%).w`)& is.na(df_final$`Peptides(95%).u`))
      test29 <- with(test29,  test29[order(-test29$`Peptides(95%).w`) , ])

      test30 <- filter(df_final, is.na(df_final$`Peptides(95%).x`)& is.na(df_final$`Peptides(95%).y`)& is.na(df_final$`Peptides(95%).z`)& !is.na(df_final$`Peptides(95%).w`)& is.na(df_final$`Peptides(95%).u`))
      test30 <- with(test30,  test30[order(-test30$`Peptides(95%).w`) , ])

      test31 <- filter(df_final, is.na(df_final$`Peptides(95%).x`)& is.na(df_final$`Peptides(95%).y`)& is.na(df_final$`Peptides(95%).z`)& is.na(df_final$`Peptides(95%).w`)& !is.na(df_final$`Peptides(95%).u`))
      test31 <- with(test31,  test31[order(-test31$`Peptides(95%).w`) , ])

      #Merge data. 
      Comunes_ALL<-rbind(test1, test2, test3, test4, test5, test6, test7, test8, 
      test9, test10, test11, test12, test13, test14, test15,test16,test17,test18,test19,test20,
      test21,test22,test23,test24,test25,test26,test27,test28,test29,test30,test31)

      Comunes5a5 <- test1
      Comunes4a4 <- rbind(test2,test3,test4,test5,test6)
      Comunes3a3 <- rbind(test7, test8, test9, test10, test11, test12, test13, test14, test15, test16,test17)
      Comunes2a2 <- rbind(test18, test19, test20, test21, test22, test23,test24,test25,test26)
      Comunes1a1 <- rbind(test27, test28, test29, test30, test31)


      #Check variables that have comunes
      mylist_comunes <- mget(ls(pattern = "Comunes*"))
      #Create a empty list
      dfList <- list()

      for(i in 1:length(mylist_comunes)){

        #Load dataframe in test_final
        test_final <- mylist_comunes[[i]]

        #Añadimos columna N 
        test_final2 <- data.frame(cbind(N = 1:nrow(test_final), test_final))
        test_final2$N.x <- NULL
        test_final2$N.y <- NULL 
        test_final2$N.z <- NULL
        test_final2$N.w <- NULL



        #Cambiamos nombres
        #colnames(test_final2) <- c("N","Name","Accesion","Score","%Cov(95)","Peptides(95%)","Species","Score","%Cov(95)","Peptides(95%)","Species","Score","%Cov(95)","Peptides(95%)","Species")

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

        names(dfList)<-sprintf(paste(gsub(".*","",mylist_comunes[1:length(dfList)]), "",names(mylist_comunes)[1:length(dfList)], sep = "", na=""),1:length(dfList))

      }

      #Export data frame to table.

      WriteXLS(dfList, ExcelFileName = "Multiconsenso_5.xls", names(dfList))

      #Function ven.diagram and grid.
      tiff( width=10, height=10, units="in",
      pointsize=8, compression="lzw", bg="white", res=600,
      restoreConsole=TRUE,"Venn_Diagram_5.tiff")

      futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

      v <- venn.diagram(list(Muestra1=df1$Accession,Muestra2=df2$Accession, Muestra3=df3$Accession, Muestra4=df4$Accession, Muestra5=df5$Accession),
      fill = c("red", "blue", "green", "yellow", "orange"),
      cat.cex = 1.5, cex=1.5,cat.pos=0,
      filename=NULL)


      # have a look at the default plot
      grid.newpage()
      grid.draw(v)
      garbage <- dev.off()

      break
    }
  }
}

