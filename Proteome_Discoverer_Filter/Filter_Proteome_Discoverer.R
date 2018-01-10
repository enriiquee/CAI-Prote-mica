######################################################################
#FILTER_Proteome_Discoverer: 
#This script filter the results from Proteome Discoverer 2.2
######################################################################


list.of.packages <- c("readr", "readxl","WriteXLS", "openxlsx")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library('readr'); library('readxl'); library('WriteXLS'); library("openxlsx")

#Suppress warnings globally
options(warn = -1)
print("Ejecutando... Puede tardar un poco. ")

readline(prompt = "Recuerda que los archivos deben estar en la misma carpeta que el programa con nombre finalizado en _PeptideSummary y _ProteinSummary \n
         PULSE ENTER [ENTER]")
#Count files
files_glob_peptides <- (Sys.glob("*_PeptideSummary.txt")) 
files_glob_proteins <- (Sys.glob("*_ProteinSummary.txt")) 


if (length(files_glob_peptides) != length(files_glob_proteins)) {
  print("WARNING!! Different number of files: Protein and Peptides")
} else {
  i <- 1
  while(i <= length(files_glob_peptides)){
    

    #Suppress warnings globally
    options(warn = -1)
    print("Ejecutando... Puede tardar un poco. ")
    
    
    #Loading data
    # Proteins_PP <- read_delim("Y:/ENRIQUE/Filtro PP en R/17-42_Muestra 161_NCBISmMt_ProteinSummary.txt", 
    #                           "\t", escape_double = FALSE, trim_ws = TRUE)
    # 
    # Peptidos_PP <- read_delim("Y:/ENRIQUE/Filtro PP en R/17-42_Muestra 161_NCBISmMt_PeptideSummary.txt", 
    #                           "\t", escape_double = FALSE, trim_ws = TRUE)
    
    
    Proteins_PP <- read_delim(files_glob_proteins[i], 
                              "\t", escape_double = FALSE, trim_ws = TRUE, col_types = cols())
    #Eliminamos los NA de la columna accesion
    Proteins_PP <- Proteins_PP[!is.na(Proteins_PP$Accession), ]
    
    #Proteins_PP <- data.frame/na.omit(Proteins_PP$Accession)
    
    Peptidos_PP <- read_delim(files_glob_peptides[1],
                              "\t", escape_double = FALSE, trim_ws = TRUE, col_types = cols())
    
    #Eliminamos los NA de la columna accesion
    Peptidos_PP <- Peptidos_PP[!is.na(Peptidos_PP$`Master Protein Accessions`), ]
    

    #choose the columns that we want. 
    Proteins_PP2 <- subset( Proteins_PP, select = -c(1,2))

    Peptidos_PP2 <- subset( Peptidos_PP, select = -c(1,2,8))
    
    #Add Master Protein Accesions 
    Peptidos_PP2 <- data.frame("Accesion"=gsub("\\;.*","", Peptidos_PP$`Master Protein Accessions`), Peptidos_PP2, stringsAsFactors = FALSE)
    
    ##############
    # #comparar tablas test: 
    # 
    duplicate <- data.frame("Accesion"=gsub("\\;.*","", Peptidos_PP2$Accesion), stringsAsFactors = FALSE)
    duplicate2 <- aggregate(list(numdup=rep(1,nrow(duplicate))), duplicate, length)
    

    duplicate3 <- with(duplicate2,  duplicate2[order(duplicate2$Accesion) , ])
    # 
    # #Protein accesion sorted
    duplicate4 <- with(Proteins_PP2,  Proteins_PP2[order(Proteins_PP2$Accession) , ])
    # #Replace in protein table
    
    #If we need to remove specific protein from Proteins file you can use this. 
    #duplicate4 <- duplicate4[ ! ( ( duplicate4$N ==192)) , ] 
    
    #In order to remove this duplicates, we use if lenght to compare. Then we create 
    #a new table with the reverse elements in peptides, remove letters after ; and duplicates
    # 
    
    # if (length(duplicate4$`Peptides(95%)`)!= length(duplicate3$numdup)){
    #   tabla_reversed <- data.frame(Peptidos_PP2[grepl("REVERSED", Peptidos_PP2$Peptidos_PP.Names),])
    #   tabla_reversed$Peptidos_PP.Accessions <- gsub("\\;.*","", tabla_reversed$Peptidos_PP.Accessions)
    #   tabla_reversed <- unique(tabla_reversed$Peptidos_PP.Accessions)
    #   
    #   #Remove element from the list
    #   duplicate4 <- duplicate4[ ! duplicate4$Accession %in% tabla_reversed, ]
    #   
    # }
    
    #If we want to have the same number of proteins and peptides, we've to remove those who are in proteins but 
    #not in peptides, for that reason we are going to create two new tables. 
    
    #We take the data that we want and remove everthing that is after the first ";".  
    df1 <- data.frame("Accesion"=gsub("\\;.*","", duplicate4$Accession), stringsAsFactors = FALSE)
    df2 <- data.frame("Accesion2"=gsub("\\;.*","", duplicate3$Accesion), stringsAsFactors = FALSE)
    

    #Create a list
    lst <- list(df1, df2)
    
    #Create a list where we can test similars
    alltests <- unique(trimws(unlist(lst, recursive = TRUE)))
    df_final <- as.data.frame(
      setNames(lapply(lst, function(a) alltests[ match(alltests, a[,1]) ]),
               sapply(lst, names)),
      stringsAsFactors = FALSE
    )
    #Create a data frame with only element of the column that have
    df_final_NA <- df_final[is.na(df_final$Accesion2),][1]
    

    for (i in df_final_NA$Accesion){
      #print(i)
      duplicate4 <- duplicate4[!duplicate4$Accession == as.character(i), ]
    }
    
    
    # duplicate4$`Peptides(95%)` <- duplicate3$numdup
    # duplicate4$Accession <- gsub("\\;.*","", duplicate4$Accession)
    
    Proteins_PP3 <- data.frame(duplicate4)
    
    Proteins_PP3 <- Proteins_PP3[ order(-Proteins_PP3[,14]), ]
    Numeric_table <- data.frame(N=1:nrow(Proteins_PP3))
    
    Proteins_PP3 <- data.frame("N"=Numeric_table$N, Proteins_PP3)
    
    #Proteins_PP6$Species <- gsub('\\d+\\.?','', Proteins_PP6$Species)

  #Ordenamos los peptidos y prote�nas por accesion
    Proteins_PP4 <- with(Proteins_PP3,  Proteins_PP3[order(Proteins_PP3$Accession) , ])
    Peptidos_PP4 <- with(Peptidos_PP2,  Peptidos_PP2[order(Peptidos_PP2$Accesion) , ])
    
    #Peptidos_PP4$N <- 1:nrow(Peptidos_PP4) 
    
    #Creamos el data frame de la enumeraci�n. 
    test1 <- data.frame("N"=Proteins_PP4$N, "Accession"=Proteins_PP4$Accession)
    test2 <- data.frame("Accession"=Peptidos_PP4$Accesion)
    
    test3 <- data.frame(merge(test2,test1,all.x = T))
    Peptidos_PP4 <- data.frame("N"=test3$N, Peptidos_PP4)
    
    #Ordenamos los peptidos
    Proteins_PP4 <- with(Proteins_PP4,  Proteins_PP4[order(Proteins_PP4$N) , ])
    Peptidos_PP4 <- with(Peptidos_PP4[  Peptidos_PP4[order(Peptidos_PP4$N) , ])

    ########################################
    
    test_value <- data.frame(expand.grid(rep(list(0:1),4))) 
    test_value <- test_value[order(nrow(test_value):1),]
    test_value <- test_value[-nrow(test_value),]
    
    #Clasification. Se filtran los resultados. Podemos hacer hasta 6. 
    
    if(){
      
      test1 <- data.frame(Peptidos_PP4[!is.na(Peptidos_PP4$Ions.Score.A4.Mascot) & !is.na(Peptidos_PP4$Ions.Score.B4.Mascot),])
      test2 <- data.frame(Peptidos_PP4[!is.na(Peptidos_PP4$Ions.Score.A4.Mascot) & is.na(Peptidos_PP4$Ions.Score.B4.Mascot),])
      test3 <- data.frame(Peptidos_PP4[is.na(Peptidos_PP4$Ions.Score.A4.Mascot) & !is.na(Peptidos_PP4$Ions.Score.B4.Mascot),])
      
      Comunes_ALL <- rbind(test1, test2, test3)
      Comunes2x2 <- rbind(test1)
      Comunes1x1 <- rbind(test2,test3)
      
    }else if{
      
      #Filter with the same peptides
      test1 <- data.frame(Peptidos_PP4[  !is.na(Peptidos_PP4$Ions.Score.A4.Mascot) & !is.na(Peptidos_PP4$Ions.Score.B4.Mascot`)& !is.na(Peptidos_PP4$Ions.Score.C4.Mascot),])
      test1 <- with(test1,  test1[order(-test1$Ions.Score.A4.Mascot) , ])
      #Different peptides between 1º and the 2º. 
      test2 <- data.frame(Peptidos_PP4[  !is.na(Peptidos_PP4$Ions.Score.A4.Mascot) & !is.na(Peptidos_PP4$Ions.Score.B4.Mascot`)& is.na(Peptidos_PP4$Ions.Score.C4.Mascot),])
      test2 <- with(test2,  test2[order(-test2$Ions.Score.A4.Mascot) , ])
      #Different peptides between 2º and the 1º.
      test3 <- data.frame(Peptidos_PP4[  is.na(Peptidos_PP4$Ions.Score.A4.Mascot) & !is.na(Peptidos_PP4$Ions.Score.B4.Mascot`)& !is.na(Peptidos_PP4$Ions.Score.C4.Mascot),])
      test3 <- with(test3,  test3[order(-test3$Ions.Score.B4.Mascot`) , ])
      #Different between 3º and 2º,1º
      test4 <- data.frame(Peptidos_PP4[  !is.na(Peptidos_PP4$Ions.Score.A4.Mascot) & is.na(Peptidos_PP4$Ions.Score.B4.Mascot`)& !is.na(Peptidos_PP4$Ions.Score.C4.Mascot),])
      test4 <- with(test4,  test4[order(-test4$Ions.Score.C4.Mascot) , ])
      
      test5 <- data.frame(Peptidos_PP4[  !is.na(Peptidos_PP4$Ions.Score.A4.Mascot) & is.na(Peptidos_PP4$Ions.Score.B4.Mascot`)& is.na(Peptidos_PP4$Ions.Score.C4.Mascot),])
      test5 <- with(test5,  test5[order(-test5$Ions.Score.A4.Mascot) , ])
      
      test6 <- data.frame(Peptidos_PP4[  is.na(Peptidos_PP4$Ions.Score.A4.Mascot) & !is.na(Peptidos_PP4$Ions.Score.B4.Mascot`)& is.na(Peptidos_PP4$Ions.Score.C4.Mascot),])
      test6 <- with(test6,  test6[order(-test6$Ions.Score.B4.Mascot`) , ])
      
      test7 <- data.frame(Peptidos_PP4[  is.na(Peptidos_PP4$Ions.Score.A4.Mascot) &  is.na(Peptidos_PP4$Ions.Score.B4.Mascot`)& !is.na(Peptidos_PP4$Ions.Score.C4.Mascot),])
      test7 <- with(test7,  test7[order(-test7$Ions.Score.C4.Mascot) , ])
      
      #Lista de comunes
      Comunes_ALL <- rbind(test1,test2,test3,test4, test5, test6,test7)
      Comunes1a1 <- rbind(test5, test6,test7)
      Comunes2a2 <- rbind(test2,test3,test4)
      Comunes3a3 <- test1
      
    }else if(){
      
      #4
      
      #Clasification. Se filtran los resultados
      #ALL
      test1 <- data.frame(Peptidos_PP4[ !is.na(Peptidos_PP4$Ions.Score.A4.Mascot)& !is.na(Peptidos_PP4$Ions.Score.B4.Mascot)& !is.na(Peptidos_PP4$Ions.Score.C4.Mascot)& !is.na(Peptidos_PP4$Ions.Score.D4.Mascot),])
      test1 <- with(test1,  test1[order(-test1$Ions.Score.A4.Mascot) , ])
      #Different peptides between 1º and the 2º. 
      test2 <- data.frame(Peptidos_PP4[ !is.na(Peptidos_PP4$Ions.Score.A4.Mascot)& !is.na(Peptidos_PP4$Ions.Score.B4.Mascot)& !is.na(Peptidos_PP4$Ions.Score.C4.Mascot)& is.na(Peptidos_PP4$Ions.Score.D4.Mascot),])
      test2 <- with(test2,  test2[order(-test2$Ions.Score.A4.Mascot) , ])
      #Different peptides between 2º and the 1º.
      test3 <- data.frame(Peptidos_PP4[ !is.na(Peptidos_PP4$Ions.Score.A4.Mascot)& !is.na(Peptidos_PP4$Ions.Score.B4.Mascot)& is.na(Peptidos_PP4$Ions.Score.C4.Mascot)& !is.na(Peptidos_PP4$Ions.Score.D4.Mascot),])
      test3 <- with(test3,  test3[order(-test3$Ions.Score.A4.Mascot) , ])
      #Different between 3º and 2º,1º
      test4 <- data.frame(Peptidos_PP4[ !is.na(Peptidos_PP4$Ions.Score.A4.Mascot)& is.na(Peptidos_PP4$Ions.Score.B4.Mascot)& !is.na(Peptidos_PP4$Ions.Score.C4.Mascot)& !is.na(Peptidos_PP4$Ions.Score.D4.Mascot),])
      test4 <- with(test4,  test4[order(-test4$Ions.Score.B4.Mascot) , ])
      
      test5 <- data.frame(Peptidos_PP4[ is.na(Peptidos_PP4$Ions.Score.A4.Mascot)& !is.na(Peptidos_PP4$Ions.Score.B4.Mascot)& !is.na(Peptidos_PP4$Ions.Score.C4.Mascot)& !is.na(Peptidos_PP4$Ions.Score.D4.Mascot),])        
      test5 <- with(test5,  test5[order(-test5$Ions.Score.B4.Mascot) , ])
      
      test6 <- data.frame(Peptidos_PP4[ !is.na(Peptidos_PP4$Ions.Score.A4.Mascot)& !is.na(Peptidos_PP4$Ions.Score.B4.Mascot)& is.na(Peptidos_PP4$Ions.Score.C4.Mascot)& is.na(Peptidos_PP4$Ions.Score.D4.Mascot),])        
      test6 <- with(test6,  test6[order(-test6$Ions.Score.A4.Mascot) , ])
      
      test7 <- data.frame(Peptidos_PP4[ !is.na(Peptidos_PP4$Ions.Score.A4.Mascot)& is.na(Peptidos_PP4$Ions.Score.B4.Mascot)& !is.na(Peptidos_PP4$Ions.Score.C4.Mascot)& is.na(Peptidos_PP4$Ions.Score.D4.Mascot),])
      test7 <- with(test7,  test7[order(-test7$Ions.Score.A4.Mascot) , ])
      
      test8 <- data.frame(Peptidos_PP4[ !is.na(Peptidos_PP4$Ions.Score.A4.Mascot)& is.na(Peptidos_PP4$Ions.Score.B4.Mascot)& is.na(Peptidos_PP4$Ions.Score.C4.Mascot)& !is.na(Peptidos_PP4$Ions.Score.D4.Mascot),])
      test8 <- with(test8,  test8[order(-test8$Ions.Score.A4.Mascot) , ])
      
      test9 <- data.frame(Peptidos_PP4[ is.na(Peptidos_PP4$Ions.Score.A4.Mascot)& !is.na(Peptidos_PP4$Ions.Score.B4.Mascot)& !is.na(Peptidos_PP4$Ions.Score.C4.Mascot)& is.na(Peptidos_PP4$Ions.Score.D4.Mascot),])
      test9 <- with(test9,  test9[order(-test9$Ions.Score.B4.Mascot) , ])
      
      test10 <- data.frame(Peptidos_PP4[ is.na(Peptidos_PP4$Ions.Score.A4.Mascot)& !is.na(Peptidos_PP4$Ions.Score.B4.Mascot)& is.na(Peptidos_PP4$Ions.Score.C4.Mascot)& !is.na(Peptidos_PP4$Ions.Score.D4.Mascot),])
      test10 <- with(test10,  test10[order(-test10$Ions.Score.B4.Mascot) , ])
      
      test11 <- data.frame(Peptidos_PP4[ is.na(Peptidos_PP4$Ions.Score.A4.Mascot)& is.na(Peptidos_PP4$Ions.Score.B4.Mascot)& !is.na(Peptidos_PP4$Ions.Score.C4.Mascot)& !is.na(Peptidos_PP4$Ions.Score.D4.Mascot),])
      test11 <- with(test11,  test11[order(-test11$Ions.Score.C4.Mascot) , ])
      
      test12 <- data.frame(Peptidos_PP4[ !is.na(Peptidos_PP4$Ions.Score.A4.Mascot)& is.na(Peptidos_PP4$Ions.Score.B4.Mascot)& is.na(Peptidos_PP4$Ions.Score.C4.Mascot)& is.na(Peptidos_PP4$Ions.Score.D4.Mascot),])
      test12 <- with(test12,  test12[order(-test12$Ions.Score.A4.Mascot) , ])
      
      test13 <- data.frame(Peptidos_PP4[ is.na(Peptidos_PP4$Ions.Score.A4.Mascot)& !is.na(Peptidos_PP4$Ions.Score.B4.Mascot)& is.na(Peptidos_PP4$Ions.Score.C4.Mascot)& is.na(Peptidos_PP4$Ions.Score.D4.Mascot),])
      test13 <- with(test13,  test13[order(-test13$Ions.Score.B4.Mascot) , ])
      
      test14 <- data.frame(Peptidos_PP4[ is.na(Peptidos_PP4$Ions.Score.A4.Mascot)& is.na(Peptidos_PP4$Ions.Score.B4.Mascot)& !is.na(Peptidos_PP4$Ions.Score.C4.Mascot)& is.na(Peptidos_PP4$Ions.Score.D4.Mascot),])
      test14 <- with(test14,  test14[order(-test14$Ions.Score.C4.Mascot) , ])
      
      test15 <- data.frame(Peptidos_PP4[ is.na(Peptidos_PP4$Ions.Score.A4.Mascot)& is.na(Peptidos_PP4$Ions.Score.B4.Mascot)& is.na(Peptidos_PP4$Ions.Score.C4.Mascot)& !is.na(Peptidos_PP4$Ions.Score.D4.Mascot),])
      test15 <- with(test15,  test15[order(-test15$Ions.Score.D4.Mascot) , ])
      
      
      #Merge data. 
      Comunes_ALL<-rbind(test1, test2, test3, test4, test5, test6, test7, test8, test9, test10, test11, test12, test13, test14, test15)
      
      Comunes4a4 <- test1
      Comunes3a3 <- rbind(test2,test3,test4,test5)
      Comunes2a2 <- rbind(test6, test7, test8, test9, test10, test11)
      Comunes1a1 <- rbind(test12, test13, test14, test15)
      
      }else if{
      	#5
        #Clasification. Se filtran los resultados
        #ALL
        test1 <- data.frame(Peptidos_PP4[ !is.na(Peptidos_PP4$Ions.Score.A4.Mascot)& !is.na(Peptidos_PP4$Ions.Score.B4.Mascot)& !is.na(Peptidos_PP4$Ions.Score.C4.Mascot)& !is.na(Peptidos_PP4$Ions.Score.D4.Mascot)& !is.na(Peptidos_PP4$Ions.Score.E4.Mascot),])
        test1 <- with(test1,  test1[order(-test1$Ions.Score.A4.Mascot) , ])
        #Different peptides between 1º and the 2º. 
        test2 <- data.frame(Peptidos_PP4[ !is.na(Peptidos_PP4$Ions.Score.A4.Mascot)& !is.na(Peptidos_PP4$Ions.Score.B4.Mascot)& !is.na(Peptidos_PP4$Ions.Score.C4.Mascot)& !is.na(Peptidos_PP4$Ions.Score.D4.Mascot)& is.na(Peptidos_PP4$Ions.Score.E4.Mascot),])
        test2 <- with(test2,  test2[order(-test2$Ions.Score.A4.Mascot) , ])
        #Different peptides between 2º and the 1º.
        test3 <- data.frame(Peptidos_PP4[ !is.na(Peptidos_PP4$Ions.Score.A4.Mascot)& !is.na(Peptidos_PP4$Ions.Score.B4.Mascot)& !is.na(Peptidos_PP4$Ions.Score.C4.Mascot)& is.na(Peptidos_PP4$Ions.Score.D4.Mascot)& !is.na(Peptidos_PP4$Ions.Score.E4.Mascot),])
        test3 <- with(test3,  test3[order(-test3$Ions.Score.A4.Mascot) , ])
        #Different between 3º and 2º,1º
        test4 <- data.frame(Peptidos_PP4[ !is.na(Peptidos_PP4$Ions.Score.A4.Mascot)& !is.na(Peptidos_PP4$Ions.Score.B4.Mascot)& is.na(Peptidos_PP4$Ions.Score.C4.Mascot)& !is.na(Peptidos_PP4$Ions.Score.D4.Mascot)& !is.na(Peptidos_PP4$Ions.Score.E4.Mascot),])
        test4 <- with(test4,  test4[order(-test4$Ions.Score.B4.Mascot) , ])
        
        test5 <- data.frame(Peptidos_PP4[ !is.na(Peptidos_PP4$Ions.Score.A4.Mascot)& is.na(Peptidos_PP4$Ions.Score.B4.Mascot)& !is.na(Peptidos_PP4$Ions.Score.C4.Mascot)& !is.na(Peptidos_PP4$Ions.Score.D4.Mascot)& !is.na(Peptidos_PP4$Ions.Score.E4.Mascot),])        
        test5 <- with(test5,  test5[order(-test5$Ions.Score.B4.Mascot) , ])
        
        test6 <- data.frame(Peptidos_PP4[ is.na(Peptidos_PP4$Ions.Score.A4.Mascot)& !is.na(Peptidos_PP4$Ions.Score.B4.Mascot)& !is.na(Peptidos_PP4$Ions.Score.C4.Mascot)& !is.na(Peptidos_PP4$Ions.Score.D4.Mascot)& !is.na(Peptidos_PP4$Ions.Score.E4.Mascot),])        
        test6 <- with(test6,  test6[order(-test6$Ions.Score.A4.Mascot) , ])
        
        test7 <- data.frame(Peptidos_PP4[ !is.na(Peptidos_PP4$Ions.Score.A4.Mascot)& !is.na(Peptidos_PP4$Ions.Score.B4.Mascot)& !is.na(Peptidos_PP4$Ions.Score.C4.Mascot)& is.na(Peptidos_PP4$Ions.Score.D4.Mascot)& is.na(Peptidos_PP4$Ions.Score.E4.Mascot),])
        test7 <- with(test7,  test7[order(-test7$Ions.Score.A4.Mascot) , ])
        
        test8 <- data.frame(Peptidos_PP4[ !is.na(Peptidos_PP4$Ions.Score.A4.Mascot)& !is.na(Peptidos_PP4$Ions.Score.B4.Mascot)& is.na(Peptidos_PP4$Ions.Score.C4.Mascot)& !is.na(Peptidos_PP4$Ions.Score.D4.Mascot)& is.na(Peptidos_PP4$Ions.Score.E4.Mascot),])
        test8 <- with(test8,  test8[order(-test8$Ions.Score.A4.Mascot) , ])
        
        test9 <- data.frame(Peptidos_PP4[ !is.na(Peptidos_PP4$Ions.Score.A4.Mascot)& !is.na(Peptidos_PP4$Ions.Score.B4.Mascot)& is.na(Peptidos_PP4$Ions.Score.C4.Mascot)& is.na(Peptidos_PP4$Ions.Score.D4.Mascot)& !is.na(Peptidos_PP4$Ions.Score.E4.Mascot),])
        test9 <- with(test9,  test9[order(-test9$Ions.Score.B4.Mascot) , ])
        
        test10 <- data.frame(Peptidos_PP4[ !is.na(Peptidos_PP4$Ions.Score.A4.Mascot)& is.na(Peptidos_PP4$Ions.Score.B4.Mascot)& !is.na(Peptidos_PP4$Ions.Score.C4.Mascot)& !is.na(Peptidos_PP4$Ions.Score.D4.Mascot)& is.na(Peptidos_PP4$Ions.Score.E4.Mascot),])
        test10 <- with(test10,  test10[order(-test10$Ions.Score.B4.Mascot) , ])
        
        test11 <- data.frame(Peptidos_PP4[ !is.na(Peptidos_PP4$Ions.Score.A4.Mascot)& is.na(Peptidos_PP4$Ions.Score.B4.Mascot)& !is.na(Peptidos_PP4$Ions.Score.C4.Mascot)& is.na(Peptidos_PP4$Ions.Score.D4.Mascot)& !is.na(Peptidos_PP4$Ions.Score.E4.Mascot),])
        test11 <- with(test11,  test11[order(-test11$Ions.Score.C4.Mascot) , ])
        
        test12 <- data.frame(Peptidos_PP4[ !is.na(Peptidos_PP4$Ions.Score.A4.Mascot)& is.na(Peptidos_PP4$Ions.Score.B4.Mascot)& is.na(Peptidos_PP4$Ions.Score.C4.Mascot)& !is.na(Peptidos_PP4$Ions.Score.D4.Mascot)& !is.na(Peptidos_PP4$Ions.Score.E4.Mascot),])
        test12 <- with(test12,  test12[order(-test12$Ions.Score.A4.Mascot) , ])
        
        test13 <- data.frame(Peptidos_PP4[ is.na(Peptidos_PP4$Ions.Score.A4.Mascot)& !is.na(Peptidos_PP4$Ions.Score.B4.Mascot)& !is.na(Peptidos_PP4$Ions.Score.C4.Mascot)& !is.na(Peptidos_PP4$Ions.Score.D4.Mascot)& is.na(Peptidos_PP4$Ions.Score.E4.Mascot),])
        test13 <- with(test13,  test13[order(-test13$Ions.Score.B4.Mascot) , ])
        
        test14 <- data.frame(Peptidos_PP4[ is.na(Peptidos_PP4$Ions.Score.A4.Mascot)& !is.na(Peptidos_PP4$Ions.Score.B4.Mascot)& is.na(Peptidos_PP4$Ions.Score.C4.Mascot)& !is.na(Peptidos_PP4$Ions.Score.D4.Mascot)& !is.na(Peptidos_PP4$Ions.Score.E4.Mascot),])
        test14 <- with(test14,  test14[order(-test14$Ions.Score.C4.Mascot) , ])
        
        test15 <- data.frame(Peptidos_PP4[ is.na(Peptidos_PP4$Ions.Score.A4.Mascot)& !is.na(Peptidos_PP4$Ions.Score.B4.Mascot)& !is.na(Peptidos_PP4$Ions.Score.C4.Mascot)& is.na(Peptidos_PP4$Ions.Score.D4.Mascot)& !is.na(Peptidos_PP4$Ions.Score.E4.Mascot),])
        test15 <- with(test15,  test15[order(-test15$Ions.Score.D4.Mascot) , ])
        
        test16 <- data.frame(Peptidos_PP4[ is.na(Peptidos_PP4$Ions.Score.A4.Mascot)& is.na(Peptidos_PP4$Ions.Score.B4.Mascot)& !is.na(Peptidos_PP4$Ions.Score.C4.Mascot)& !is.na(Peptidos_PP4$Ions.Score.D4.Mascot)& !is.na(Peptidos_PP4$Ions.Score.E4.Mascot),])
        test16 <- with(test16,  test16[order(-test16$Ions.Score.D4.Mascot) , ])
        
        test17 <- data.frame(Peptidos_PP4[ !is.na(Peptidos_PP4$Ions.Score.A4.Mascot)& is.na(Peptidos_PP4$Ions.Score.B4.Mascot)& is.na(Peptidos_PP4$Ions.Score.C4.Mascot)& !is.na(Peptidos_PP4$Ions.Score.D4.Mascot)& is.na(Peptidos_PP4$Ions.Score.E4.Mascot),])
        test17 <- with(test17,  test17[order(-test17$Ions.Score.D4.Mascot) , ])
        
        test18 <- data.frame(Peptidos_PP4[ !is.na(Peptidos_PP4$Ions.Score.A4.Mascot)& is.na(Peptidos_PP4$Ions.Score.B4.Mascot)& is.na(Peptidos_PP4$Ions.Score.C4.Mascot)& is.na(Peptidos_PP4$Ions.Score.D4.Mascot)& !is.na(Peptidos_PP4$Ions.Score.E4.Mascot),])
        test18 <- with(test18,  test18[order(-test18$Ions.Score.D4.Mascot) , ])
        
        test19 <- data.frame(Peptidos_PP4[ !is.na(Peptidos_PP4$Ions.Score.A4.Mascot)& is.na(Peptidos_PP4$Ions.Score.B4.Mascot)& !is.na(Peptidos_PP4$Ions.Score.C4.Mascot)& is.na(Peptidos_PP4$Ions.Score.D4.Mascot)& is.na(Peptidos_PP4$Ions.Score.E4.Mascot),])
        test19 <- with(test19,  test19[order(-test19$Ions.Score.D4.Mascot) , ])
        
        test20 <- data.frame(Peptidos_PP4[ !is.na(Peptidos_PP4$Ions.Score.A4.Mascot)& !is.na(Peptidos_PP4$Ions.Score.B4.Mascot)& is.na(Peptidos_PP4$Ions.Score.C4.Mascot)& is.na(Peptidos_PP4$Ions.Score.D4.Mascot)& is.na(Peptidos_PP4$Ions.Score.E4.Mascot),])
        test20 <- with(test20,  test20[order(-test20$Ions.Score.D4.Mascot) , ])
        
        test21 <- data.frame(Peptidos_PP4[ is.na(Peptidos_PP4$Ions.Score.A4.Mascot)& !is.na(Peptidos_PP4$Ions.Score.B4.Mascot)& is.na(Peptidos_PP4$Ions.Score.C4.Mascot)& !is.na(Peptidos_PP4$Ions.Score.D4.Mascot)& is.na(Peptidos_PP4$Ions.Score.E4.Mascot),])
        test21 <- with(test21,  test21[order(-test21$Ions.Score.D4.Mascot) , ])
        
        test22 <- data.frame(Peptidos_PP4[ is.na(Peptidos_PP4$Ions.Score.A4.Mascot)& !is.na(Peptidos_PP4$Ions.Score.B4.Mascot)& !is.na(Peptidos_PP4$Ions.Score.C4.Mascot)& is.na(Peptidos_PP4$Ions.Score.D4.Mascot)& is.na(Peptidos_PP4$Ions.Score.E4.Mascot),])
        test22 <- with(test22,  test22[order(-test22$Ions.Score.D4.Mascot) , ])
        
        test23 <- data.frame(Peptidos_PP4[ is.na(Peptidos_PP4$Ions.Score.A4.Mascot)& !is.na(Peptidos_PP4$Ions.Score.B4.Mascot)& is.na(Peptidos_PP4$Ions.Score.C4.Mascot)& is.na(Peptidos_PP4$Ions.Score.D4.Mascot)& !is.na(Peptidos_PP4$Ions.Score.E4.Mascot),])
        test23 <- with(test23,  test23[order(-test23$Ions.Score.D4.Mascot) , ])
        
        test24 <- data.frame(Peptidos_PP4[ is.na(Peptidos_PP4$Ions.Score.A4.Mascot)& is.na(Peptidos_PP4$Ions.Score.B4.Mascot)& !is.na(Peptidos_PP4$Ions.Score.C4.Mascot)& is.na(Peptidos_PP4$Ions.Score.D4.Mascot)& !is.na(Peptidos_PP4$Ions.Score.E4.Mascot),])
        test24 <- with(test24,  test24[order(-test24$Ions.Score.D4.Mascot) , ])
        
        test25 <- data.frame(Peptidos_PP4[ is.na(Peptidos_PP4$Ions.Score.A4.Mascot)& is.na(Peptidos_PP4$Ions.Score.B4.Mascot)& !is.na(Peptidos_PP4$Ions.Score.C4.Mascot)& !is.na(Peptidos_PP4$Ions.Score.D4.Mascot)& is.na(Peptidos_PP4$Ions.Score.E4.Mascot),])
        test25 <- with(test25,  test25[order(-test25$Ions.Score.D4.Mascot) , ])
        
        test26 <- data.frame(Peptidos_PP4[ is.na(Peptidos_PP4$Ions.Score.A4.Mascot)& is.na(Peptidos_PP4$Ions.Score.B4.Mascot)& is.na(Peptidos_PP4$Ions.Score.C4.Mascot)& !is.na(Peptidos_PP4$Ions.Score.D4.Mascot)& !is.na(Peptidos_PP4$Ions.Score.E4.Mascot),])
        test26 <- with(test26,  test26[order(-test26$Ions.Score.D4.Mascot) , ])
        
        test27 <- data.frame(Peptidos_PP4[ !is.na(Peptidos_PP4$Ions.Score.A4.Mascot)& is.na(Peptidos_PP4$Ions.Score.B4.Mascot)& is.na(Peptidos_PP4$Ions.Score.C4.Mascot)& is.na(Peptidos_PP4$Ions.Score.D4.Mascot)& is.na(Peptidos_PP4$Ions.Score.E4.Mascot),])
        test27 <- with(test27,  test27[order(-test27$Ions.Score.D4.Mascot) , ])
        
        test28 <- data.frame(Peptidos_PP4[ is.na(Peptidos_PP4$Ions.Score.A4.Mascot)& !is.na(Peptidos_PP4$Ions.Score.B4.Mascot)& is.na(Peptidos_PP4$Ions.Score.C4.Mascot)& is.na(Peptidos_PP4$Ions.Score.D4.Mascot)& is.na(Peptidos_PP4$Ions.Score.E4.Mascot),])
        test28 <- with(test28,  test28[order(-test28$Ions.Score.D4.Mascot) , ])
        
        test29 <- data.frame(Peptidos_PP4[ is.na(Peptidos_PP4$Ions.Score.A4.Mascot)& is.na(Peptidos_PP4$Ions.Score.B4.Mascot)& !is.na(Peptidos_PP4$Ions.Score.C4.Mascot)& is.na(Peptidos_PP4$Ions.Score.D4.Mascot)& is.na(Peptidos_PP4$Ions.Score.E4.Mascot),])
        test29 <- with(test29,  test29[order(-test29$Ions.Score.D4.Mascot) , ])
        
        test30 <- data.frame(Peptidos_PP4[ is.na(Peptidos_PP4$Ions.Score.A4.Mascot)& is.na(Peptidos_PP4$Ions.Score.B4.Mascot)& is.na(Peptidos_PP4$Ions.Score.C4.Mascot)& !is.na(Peptidos_PP4$Ions.Score.D4.Mascot)& is.na(Peptidos_PP4$Ions.Score.E4.Mascot),])
        test30 <- with(test30,  test30[order(-test30$Ions.Score.D4.Mascot) , ])
        
        test31 <- data.frame(Peptidos_PP4[ is.na(Peptidos_PP4$Ions.Score.A4.Mascot)& is.na(Peptidos_PP4$Ions.Score.B4.Mascot)& is.na(Peptidos_PP4$Ions.Score.C4.Mascot)& is.na(Peptidos_PP4$Ions.Score.D4.Mascot)& !is.na(Peptidos_PP4$Ions.Score.E4.Mascot),])
        test31 <- with(test31,  test31[order(-test31$Ions.Score.D4.Mascot) , ])
        
        #Merge data. 
        Comunes_ALL<-rbind(test1, test2, test3, test4, test5, test6, test7, test8, 
                           test9, test10, test11, test12, test13, test14, test15,test16,test17,test18,test19,test20,
                           test21,test22,test23,test24,test25,test26,test27,test28,test29,test30,test31)
        
        Comunes5a5 <- test1
        Comunes4a4 <- rbind(test2,test3,test4,test5,test6)
        Comunes3a3 <- rbind(test7, test8, test9, test10, test11, test12, test13, test14, test15, test16,test17)
        Comunes2a2 <- rbind(test18, test19, test20, test21, test22, test23,test24,test25,test26)
      	
      }else if{

      	#5
        #Clasification. Se filtran los resultados
        #ALL
        test1 <- data.frame(Peptidos_PP4[ !is.na(Peptidos_PP4$Ions.Score.A4.Mascot)& !is.na(Peptidos_PP4$Ions.Score.B4.Mascot)& !is.na(Peptidos_PP4$Ions.Score.C4.Mascot)& !is.na(Peptidos_PP4$Ions.Score.D4.Mascot)& !is.na(Peptidos_PP4$Ions.Score.E4.Mascot)& !is.na(Peptidos_PP4$Ions.Score.F4.Mascot),])
        test1 <- with(test1,  test1[order(-test1$Ions.Score.A4.Mascot) , ]) #6
        #Different peptides between 1º and the 2º. 
        test2 <- data.frame(Peptidos_PP4[ is.na(Peptidos_PP4$Ions.Score.A4.Mascot)& is.na(Peptidos_PP4$Ions.Score.B4.Mascot)& is.na(Peptidos_PP4$Ions.Score.C4.Mascot)& is.na(Peptidos_PP4$Ions.Score.D4.Mascot)& is.na(Peptidos_PP4$Ions.Score.E4.Mascot)& is.na(Peptidos_PP4$Ions.Score.F4.Mascot),])
        test2 <- with(test2,  test2[order(-test2$Ions.Score.A4.Mascot) , ])
        #Different peptides between 2º and the 1º.
        test3 <- data.frame(Peptidos_PP4[ is.na(Peptidos_PP4$Ions.Score.A4.Mascot)& is.na(Peptidos_PP4$Ions.Score.B4.Mascot)& is.na(Peptidos_PP4$Ions.Score.C4.Mascot)& is.na(Peptidos_PP4$Ions.Score.D4.Mascot)& is.na(Peptidos_PP4$Ions.Score.E4.Mascot)& is.na(Peptidos_PP4$Ions.Score.F4.Mascot),])
        test3 <- with(test3,  test3[order(-test3$Ions.Score.A4.Mascot) , ])
        #Different between 3º and 2º,1º
        test4 <- data.frame(Peptidos_PP4[ is.na(Peptidos_PP4$Ions.Score.A4.Mascot)& is.na(Peptidos_PP4$Ions.Score.B4.Mascot)& is.na(Peptidos_PP4$Ions.Score.C4.Mascot)& is.na(Peptidos_PP4$Ions.Score.D4.Mascot)& is.na(Peptidos_PP4$Ions.Score.E4.Mascot)& is.na(Peptidos_PP4$Ions.Score.F4.Mascot),])
        test4 <- with(test4,  test4[order(-test4$Ions.Score.B4.Mascot) , ])
        
        test5 <- data.frame(Peptidos_PP4[ is.na(Peptidos_PP4$Ions.Score.A4.Mascot)& is.na(Peptidos_PP4$Ions.Score.B4.Mascot)& is.na(Peptidos_PP4$Ions.Score.C4.Mascot)& is.na(Peptidos_PP4$Ions.Score.D4.Mascot)& is.na(Peptidos_PP4$Ions.Score.E4.Mascot)& is.na(Peptidos_PP4$Ions.Score.F4.Mascot),])        
        test5 <- with(test5,  test5[order(-test5$Ions.Score.B4.Mascot) , ])
        
        test6 <- data.frame(Peptidos_PP4[ is.na(Peptidos_PP4$Ions.Score.A4.Mascot)& is.na(Peptidos_PP4$Ions.Score.B4.Mascot)& is.na(Peptidos_PP4$Ions.Score.C4.Mascot)& is.na(Peptidos_PP4$Ions.Score.D4.Mascot)& is.na(Peptidos_PP4$Ions.Score.E4.Mascot)& is.na(Peptidos_PP4$Ions.Score.F4.Mascot),])        
        test6 <- with(test6,  test6[order(-test6$Ions.Score.A4.Mascot) , ])
        
        test7 <- data.frame(Peptidos_PP4[ is.na(Peptidos_PP4$Ions.Score.A4.Mascot)& is.na(Peptidos_PP4$Ions.Score.B4.Mascot)& is.na(Peptidos_PP4$Ions.Score.C4.Mascot)& is.na(Peptidos_PP4$Ions.Score.D4.Mascot)& is.na(Peptidos_PP4$Ions.Score.E4.Mascot)& is.na(Peptidos_PP4$Ions.Score.F4.Mascot),])
        test7 <- with(test7,  test7[order(-test7$Ions.Score.A4.Mascot) , ])
        
        test8 <- data.frame(Peptidos_PP4[ is.na(Peptidos_PP4$Ions.Score.A4.Mascot)& is.na(Peptidos_PP4$Ions.Score.B4.Mascot)& is.na(Peptidos_PP4$Ions.Score.C4.Mascot)& is.na(Peptidos_PP4$Ions.Score.D4.Mascot)& is.na(Peptidos_PP4$Ions.Score.E4.Mascot)& is.na(Peptidos_PP4$Ions.Score.F4.Mascot),])
        test8 <- with(test8,  test8[order(-test8$Ions.Score.A4.Mascot) , ])
        
        test9 <- data.frame(Peptidos_PP4[ is.na(Peptidos_PP4$Ions.Score.A4.Mascot)& is.na(Peptidos_PP4$Ions.Score.B4.Mascot)& is.na(Peptidos_PP4$Ions.Score.C4.Mascot)& is.na(Peptidos_PP4$Ions.Score.D4.Mascot)& is.na(Peptidos_PP4$Ions.Score.E4.Mascot)& is.na(Peptidos_PP4$Ions.Score.F4.Mascot),])
        test9 <- with(test9,  test9[order(-test9$Ions.Score.B4.Mascot) , ])
        
        test10 <- data.frame(Peptidos_PP4[ is.na(Peptidos_PP4$Ions.Score.A4.Mascot)& is.na(Peptidos_PP4$Ions.Score.B4.Mascot)& is.na(Peptidos_PP4$Ions.Score.C4.Mascot)& is.na(Peptidos_PP4$Ions.Score.D4.Mascot)& is.na(Peptidos_PP4$Ions.Score.E4.Mascot)& is.na(Peptidos_PP4$Ions.Score.F4.Mascot),])
        test10 <- with(test10,  test10[order(-test10$Ions.Score.B4.Mascot) , ])
        
        test11 <- data.frame(Peptidos_PP4[ is.na(Peptidos_PP4$Ions.Score.A4.Mascot)& is.na(Peptidos_PP4$Ions.Score.B4.Mascot)& is.na(Peptidos_PP4$Ions.Score.C4.Mascot)& is.na(Peptidos_PP4$Ions.Score.D4.Mascot)& is.na(Peptidos_PP4$Ions.Score.E4.Mascot)& is.na(Peptidos_PP4$Ions.Score.F4.Mascot),])
        test11 <- with(test11,  test11[order(-test11$Ions.Score.C4.Mascot) , ])
        
        test12 <- data.frame(Peptidos_PP4[ is.na(Peptidos_PP4$Ions.Score.A4.Mascot)& is.na(Peptidos_PP4$Ions.Score.B4.Mascot)& is.na(Peptidos_PP4$Ions.Score.C4.Mascot)& is.na(Peptidos_PP4$Ions.Score.D4.Mascot)& is.na(Peptidos_PP4$Ions.Score.E4.Mascot)& is.na(Peptidos_PP4$Ions.Score.F4.Mascot),])
        test12 <- with(test12,  test12[order(-test12$Ions.Score.A4.Mascot) , ])
        
        test13 <- data.frame(Peptidos_PP4[ is.na(Peptidos_PP4$Ions.Score.A4.Mascot)& is.na(Peptidos_PP4$Ions.Score.B4.Mascot)& is.na(Peptidos_PP4$Ions.Score.C4.Mascot)& is.na(Peptidos_PP4$Ions.Score.D4.Mascot)& is.na(Peptidos_PP4$Ions.Score.E4.Mascot)& is.na(Peptidos_PP4$Ions.Score.F4.Mascot),])
        test13 <- with(test13,  test13[order(-test13$Ions.Score.B4.Mascot) , ])
        
        test14 <- data.frame(Peptidos_PP4[ is.na(Peptidos_PP4$Ions.Score.A4.Mascot)& is.na(Peptidos_PP4$Ions.Score.B4.Mascot)& is.na(Peptidos_PP4$Ions.Score.C4.Mascot)& is.na(Peptidos_PP4$Ions.Score.D4.Mascot)& is.na(Peptidos_PP4$Ions.Score.E4.Mascot)& is.na(Peptidos_PP4$Ions.Score.F4.Mascot),])
        test14 <- with(test14,  test14[order(-test14$Ions.Score.C4.Mascot) , ])
        
        test15 <- data.frame(Peptidos_PP4[ is.na(Peptidos_PP4$Ions.Score.A4.Mascot)& is.na(Peptidos_PP4$Ions.Score.B4.Mascot)& is.na(Peptidos_PP4$Ions.Score.C4.Mascot)& is.na(Peptidos_PP4$Ions.Score.D4.Mascot)& is.na(Peptidos_PP4$Ions.Score.E4.Mascot)& is.na(Peptidos_PP4$Ions.Score.F4.Mascot),])
        test15 <- with(test15,  test15[order(-test15$Ions.Score.D4.Mascot) , ])
        
        test16 <- data.frame(Peptidos_PP4[ is.na(Peptidos_PP4$Ions.Score.A4.Mascot)& is.na(Peptidos_PP4$Ions.Score.B4.Mascot)& is.na(Peptidos_PP4$Ions.Score.C4.Mascot)& is.na(Peptidos_PP4$Ions.Score.D4.Mascot)& is.na(Peptidos_PP4$Ions.Score.E4.Mascot)& is.na(Peptidos_PP4$Ions.Score.F4.Mascot),])
        test16 <- with(test16,  test16[order(-test16$Ions.Score.D4.Mascot) , ])
        
        test17 <- data.frame(Peptidos_PP4[ is.na(Peptidos_PP4$Ions.Score.A4.Mascot)& is.na(Peptidos_PP4$Ions.Score.B4.Mascot)& is.na(Peptidos_PP4$Ions.Score.C4.Mascot)& is.na(Peptidos_PP4$Ions.Score.D4.Mascot)& is.na(Peptidos_PP4$Ions.Score.E4.Mascot)& is.na(Peptidos_PP4$Ions.Score.F4.Mascot),])
        test17 <- with(test17,  test17[order(-test17$Ions.Score.D4.Mascot) , ])
        
        test18 <- data.frame(Peptidos_PP4[ is.na(Peptidos_PP4$Ions.Score.A4.Mascot)& is.na(Peptidos_PP4$Ions.Score.B4.Mascot)& is.na(Peptidos_PP4$Ions.Score.C4.Mascot)& is.na(Peptidos_PP4$Ions.Score.D4.Mascot)& is.na(Peptidos_PP4$Ions.Score.E4.Mascot)& is.na(Peptidos_PP4$Ions.Score.F4.Mascot),])
        test18 <- with(test18,  test18[order(-test18$Ions.Score.D4.Mascot) , ])
        
        test19 <- data.frame(Peptidos_PP4[ is.na(Peptidos_PP4$Ions.Score.A4.Mascot)& is.na(Peptidos_PP4$Ions.Score.B4.Mascot)& is.na(Peptidos_PP4$Ions.Score.C4.Mascot)& is.na(Peptidos_PP4$Ions.Score.D4.Mascot)& is.na(Peptidos_PP4$Ions.Score.E4.Mascot)& is.na(Peptidos_PP4$Ions.Score.F4.Mascot),])
        test19 <- with(test19,  test19[order(-test19$Ions.Score.D4.Mascot) , ])
        
        test20 <- data.frame(Peptidos_PP4[ is.na(Peptidos_PP4$Ions.Score.A4.Mascot)& is.na(Peptidos_PP4$Ions.Score.B4.Mascot)& is.na(Peptidos_PP4$Ions.Score.C4.Mascot)& is.na(Peptidos_PP4$Ions.Score.D4.Mascot)& is.na(Peptidos_PP4$Ions.Score.E4.Mascot)& is.na(Peptidos_PP4$Ions.Score.F4.Mascot),])
        test20 <- with(test20,  test20[order(-test20$Ions.Score.D4.Mascot) , ])
        
        test21 <- data.frame(Peptidos_PP4[ is.na(Peptidos_PP4$Ions.Score.A4.Mascot)& is.na(Peptidos_PP4$Ions.Score.B4.Mascot)& is.na(Peptidos_PP4$Ions.Score.C4.Mascot)& is.na(Peptidos_PP4$Ions.Score.D4.Mascot)& is.na(Peptidos_PP4$Ions.Score.E4.Mascot)& is.na(Peptidos_PP4$Ions.Score.F4.Mascot),])
        test21 <- with(test21,  test21[order(-test21$Ions.Score.D4.Mascot) , ])
        
        test22 <- data.frame(Peptidos_PP4[ is.na(Peptidos_PP4$Ions.Score.A4.Mascot)& is.na(Peptidos_PP4$Ions.Score.B4.Mascot)& is.na(Peptidos_PP4$Ions.Score.C4.Mascot)& is.na(Peptidos_PP4$Ions.Score.D4.Mascot)& is.na(Peptidos_PP4$Ions.Score.E4.Mascot)& is.na(Peptidos_PP4$Ions.Score.F4.Mascot),])
        test22 <- with(test22,  test22[order(-test22$Ions.Score.D4.Mascot) , ])
        
        test23 <- data.frame(Peptidos_PP4[ is.na(Peptidos_PP4$Ions.Score.A4.Mascot)& is.na(Peptidos_PP4$Ions.Score.B4.Mascot)& is.na(Peptidos_PP4$Ions.Score.C4.Mascot)& is.na(Peptidos_PP4$Ions.Score.D4.Mascot)& is.na(Peptidos_PP4$Ions.Score.E4.Mascot)& is.na(Peptidos_PP4$Ions.Score.F4.Mascot),])
        test23 <- with(test23,  test23[order(-test23$Ions.Score.D4.Mascot) , ])
        
        test24 <- data.frame(Peptidos_PP4[ is.na(Peptidos_PP4$Ions.Score.A4.Mascot)& is.na(Peptidos_PP4$Ions.Score.B4.Mascot)& is.na(Peptidos_PP4$Ions.Score.C4.Mascot)& is.na(Peptidos_PP4$Ions.Score.D4.Mascot)& is.na(Peptidos_PP4$Ions.Score.E4.Mascot)& is.na(Peptidos_PP4$Ions.Score.F4.Mascot),])
        test24 <- with(test24,  test24[order(-test24$Ions.Score.D4.Mascot) , ])
        
        test25 <- data.frame(Peptidos_PP4[ is.na(Peptidos_PP4$Ions.Score.A4.Mascot)& is.na(Peptidos_PP4$Ions.Score.B4.Mascot)& is.na(Peptidos_PP4$Ions.Score.C4.Mascot)& is.na(Peptidos_PP4$Ions.Score.D4.Mascot)& is.na(Peptidos_PP4$Ions.Score.E4.Mascot)& is.na(Peptidos_PP4$Ions.Score.F4.Mascot),])
        test25 <- with(test25,  test25[order(-test25$Ions.Score.D4.Mascot) , ])
        
        test26 <- data.frame(Peptidos_PP4[ is.na(Peptidos_PP4$Ions.Score.A4.Mascot)& is.na(Peptidos_PP4$Ions.Score.B4.Mascot)& is.na(Peptidos_PP4$Ions.Score.C4.Mascot)& is.na(Peptidos_PP4$Ions.Score.D4.Mascot)& is.na(Peptidos_PP4$Ions.Score.E4.Mascot)& is.na(Peptidos_PP4$Ions.Score.F4.Mascot),])
        test26 <- with(test26,  test26[order(-test26$Ions.Score.D4.Mascot) , ])
        
        test27 <- data.frame(Peptidos_PP4[ is.na(Peptidos_PP4$Ions.Score.A4.Mascot)& is.na(Peptidos_PP4$Ions.Score.B4.Mascot)& is.na(Peptidos_PP4$Ions.Score.C4.Mascot)& is.na(Peptidos_PP4$Ions.Score.D4.Mascot)& is.na(Peptidos_PP4$Ions.Score.E4.Mascot)& is.na(Peptidos_PP4$Ions.Score.F4.Mascot),])
        test27 <- with(test27,  test27[order(-test27$Ions.Score.D4.Mascot) , ])
        
        test28 <- data.frame(Peptidos_PP4[ is.na(Peptidos_PP4$Ions.Score.A4.Mascot)& is.na(Peptidos_PP4$Ions.Score.B4.Mascot)& is.na(Peptidos_PP4$Ions.Score.C4.Mascot)& is.na(Peptidos_PP4$Ions.Score.D4.Mascot)& is.na(Peptidos_PP4$Ions.Score.E4.Mascot)& is.na(Peptidos_PP4$Ions.Score.F4.Mascot),])
        test28 <- with(test28,  test28[order(-test28$Ions.Score.D4.Mascot) , ])
        
        test29 <- data.frame(Peptidos_PP4[ is.na(Peptidos_PP4$Ions.Score.A4.Mascot)& is.na(Peptidos_PP4$Ions.Score.B4.Mascot)& is.na(Peptidos_PP4$Ions.Score.C4.Mascot)& is.na(Peptidos_PP4$Ions.Score.D4.Mascot)& is.na(Peptidos_PP4$Ions.Score.E4.Mascot)& is.na(Peptidos_PP4$Ions.Score.F4.Mascot),])
        test29 <- with(test29,  test29[order(-test29$Ions.Score.D4.Mascot) , ])

        test30 <- data.frame(Peptidos_PP4[ is.na(Peptidos_PP4$Ions.Score.A4.Mascot)& is.na(Peptidos_PP4$Ions.Score.B4.Mascot)& is.na(Peptidos_PP4$Ions.Score.C4.Mascot)& is.na(Peptidos_PP4$Ions.Score.D4.Mascot)& is.na(Peptidos_PP4$Ions.Score.E4.Mascot)& is.na(Peptidos_PP4$Ions.Score.F4.Mascot),])

        test30 <- with(test30,  test30[order(-test30$Ions.Score.D4.Mascot) , ])

        test31 <- data.frame(Peptidos_PP4[ is.na(Peptidos_PP4$Ions.Score.A4.Mascot)& is.na(Peptidos_PP4$Ions.Score.B4.Mascot)& is.na(Peptidos_PP4$Ions.Score.C4.Mascot)& is.na(Peptidos_PP4$Ions.Score.D4.Mascot)& is.na(Peptidos_PP4$Ions.Score.E4.Mascot)& is.na(Peptidos_PP4$Ions.Score.F4.Mascot),])

        test31 <- with(test31,  test31[order(-test31$Ions.Score.D4.Mascot) , ])

        test32 <- data.frame(Peptidos_PP4[ is.na(Peptidos_PP4$Ions.Score.A4.Mascot)& is.na(Peptidos_PP4$Ions.Score.B4.Mascot)& is.na(Peptidos_PP4$Ions.Score.C4.Mascot)& is.na(Peptidos_PP4$Ions.Score.D4.Mascot)& is.na(Peptidos_PP4$Ions.Score.E4.Mascot)& is.na(Peptidos_PP4$Ions.Score.F4.Mascot),])

        test33 <- data.frame(Peptidos_PP4[ is.na(Peptidos_PP4$Ions.Score.A4.Mascot)& is.na(Peptidos_PP4$Ions.Score.B4.Mascot)& is.na(Peptidos_PP4$Ions.Score.C4.Mascot)& is.na(Peptidos_PP4$Ions.Score.D4.Mascot)& is.na(Peptidos_PP4$Ions.Score.E4.Mascot)& is.na(Peptidos_PP4$Ions.Score.F4.Mascot),])

        test34 <- data.frame(Peptidos_PP4[ is.na(Peptidos_PP4$Ions.Score.A4.Mascot)& is.na(Peptidos_PP4$Ions.Score.B4.Mascot)& is.na(Peptidos_PP4$Ions.Score.C4.Mascot)& is.na(Peptidos_PP4$Ions.Score.D4.Mascot)& is.na(Peptidos_PP4$Ions.Score.E4.Mascot)& is.na(Peptidos_PP4$Ions.Score.F4.Mascot),])

        test35 <- data.frame(Peptidos_PP4[ is.na(Peptidos_PP4$Ions.Score.A4.Mascot)& is.na(Peptidos_PP4$Ions.Score.B4.Mascot)& is.na(Peptidos_PP4$Ions.Score.C4.Mascot)& is.na(Peptidos_PP4$Ions.Score.D4.Mascot)& is.na(Peptidos_PP4$Ions.Score.E4.Mascot)& is.na(Peptidos_PP4$Ions.Score.F4.Mascot),])

        test36 <- data.frame(Peptidos_PP4[ is.na(Peptidos_PP4$Ions.Score.A4.Mascot)& is.na(Peptidos_PP4$Ions.Score.B4.Mascot)& is.na(Peptidos_PP4$Ions.Score.C4.Mascot)& is.na(Peptidos_PP4$Ions.Score.D4.Mascot)& is.na(Peptidos_PP4$Ions.Score.E4.Mascot)& is.na(Peptidos_PP4$Ions.Score.F4.Mascot),])

        test37 <- data.frame(Peptidos_PP4[ is.na(Peptidos_PP4$Ions.Score.A4.Mascot)& is.na(Peptidos_PP4$Ions.Score.B4.Mascot)& is.na(Peptidos_PP4$Ions.Score.C4.Mascot)& is.na(Peptidos_PP4$Ions.Score.D4.Mascot)& is.na(Peptidos_PP4$Ions.Score.E4.Mascot)& is.na(Peptidos_PP4$Ions.Score.F4.Mascot),])

        test38 <- data.frame(Peptidos_PP4[ is.na(Peptidos_PP4$Ions.Score.A4.Mascot)& is.na(Peptidos_PP4$Ions.Score.B4.Mascot)& is.na(Peptidos_PP4$Ions.Score.C4.Mascot)& is.na(Peptidos_PP4$Ions.Score.D4.Mascot)& is.na(Peptidos_PP4$Ions.Score.E4.Mascot)& is.na(Peptidos_PP4$Ions.Score.F4.Mascot),])

        test39 <- data.frame(Peptidos_PP4[ is.na(Peptidos_PP4$Ions.Score.A4.Mascot)& is.na(Peptidos_PP4$Ions.Score.B4.Mascot)& is.na(Peptidos_PP4$Ions.Score.C4.Mascot)& is.na(Peptidos_PP4$Ions.Score.D4.Mascot)& is.na(Peptidos_PP4$Ions.Score.E4.Mascot)& is.na(Peptidos_PP4$Ions.Score.F4.Mascot),])

        test40 <- data.frame(Peptidos_PP4[ is.na(Peptidos_PP4$Ions.Score.A4.Mascot)& is.na(Peptidos_PP4$Ions.Score.B4.Mascot)& is.na(Peptidos_PP4$Ions.Score.C4.Mascot)& is.na(Peptidos_PP4$Ions.Score.D4.Mascot)& is.na(Peptidos_PP4$Ions.Score.E4.Mascot)& is.na(Peptidos_PP4$Ions.Score.F4.Mascot),])

        test41 <- data.frame(Peptidos_PP4[ is.na(Peptidos_PP4$Ions.Score.A4.Mascot)& is.na(Peptidos_PP4$Ions.Score.B4.Mascot)& is.na(Peptidos_PP4$Ions.Score.C4.Mascot)& is.na(Peptidos_PP4$Ions.Score.D4.Mascot)& is.na(Peptidos_PP4$Ions.Score.E4.Mascot)& is.na(Peptidos_PP4$Ions.Score.F4.Mascot),])

        test42<- data.frame(Peptidos_PP4[ is.na(Peptidos_PP4$Ions.Score.A4.Mascot)& is.na(Peptidos_PP4$Ions.Score.B4.Mascot)& is.na(Peptidos_PP4$Ions.Score.C4.Mascot)& is.na(Peptidos_PP4$Ions.Score.D4.Mascot)& is.na(Peptidos_PP4$Ions.Score.E4.Mascot)& is.na(Peptidos_PP4$Ions.Score.F4.Mascot),])

        test43 <- data.frame(Peptidos_PP4[ is.na(Peptidos_PP4$Ions.Score.A4.Mascot)& is.na(Peptidos_PP4$Ions.Score.B4.Mascot)& is.na(Peptidos_PP4$Ions.Score.C4.Mascot)& is.na(Peptidos_PP4$Ions.Score.D4.Mascot)& is.na(Peptidos_PP4$Ions.Score.E4.Mascot)& is.na(Peptidos_PP4$Ions.Score.F4.Mascot),])

        test44 <- data.frame(Peptidos_PP4[ is.na(Peptidos_PP4$Ions.Score.A4.Mascot)& is.na(Peptidos_PP4$Ions.Score.B4.Mascot)& is.na(Peptidos_PP4$Ions.Score.C4.Mascot)& is.na(Peptidos_PP4$Ions.Score.D4.Mascot)& is.na(Peptidos_PP4$Ions.Score.E4.Mascot)& is.na(Peptidos_PP4$Ions.Score.F4.Mascot),])

        test45 <- data.frame(Peptidos_PP4[ is.na(Peptidos_PP4$Ions.Score.A4.Mascot)& is.na(Peptidos_PP4$Ions.Score.B4.Mascot)& is.na(Peptidos_PP4$Ions.Score.C4.Mascot)& is.na(Peptidos_PP4$Ions.Score.D4.Mascot)& is.na(Peptidos_PP4$Ions.Score.E4.Mascot)& is.na(Peptidos_PP4$Ions.Score.F4.Mascot),])

        test46 <- data.frame(Peptidos_PP4[ is.na(Peptidos_PP4$Ions.Score.A4.Mascot)& is.na(Peptidos_PP4$Ions.Score.B4.Mascot)& is.na(Peptidos_PP4$Ions.Score.C4.Mascot)& is.na(Peptidos_PP4$Ions.Score.D4.Mascot)& is.na(Peptidos_PP4$Ions.Score.E4.Mascot)& is.na(Peptidos_PP4$Ions.Score.F4.Mascot),])

        test47 <- data.frame(Peptidos_PP4[ is.na(Peptidos_PP4$Ions.Score.A4.Mascot)& is.na(Peptidos_PP4$Ions.Score.B4.Mascot)& is.na(Peptidos_PP4$Ions.Score.C4.Mascot)& is.na(Peptidos_PP4$Ions.Score.D4.Mascot)& is.na(Peptidos_PP4$Ions.Score.E4.Mascot)& is.na(Peptidos_PP4$Ions.Score.F4.Mascot),])

        test48 <- data.frame(Peptidos_PP4[ is.na(Peptidos_PP4$Ions.Score.A4.Mascot)& is.na(Peptidos_PP4$Ions.Score.B4.Mascot)& is.na(Peptidos_PP4$Ions.Score.C4.Mascot)& is.na(Peptidos_PP4$Ions.Score.D4.Mascot)& is.na(Peptidos_PP4$Ions.Score.E4.Mascot)& is.na(Peptidos_PP4$Ions.Score.F4.Mascot),])

        test49 <- data.frame(Peptidos_PP4[ is.na(Peptidos_PP4$Ions.Score.A4.Mascot)& is.na(Peptidos_PP4$Ions.Score.B4.Mascot)& is.na(Peptidos_PP4$Ions.Score.C4.Mascot)& is.na(Peptidos_PP4$Ions.Score.D4.Mascot)& is.na(Peptidos_PP4$Ions.Score.E4.Mascot)& is.na(Peptidos_PP4$Ions.Score.F4.Mascot),])

        test50 <- data.frame(Peptidos_PP4[ is.na(Peptidos_PP4$Ions.Score.A4.Mascot)& is.na(Peptidos_PP4$Ions.Score.B4.Mascot)& is.na(Peptidos_PP4$Ions.Score.C4.Mascot)& is.na(Peptidos_PP4$Ions.Score.D4.Mascot)& is.na(Peptidos_PP4$Ions.Score.E4.Mascot)& is.na(Peptidos_PP4$Ions.Score.F4.Mascot),])

        test51 <- data.frame(Peptidos_PP4[ is.na(Peptidos_PP4$Ions.Score.A4.Mascot)& is.na(Peptidos_PP4$Ions.Score.B4.Mascot)& is.na(Peptidos_PP4$Ions.Score.C4.Mascot)& is.na(Peptidos_PP4$Ions.Score.D4.Mascot)& is.na(Peptidos_PP4$Ions.Score.E4.Mascot)& is.na(Peptidos_PP4$Ions.Score.F4.Mascot),])

        test52 <- data.frame(Peptidos_PP4[ is.na(Peptidos_PP4$Ions.Score.A4.Mascot)& is.na(Peptidos_PP4$Ions.Score.B4.Mascot)& is.na(Peptidos_PP4$Ions.Score.C4.Mascot)& is.na(Peptidos_PP4$Ions.Score.D4.Mascot)& is.na(Peptidos_PP4$Ions.Score.E4.Mascot)& is.na(Peptidos_PP4$Ions.Score.F4.Mascot),])

        test53 <- data.frame(Peptidos_PP4[ is.na(Peptidos_PP4$Ions.Score.A4.Mascot)& is.na(Peptidos_PP4$Ions.Score.B4.Mascot)& is.na(Peptidos_PP4$Ions.Score.C4.Mascot)& is.na(Peptidos_PP4$Ions.Score.D4.Mascot)& is.na(Peptidos_PP4$Ions.Score.E4.Mascot)& is.na(Peptidos_PP4$Ions.Score.F4.Mascot),])

        test54 <- data.frame(Peptidos_PP4[ is.na(Peptidos_PP4$Ions.Score.A4.Mascot)& is.na(Peptidos_PP4$Ions.Score.B4.Mascot)& is.na(Peptidos_PP4$Ions.Score.C4.Mascot)& is.na(Peptidos_PP4$Ions.Score.D4.Mascot)& is.na(Peptidos_PP4$Ions.Score.E4.Mascot)& is.na(Peptidos_PP4$Ions.Score.F4.Mascot),])

        test55 <- data.frame(Peptidos_PP4[ is.na(Peptidos_PP4$Ions.Score.A4.Mascot)& is.na(Peptidos_PP4$Ions.Score.B4.Mascot)& is.na(Peptidos_PP4$Ions.Score.C4.Mascot)& is.na(Peptidos_PP4$Ions.Score.D4.Mascot)& is.na(Peptidos_PP4$Ions.Score.E4.Mascot)& is.na(Peptidos_PP4$Ions.Score.F4.Mascot),])

        test56 <- data.frame(Peptidos_PP4[ is.na(Peptidos_PP4$Ions.Score.A4.Mascot)& is.na(Peptidos_PP4$Ions.Score.B4.Mascot)& is.na(Peptidos_PP4$Ions.Score.C4.Mascot)& is.na(Peptidos_PP4$Ions.Score.D4.Mascot)& is.na(Peptidos_PP4$Ions.Score.E4.Mascot)& is.na(Peptidos_PP4$Ions.Score.F4.Mascot),])

        test57 <- data.frame(Peptidos_PP4[ is.na(Peptidos_PP4$Ions.Score.A4.Mascot)& is.na(Peptidos_PP4$Ions.Score.B4.Mascot)& is.na(Peptidos_PP4$Ions.Score.C4.Mascot)& is.na(Peptidos_PP4$Ions.Score.D4.Mascot)& is.na(Peptidos_PP4$Ions.Score.E4.Mascot)& is.na(Peptidos_PP4$Ions.Score.F4.Mascot),])

        test58 <- data.frame(Peptidos_PP4[ is.na(Peptidos_PP4$Ions.Score.A4.Mascot)& is.na(Peptidos_PP4$Ions.Score.B4.Mascot)& is.na(Peptidos_PP4$Ions.Score.C4.Mascot)& is.na(Peptidos_PP4$Ions.Score.D4.Mascot)& is.na(Peptidos_PP4$Ions.Score.E4.Mascot)& is.na(Peptidos_PP4$Ions.Score.F4.Mascot),])

        test59 <- data.frame(Peptidos_PP4[ is.na(Peptidos_PP4$Ions.Score.A4.Mascot)& is.na(Peptidos_PP4$Ions.Score.B4.Mascot)& is.na(Peptidos_PP4$Ions.Score.C4.Mascot)& is.na(Peptidos_PP4$Ions.Score.D4.Mascot)& is.na(Peptidos_PP4$Ions.Score.E4.Mascot)& is.na(Peptidos_PP4$Ions.Score.F4.Mascot),])

        test60 <- data.frame(Peptidos_PP4[ is.na(Peptidos_PP4$Ions.Score.A4.Mascot)& is.na(Peptidos_PP4$Ions.Score.B4.Mascot)& is.na(Peptidos_PP4$Ions.Score.C4.Mascot)& is.na(Peptidos_PP4$Ions.Score.D4.Mascot)& is.na(Peptidos_PP4$Ions.Score.E4.Mascot)& is.na(Peptidos_PP4$Ions.Score.F4.Mascot),])

        test61 <- data.frame(Peptidos_PP4[ is.na(Peptidos_PP4$Ions.Score.A4.Mascot)& is.na(Peptidos_PP4$Ions.Score.B4.Mascot)& is.na(Peptidos_PP4$Ions.Score.C4.Mascot)& is.na(Peptidos_PP4$Ions.Score.D4.Mascot)& is.na(Peptidos_PP4$Ions.Score.E4.Mascot)& is.na(Peptidos_PP4$Ions.Score.F4.Mascot),])

        test62 <- data.frame(Peptidos_PP4[ is.na(Peptidos_PP4$Ions.Score.A4.Mascot)& is.na(Peptidos_PP4$Ions.Score.B4.Mascot)& is.na(Peptidos_PP4$Ions.Score.C4.Mascot)& is.na(Peptidos_PP4$Ions.Score.D4.Mascot)& is.na(Peptidos_PP4$Ions.Score.E4.Mascot)& is.na(Peptidos_PP4$Ions.Score.F4.Mascot),])

        test63 <- data.frame(Peptidos_PP4[ is.na(Peptidos_PP4$Ions.Score.A4.Mascot)& is.na(Peptidos_PP4$Ions.Score.B4.Mascot)& is.na(Peptidos_PP4$Ions.Score.C4.Mascot)& is.na(Peptidos_PP4$Ions.Score.D4.Mascot)& is.na(Peptidos_PP4$Ions.Score.E4.Mascot)& !is.na(Peptidos_PP4$Ions.Score.F4.Mascot),])



      }else{
      	#Print ya no
      }
    

    

    x <- data.frame(Peptidos_PP4[!is.na(Peptidos_PP4$Ions.Score.A4.Mascot),] & Peptidos_PP4[!is.na(Peptidos_PP4$Ions.Score.B4.Mascot),] )
    
    ##
    x <- data.frame(Peptidos_PP4[!is.na(Peptidos_PP4$Ions.Score.A4.Mascot) & !is.na(Peptidos_PP4$Ions.Score.B4.Mascot),])
    
    
    #Filter with the same peptides
    test1 <- filter(df_final,  !is.na(df_final$Ions.Score.A4.Mascot) & !is.na(df_final$Ions.Score.B4.Mascot)& !is.na(df_final$`Peptides(95%)`))
    test1 <- with(test1,  test1[order(-test1$Ions.Score.A4.Mascot) , ])
    #Different peptides between 1º and the 2º. 
    test2 <- filter(df_final,  !is.na(df_final$Ions.Score.A4.Mascot) & !is.na(df_final$Ions.Score.B4.Mascot)& is.na(df_final$`Peptides(95%)`))
    test2 <- with(test2,  test2[order(-test2$Ions.Score.A4.Mascot) , ])
    #Different peptides between 2º and the 1º.
    test3 <- filter(df_final,  is.na(df_final$Ions.Score.A4.Mascot) & !is.na(df_final$Ions.Score.B4.Mascot)& !is.na(df_final$`Peptides(95%)`))
    test3 <- with(test3,  test3[order(-test3$Ions.Score.B4.Mascot) , ])
    #Different between 3º and 2º,1º
    test4 <- filter(df_final,  !is.na(df_final$Ions.Score.A4.Mascot) & is.na(df_final$Ions.Score.B4.Mascot)& !is.na(df_final$`Peptides(95%)`))
    test4 <- with(test4,  test4[order(-test4$`Peptides(95%)`) , ])
    
    test5 <- filter(df_final,  !is.na(df_final$Ions.Score.A4.Mascot) & is.na(df_final$Ions.Score.B4.Mascot)& is.na(df_final$`Peptides(95%)`))
    test5 <- with(test5,  test5[order(-test5$Ions.Score.A4.Mascot) , ])
    
    test6 <- filter(df_final,  is.na(df_final$Ions.Score.A4.Mascot) & !is.na(df_final$Ions.Score.B4.Mascot)& is.na(df_final$`Peptides(95%)`))
    test6 <- with(test6,  test6[order(-test6$Ions.Score.B4.Mascot) , ])
    
    test7 <- filter(df_final,  is.na(df_final$Ions.Score.A4.Mascot) &  is.na(df_final$Ions.Score.B4.Mascot)& !is.na(df_final$`Peptides(95%)`))
    test7 <- with(test7,  test7[order(-test7$`Peptides(95%)`) , ])
    
    
    
    ########################################
    
    
    
    
    x <- list(Proteins = data.frame(Proteins_PP6), Peptides = data.frame(Peptidos_PP6))
    #WriteXLS(x, paste(gsub("*?PeptideSummary.txt","",files_glob_peptides[i]), "Summary.xlsx", sep = "", na=""), names(x))
    write.xlsx(x, file = paste(gsub("*?PeptideSummary.txt","",files_glob_peptides[i]), "Summary.xlsx", sep = "", na=""))
    
    
    i <- i + 1
  }
}

