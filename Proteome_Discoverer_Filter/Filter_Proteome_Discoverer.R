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
    Peptidos_PP <- read_delim(files_glob_peptides[i],
                              "\t", escape_double = FALSE, trim_ws = TRUE, col_types = cols())
    
    
    Peptidos_PP2 <- data.frame("Master Protein Accesions"=Peptidos_PP$`Master Protein Accessions`, Sequence=Peptidos_PP$Sequence, Modifications=Peptidos_PP$Modifications, "Protein Groups"=subset(Peptidos_PP, select = c(5)), Proteins=subset(Peptidos_PP, select = c(6)),"PSMs"=subset(Peptidos_PP, select = c(7)), "Missed Clavages"=subset(Peptidos_PP, select = c(9)), "Theo.MH+[Da]"=Peptidos_PP$`Theo. MH+ [Da]` )
    
    #choose the columns that we want. 
    Proteins_PP2 <- subset( Proteins_PP, select = -c(1,2))

    Peptidos_PP2 <- subset( Peptidos_PP, select = -c(1,2,8))
    Peptidos_PP2 <- data.frame("Master Protein Accesion"=Peptidos_PP$`Master Protein Accessions`, Peptidos_PP2)
    
    #Remove rows with REVERSED as string
    ###### Añadir aqui nuevo código para que te quite las Low en ambos lados. 
    Proteins_PP3 <- Proteins_PP2[!grepl("Low", Proteins_PP2$Proteins_PP.Name),]
    Peptidos_PP3 <- Peptidos_PP2[!grepl("REVERSED", Peptidos_PP2$Peptidos_PP.Names),]
    
    #Proteins_PP2[- grep("^REVERSED ", Proteins_PP2$Name),]
    #Proteins_PP2[ grep("REVERSED ", Proteins_PP2$Name, invert = TRUE) , ]
    
    # #Remove rows with Unised value => 1,3
    # Proteins_PP4 <- subset(Proteins_PP3, Proteins_PP3$Proteins_PP.Unused > 1.3 ) 
    # Proteins_PP4 <- subset(Proteins_PP4, Proteins_PP4$Proteins_PP...Cov.95.. > 0) 
    # Proteins_PP4 <- subset(Proteins_PP4, Proteins_PP4$Proteins_PP..Peptides.95... > 0)
    # 
    # 
    # #Filter peptide and filters
    # Peptidos_PP4 <- subset(Peptidos_PP3, Peptidos_PP3$Peptidos_PP.Unused > 1.3)
    # Peptidos_PP4 <- subset(Peptidos_PP4, Peptidos_PP4$Peptidos_PP.Contrib > 0) 
    # Peptidos_PP4 <- subset(Peptidos_PP4, Peptidos_PP4$Peptidos_PP.Conf > 95) 
    # Peptidos_PP4 <- subset(Peptidos_PP4, Peptidos_PP4$Peptidos_PP...Cov.95.. > 0) 
    # 
    # #Round numbers
    # colnames(Proteins_PP4) <- c("N", "Unused", "%Cov(95)", "Peptides(95%)", "Accession","Name","Species")
    # colnames(Peptidos_PP4) <- c("N","Unused",	"%Cov(95)",	"Accessions",	"Names",	"Contrib",	"Conf",	"Sequence",	"Modifications",	"Cleavages",	"dMass",	"Prec MW",	"Prec m/z",	"Theor MW",	"Theor m/z",	"Theor z",	"Sc",	"Spectrum",	"Time")
    # 
    # round_df <- function(df, digits) {
    #   nums <- vapply(df, is.numeric, FUN.VALUE = logical(1))
    #   
    #   df[,nums] <- round(df[,nums], digits = digits)
    #   
    #   (df)
    # }
    # 
    # Proteins_PP5 <- round_df(Proteins_PP4, digits=2)
    # Peptidos_PP5 <- round_df(Peptidos_PP4, digits=2)
    # 
    
    
    ##############
    #comparar tablas test: 
    
    duplicates <- data.frame(Peptidos_PP5$Accession)
    duplicate2 <- aggregate(list(numdup=rep(1,nrow(duplicates))), duplicates, length)
    duplicate3 <- with(duplicate2,  duplicate2[order(duplicate2$Peptidos_PP5.Accession) , ])
    
    #Protein accesion sorted
    duplicate4 <- with(Proteins_PP5,  Proteins_PP5[order(Proteins_PP5$Accession) , ])
    #Replace in protein table
    
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
    df2 <- data.frame("Accesion2"=gsub("\\;.*","", duplicate3$Peptidos_PP5.Accession), stringsAsFactors = FALSE)
    
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
      print(i)
      duplicate4 <- duplicate4[!duplicate4$Accession == as.character(i), ]
    }
    
    duplicate4$`Peptides(95%)` <- duplicate3$numdup
    duplicate4$Accession <- gsub("\\;.*","", duplicate4$Accession)
    
    
    #Sort columns by peptides
    Proteins_PP6 <- duplicate4[ order(-duplicate4[,4], -duplicate4[,2]), ]
    Proteins_PP6$N <- 1:nrow(Proteins_PP6)
    
    Proteins_PP6$Species <- gsub('\\d+\\.?','', Proteins_PP6$Species)

    #Compare and create a new column with element assigned by comparation with N peptide table. 
    
    Peptidos_PP6 <- Peptidos_PP5[order(Peptidos_PP5[,1],-Peptidos_PP5[,6]),]
    Peptidos_PP6$Accessions <- gsub("\\;.*","", Peptidos_PP6$Accessions)
    
    for (i in 1:length(Peptidos_PP6$Accessions)) {
      #print(match(prueba1_2$Peptidos_PP6.Accessions[i], prueba1$Proteins_PP6.Accession, nomatch=NA))
      Peptidos_PP6$N[i] <- match(Peptidos_PP6$Accessions[i],Proteins_PP6$Accession, nomatch=NA)
      
    } 
    
    Peptidos_PP6 <- Peptidos_PP6[order(Peptidos_PP6[,1],-Peptidos_PP6[,6]),]
    
    #Export xml
    # WriteXLS(Proteins_PP6, ExcelFileName = "Proteins_summary.xls", SheetNames = NULL, perl = "perl",
    #          verbose = FALSE, Encoding = c("UTF-8", "latin1", "cp1252"),
    #          row.names = FALSE, col.names = TRUE,
    #          AdjWidth = FALSE, AutoFilter = FALSE, BoldHeaderRow = FALSE,
    #          na = "",
    #          FreezeRow = 0, FreezeCol = 0,
    #          envir = parent.frame())
    
    
    #Check if there is any different between N and Proteins
    
    test <- data.frame(Peptidos_PP6$N)
    colnames(test) <- c("N")
    
    for (z in 1:nrow(test)){

        if (z==1){
            test$New[z] <- NA
            test$New[z] <- 1
    }
        else if (test$N[z]==test$N[z-1]){
            test$New[z] <- test$New[z-1]
        
        }
        else if (test$N[z]!=test$N[z-1]){
            test$New[z] <- test$New[z-1]+1
        }
    }
    
    Peptidos_PP6$N <- test$New
    
  
    x <- list(Proteins = data.frame(Proteins_PP6), Peptides = data.frame(Peptidos_PP6))
    #WriteXLS(x, paste(gsub("*?PeptideSummary.txt","",files_glob_peptides[i]), "Summary.xlsx", sep = "", na=""), names(x))
    write.xlsx(x, file = paste(gsub("*?PeptideSummary.txt","",files_glob_peptides[i]), "Summary.xlsx", sep = "", na=""))
    
    

    i <- i + 1
  }
}

