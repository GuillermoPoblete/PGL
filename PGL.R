
GetROIs<-function(x, test=FALSE){
  
  #Install the required packages IF they are not already installed.
  packages = c("matrixStats", "randomForest", "penalizedSVM", "rknn", "caret", "xlsx", "reshape2", "sqldf", "tidyverse","gdata")
  package.check <- lapply(
    packages,
    FUN = function(x) {
      if (!require(x, character.only = TRUE)) {
        install.packages(x, dependencies = TRUE)
        library(x, character.only = TRUE)
      }
    }
  )
  
  tryCatch( {
    
    
    FSFunction<-"W"
    print("FAST Wilcoxon")
    currentPGLPath<-paste0(dirname(rstudioapi::getSourceEditorContext()$path),'/')
    setwd(currentPGLPath)

    GenesList<-read.csv(x, header = FALSE)
    
    if(!file.exists(paste0(currentPGLPath, "HumanBrainAtlas.rda"))){
      print('Downloading Allen Human Brain Atlas....please wait....It will take a few minutes')
      url <-'https://filetransfer.io/data-package/UewGYXiP/download'
      destfile<-'HumanBrainAtlas.rda'
      download.file(url, destfile, quiet = TRUE, mode = "wb")
    }
    print('Loading Human Brain Atlas')
    load(file="HumanBrainAtlas.rda", .GlobalEnv)
    
    print("GetROIs Started")
    set.seed(1)
    #Load the input list genes
    
    #Initialize a variable where I will save the results and extra information
    Result<-list()
    #I use the SelectedProbes variable to store the result of filtering the Probe table,
    #according to the genes from the input available in the Atlas.
    
    
    SelectedProbes<-HumanBrainAtlas$Donor1$Probe[HumanBrainAtlas$Donor1$Probe$gene_symbol %in% GenesList[,1],]
    `%ni%` <- Negate(`%in%`)
    
    #Based on the ProbeIDs stored in SelectedProbes, I create a new column in Donor 1 Microarray Table,
    #tagging as "1" those rows where ProbeIDs match with the available in SelectedProbes, and "0"
    #in the rest of the rows. Then, I save this new columns with the ProbeID column in Results$ProbeIDTarget
    
    Result$ProbeIDTarget<-as.data.frame(
      rbind(cbind(HumanBrainAtlas$Donor1$Microarray$Raw[HumanBrainAtlas$Donor1$Microarray$Raw[,1] %in% SelectedProbes[,1],]$ProbeID, 1)
            ,(cbind(HumanBrainAtlas$Donor1$Microarray$Raw[HumanBrainAtlas$Donor1$Microarray$Raw[,1] %ni% SelectedProbes[,1],]$ProbeID,0))))
    
    colnames(Result$ProbeIDTarget)<-c('ProbeID','Target')
    
    #Save the list of found genes, and the original input genes list.
    Result$FoundGenes<-as.list(unique(SelectedProbes$gene_symbol))
    print(paste0("Number of Found Genes:",length(Result$FoundGenes)))
    
    Result$OriginalList<-as.list(GenesList)
    
    
    
    
    for(i in 1:6){
      #Use Donor variable for looping through the Donors Microarrays in the Allen Human Brain Atlas
      Donor<-paste('Donor',i, sep="")
      
      print(Donor)
      
      #Merge the target with each Microarray, then separate the target from the 
      
      data<-merge(Result$ProbeIDTarget, HumanBrainAtlas[[Donor]]$Microarray$Standarized, by='ProbeID')
      target<-data$Target
      #Borro las variables ProbeID y Target, las cuales no forman parte de la selecci?n de variables
      data$ProbeID<-NULL
      data$Target<-NULL
      #Renombro apropiadamente cada una de las variables, seg?n el SampleAnnot de cada donante.
      colnames(data)<-c(HumanBrainAtlas[[Donor]]$SampleAnnot$well_id)

      ################### WilconxonFS ##############################################
      
      Result$WilcoxonFS[[Donor]]<-WilcoxonFS(data,target)
      if(ncol(Result$WilcoxonFS[[Donor]])>0){
        colnames(Result$WilcoxonFS[[Donor]])[1]<-'well_id'
        colnames(Result$WilcoxonFS[[Donor]])[2]<-'p-value'
        colnames(Result$WilcoxonFS[[Donor]])[3]<-'corrected.p-value'
        colnames(Result$WilcoxonFS[[Donor]])[4]<-'conf.int.inf'
        colnames(Result$WilcoxonFS[[Donor]])[5]<-'conf.int.sup'
        colnames(Result$WilcoxonFS[[Donor]])[6]<-'difference'
        colnames(Result$WilcoxonFS[[Donor]])[7]<-'W'
        Result$WilcoxonFS[[Donor]]<-merge(  Result$WilcoxonFS[[Donor]], HumanBrainAtlas[[Donor]]$SampleAnnot, by='well_id')
        Result$WilcoxonFS[[Donor]]$Donor<-i
      }
      
    }
    
    Result<-summarizePGL(Result)
    
    result.date<-gsub(':','.', paste0(x,' PGL GetROIs result ',Sys.time(), '.csv'))
    write.csv(Result$summarized$Merged$final2, file=result.date, row.names = FALSE)
    
    return(Result)
  }
  , error = function(cond) {
    return(cond)
    
  })
  
}



####################################################################################################################
###################################     WilcoxonFS     #############################################################
####################################################################################################################



WilcoxonFS<-function(x,y){
  #library('matrixStats')
  ncolX=ncol(x)
  #Join the Microarray data with the target feature
  data<-cbind(x, y)
  #Rename the target feature
  colnames(data)[ncolX+1] <- "target"
  #Split the Microarray rows based on target
  data.class0<-subset(data, target==0)
  data.class1<-subset(data, target==1)
  
  #Delete the target feature to exclude it from the analysis
  data.class0$Target<-NULL
  data.class1$Target<-NULL
  
  #calculo la cantidad de variables predictoras, para inicializar las variables donde guardar? los resultados
  #y luego procedo a inicializarlas
  ncoldata<-ncol(data.class0)
  wilcoxtest<-list()
  wilcoxtest.pvalue<-1:ncoldata
  wilcoxtest.confint1<-1:ncoldata
  wilcoxtest.confint2<-1:ncoldata
  wilcoxtest.statistic<-1:ncoldata
  wilcoxtest.difference<-1:ncoldata
  
  
  confinter<-FALSE
  
  #Perform the Wilcoxon test on each Microarray column, comparing the rows tagged as "0" vs
  #those tagged as "1"
  
  for(i in 1:ncoldata){
    cat("\r",paste('Area:',i))
    wilcoxtest[[i]]<-try(wilcox.test(data.class1[,i], data.class0[,i], conf.int = confinter, conf.level = 0.95))
    #IF in some test there is an error, I set the result's values to "0"
    
    if(class(wilcoxtest[[i]]) == "try-error" ){
      print("ERROR")
      wilcoxtest.pvalue[i]<-0
      wilcoxtest.confint1[i]<-0
      wilcoxtest.confint2[i]<-0
      wilcoxtest.difference[i]<-0
      wilcoxtest.statistic[i]<-0
      
    }
    #Si se llam? a la funci?n WilcoxonFS con el par?metro confinter=FALSE, entonces solo guardo el p-value y coloco
    #cero en el resto de los datos.
    if(confinter==FALSE){
      wilcoxtest.pvalue[i]<-wilcoxtest[[i]]$p.value
      wilcoxtest.confint1[i]<-0
      wilcoxtest.confint2[i]<-0
      wilcoxtest.difference[i]<-0
      wilcoxtest.statistic[i]<-0
    }else{
      #Si confinter=TRUE, guardo todos los datos, incluyendo los relacionados a los intervalos de confianza.
      wilcoxtest.pvalue[i]<-wilcoxtest[[i]]$p.value
      wilcoxtest.confint1[i]<-wilcoxtest[[i]]$conf.int[1]
      wilcoxtest.confint2[i]<-wilcoxtest[[i]]$conf.int[2]
      wilcoxtest.difference[i]<-wilcoxtest[[i]]$estimate
      wilcoxtest.statistic[i]<-wilcoxtest[[i]]$statistic
    }
    
  }
  
  #Aplico la correcci?n de Bonferroni a cada uno de los p-values obtenidos
  wilcoxtest.bonferroni<- p.adjust(wilcoxtest.pvalue, "bonferroni")
  
  #Inicializo el dataframe features.selected, donde ordenar? todos los datos obtenidos sobre las variables seleccionadas.
  features.selected<-data.frame()
  
  #Recupero los nombres de las columnas del conjunto de datos original, los cuales servir?n ahora como los datos de
  #la columna que contiene los nombres de las variables seleccionadas.
  Header<-colnames(x)
  
  #Reviso el p-value corregido de cada una de los test realizados sobre cada variable. Si es menor o igual que 0.05,
  #entonces lo coloco entre las variables seleccionadas.
  for(i in 1:ncoldata){
    if(wilcoxtest.bonferroni[i]<=0.05){
      features.selected<-rbind(features.selected,cbind(as.character(Header[i]),wilcoxtest.pvalue[i], wilcoxtest.bonferroni[i], wilcoxtest.confint1[i],wilcoxtest.confint2[i],wilcoxtest.difference[i], wilcoxtest.statistic[i]))
      
      features.selected[,1]<-as.character(features.selected[,1])
    }
  }
  
  return(features.selected)
}





####################################################################################################################
####################################################################################################################
###################################    ROIs To Genes     ##################################################################
####################################################################################################################
####################################################################################################################

GetGenes<-function(genes,ROIs){
  packages = c("gdata","matrixStats", "randomForest", "penalizedSVM", "rknn", "caret", "xlsx", "reshape2", "sqldf", "tidyverse","gdata")
  
  ## Now load or install&load all
  package.check <- lapply(
    packages,
    FUN = function(x) {
      if (!require(x, character.only = TRUE)) {
        install.packages(x, dependencies = TRUE)
        library(x, character.only = TRUE)
      }
    }
  )
  tryCatch( {
    currentPGLPath<-paste0(dirname(rstudioapi::getSourceEditorContext()$path),'/')
    setwd(currentPGLPath)
    
    if(!file.exists(paste0(currentPGLPath, "HumanBrainAtlas.rda"))){
      print('Downloading Allen Human Brain Atlas....please wait....It will take a few minutes')
      url <-'https://filetransfer.io/data-package/UewGYXiP/download'
      destfile<-'HumanBrainAtlas.rda'
      download.file(url, destfile, quiet = TRUE, mode = "wb")
    }
    print('Loading Human Brain Atlas')
    load(file="HumanBrainAtlas.rda", .GlobalEnv )
    print(length(HumanBrainAtlas))
    
    
    
    
    set.seed(1)
    ROIs<-read.csv(file=ROIs, header = FALSE)
    genes<-read.csv(file=genes, header=FALSE)
    Result<-list()
    GenesList<-as.character(unlist(genes))
    
    ROIsMatrix<-GetROIsinDonors(as.data.frame(ROIs))
    #Obtengo los ProbeIDs de los genes con los que trabajar?
    
    
    
    ProbeIDs<-as.character(HumanBrainAtlas$Donor1$Probe[HumanBrainAtlas$Donor1$Probe$gene_symbol %in% GenesList,]$probe_id)
    
    #Bucle para aplicar los m?todos de selecci?n de variables en los 6 donantes
    filename<-paste(Sys.time(),'.GetGenes.Results.xlsx', sep="")
    filename<-gsub("-",".",filename)
    filename<-gsub(":",".",filename)
    filename<-gsub(" ",".",filename)
    for(i in 1:6){
      
      Donor<-paste('Donor',i, sep="")
      
      
      
      print(Donor)
      
      
      ROIsList<-(unlist(ROIsMatrix[,i]))
      
      
      ROIsList<-as.character(na.omit(ROIsList))
      
      
      Microarray.T<-HumanBrainAtlas[[Donor]]$Microarray$Standarized
      Microarray.T.colnames<-Microarray.T$ProbeID
      Microarray.T$ProbeID<-NULL
      Microarray.T<-t(Microarray.T)
      colnames(Microarray.T)<-Microarray.T.colnames
      data<-as.data.frame(Microarray.T)
      #Creo la variable target, ubicando un "1" en los casos donde los wellids coinciden con aquellos
      #pasados a la funci?n GetGenes, y un 0 en el resto.
      
      
      data$target<-0
      
      
      data[ROIsList,]$target<-1
      data<-as.data.frame(data)
      target<-as.factor(data$target)
      colnames(data)<-HumanBrainAtlas$Donor1$Probe$probe_id
      
      
      #Conservo solo aquellas columnas donde
      data<-data[,c(ProbeIDs)]
      data$target<-target
      data<-lapply(data, as.numeric)
      
      
      data$target<-as.factor(data$target)
      data<-as.data.frame(data)
      
      
      Result$RandomForestFS[[Donor]]<-GetGenes.RandomForestFS(data)
      
      
      #La funci?n RandomForestFS me entrega el ranking de genes, indic?ndome solo el ProbeID, por lo que debo
      #buscar el resto de los datos en el dataframe Probe de cada donante.
      Result$RandomForestFS[[Donor]]<-merge(Result$RandomForestFS[[Donor]], HumanBrainAtlas$Donor1$Probe, by='probe_id')
      
      Result$RandomForestFS[[Donor]]<-Result$RandomForestFS[[Donor]][order(-Result$RandomForestFS[[Donor]]$`fit$ranking`),]
      
      
      #Agrego al archivo final .xlsx la hoja con el ranking de genes, correspondiente a cada donante
      
      
      write.xlsx( Result$RandomForestFS[[Donor]], file=filename, sheetName=Donor, row.names=FALSE, append = TRUE)
      
    } 
    
    
  }, 
  error = function(cond) {
    message(paste("Ocurrio un error",cond, sep=" :"))
    })
}



GetROIsinDonors<-function(ROIs){
  toMatch <- as.vector(ROIs[,1])
  ROIs<-as.data.frame(unique(dplyr::filter(HumanBrainAtlas[["Donor1"]]$SampleAnnot, grepl(paste(toMatch,collapse="|"), structure_name)))$well_id)
  for(i in 2:6){
    Donor<-paste('Donor',i, sep="")
    ROIs<-cbindX(ROIs,as.data.frame(unique(dplyr::filter(HumanBrainAtlas[[Donor]]$SampleAnnot, grepl(paste(toMatch,collapse="|"), structure_name)))$well_id))
  }
  colnames(ROIs)<-c('Donor1','Donor2','Donor3','Donor4','Donor5','Donor6')
  return(ROIs)
}









###################### GetGenes Summary #########################




GetGenes.RandomForestFS<-function(data){
  
  library(AUCRF)
  
  
  start.time <- Sys.time()
  set.seed(1)
  
  fit <- AUCRF(target~., data=data, pdel=0.2)
  
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  time.taken
  
  fit$ranking<-as.data.frame(fit$ranking)
  rows<-rownames(fit$ranking, prefix = "row")
  fit$ranking$probe_id<-gsub("X","",rows)
  return(fit$ranking)
}






summarizePGL<-function(Result){
  print('Start Summarizing')
  Result$summarized<-list()
  ndonors<-c(1,2,3,4,5,6)
  methods<-vector()
  if(!is.null(Result$WilcoxonFS)){methods<-c(methods,'WilcoxonFS')}
  if(!is.null(Result$RandomForestFS)){methods<-c(methods,'RandomForestFS')}
  if(!is.null(Result$SVMFS)){methods<-c(methods,'SVMFS')}
  if(!is.null(Result$rKNNFS)){methods<-c(methods,'rKNNFS')}
  
  donorlist<-list()
  donorcount<-0
  for(i in 1:length(ndonors)){
    Donor<-paste('Donor',ndonors[i], sep="")
    Result$summarized[[Donor]]<-list()
    
    for(j in 1:length(methods)){
      donorcount<-donorcount+1
      method<-methods[j]
      temp<-Result[[method]][[Donor]]
      if(!is.null(temp)){
        if(method=='WilcoxonFS'){
          temp<-Result[[method]][[Donor]]
          colnames(temp)[3]<-'correctedpvalue'
          temp[,3]<-as.numeric(as.character(temp[,3]))
          colname<-'correctedpvalue'
          Result$summarized[[Donor]][[method]]<-sqldf(paste('SELECT structure_name, COUNT(structure_name), MIN(',colname,') FROM temp GROUP BY structure_name ORDER BY MIN(',colname,')', sep=""))
          
        }
        if(method=='RandomForestFS'){
          temp<-Result[[method]][[Donor]]$rfvarimp
          colname<-'rfimportance'
          temp[[colname]]<-as.numeric(as.character(temp[[colname]]))
          Result$summarized[[Donor]][[method]]<-sqldf(paste('SELECT structure_name, COUNT(structure_name), MAX(',colname,') FROM temp GROUP BY structure_name ORDER BY MAX(',colname,') DESC', sep=""))      
        }
        if(method=='SVMFS'){
          temp<-Result[[method]][[Donor]]$svmvarimp
          colname<-'svmimportance'
          temp[[colname]]<-as.numeric(as.character(temp[[colname]]))
          Result$summarized[[Donor]][[method]]<-sqldf(paste('SELECT structure_name, COUNT(structure_name), MAX(',colname,') FROM temp GROUP BY structure_name ORDER BY MAX(',colname,') DESC', sep=""))      
        }
        if(method=='rKNNFS'){
          temp<-Result[[method]][[Donor]]$knnvarimp
          colname<-'knnimportance'
          temp[[colname]]<-as.numeric(as.character(temp[[colname]]))
          Result$summarized[[Donor]][[method]]<-sqldf(paste('SELECT structure_name, COUNT(structure_name), MAX(',colname,') FROM temp GROUP BY structure_name ORDER BY MAX(',colname,') DESC', sep=""))      
        }
        colnames(Result$summarized[[Donor]][[method]])<-c('structurename','structurenamecount',paste('min',colname,sep =""))
        
        Result$summarized[[Donor]][[method]]<-Result$summarized[[Donor]][[method]][1:10,]
        Result$summarized[[Donor]][[method]]<-Result$summarized[[Donor]][[method]][complete.cases(Result$summarized[[Donor]][[method]]), ]
        Result$summarized[[Donor]][[method]]$ranking<-nrow(Result$summarized[[Donor]][[method]]):1/10
        donorlist[[donorcount]]<-as.data.frame(cbind(as.character(Result$summarized[[Donor]][[method]]$structurename), Result$summarized[[Donor]][[method]]$ranking))
        colnames(donorlist[[donorcount]])<-c('structurename', paste('Donor',i,method, "ranking", sep=""))
      }
      
      
      
    }
  }
  
  Result$summarized$Merged<-Reduce(function(x,y) merge(x = x, y = y, by = "structurename", all=TRUE, no.dups=FALSE),donorlist )     
  
  
  if (ncol(Result$summarized$Merged)==2){
    Result$summarized$Merged[,2]<-as.numeric(Result$summarized$Merged[,2])
  }else{
    Result$summarized$Merged[,2:ncol(Result$summarized$Merged)]<-Result$summarized$Merged[,2:ncol(Result$summarized$Merged)] %>%
      mutate_all(as.numeric)
  }
  
  
  Result$summarized$Merged[is.na(Result$summarized$Merged)]<-0
  
  Donor1Merged<-Result$summarized$Merged[grepl("Donor1", names(Result$summarized$Merged))]
  Donor2Merged<-Result$summarized$Merged[grepl("Donor2", names(Result$summarized$Merged))]
  Donor3Merged<-Result$summarized$Merged[grepl("Donor3", names(Result$summarized$Merged))]
  Donor4Merged<-Result$summarized$Merged[grepl("Donor4", names(Result$summarized$Merged))]
  Donor5Merged<-Result$summarized$Merged[grepl("Donor5", names(Result$summarized$Merged))]
  Donor6Merged<-Result$summarized$Merged[grepl("Donor6", names(Result$summarized$Merged))]
  
  Result$summarized$Merged$Mean1<-rowMeans(Donor1Merged)
  Result$summarized$Merged$Mean2<-rowMeans(Donor2Merged)
  Result$summarized$Merged$Mean3<-rowMeans(Donor3Merged)
  Result$summarized$Merged$Mean4<-rowMeans(Donor4Merged)
  Result$summarized$Merged$Mean5<-rowMeans(Donor5Merged)
  Result$summarized$Merged$Mean6<-rowMeans(Donor6Merged)
  
  meanbind<-as.data.frame(cbind(Result$summarized$Merged$Mean1, Result$summarized$Merged$Mean2, Result$summarized$Merged$Mean3, Result$summarized$Merged$Mean4, Result$summarized$Merged$Mean5, Result$summarized$Merged$Mean6))
  meanbind<-meanbind[colSums(!is.na(meanbind)) > 0]
  
  Result$summarized$Merged$final1<-cbind(as.character(Result$summarized$Merged$structurename), rowMeans(meanbind))
  Result$summarized$Merged$final1<-as.data.frame(Result$summarized$Merged$final1)
  colnames(Result$summarized$Merged$final1)<- c('structurename','finalranking')
  Result$summarized$Merged$final1$finalranking<-as.numeric(as.character(Result$summarized$Merged$final1$finalranking))
  Result$summarized$Merged$final1<-Result$summarized$Merged$final1[order(-Result$summarized$Merged$final1$finalranking),c(1,2)]
  #write.csv(Result$summarized$Merged$final1, file='Result.final1.csv')
  
  #############################################################################
  
  Result$summarized$Merged$Mean1<-rowMeans(Donor1Merged)*HumanBrainAtlas$Donor1$summary$BrainCompleteness
  Result$summarized$Merged$Mean2<-rowMeans(Donor2Merged)*HumanBrainAtlas$Donor2$summary$BrainCompleteness
  Result$summarized$Merged$Mean3<-rowMeans(Donor3Merged)*HumanBrainAtlas$Donor3$summary$BrainCompleteness
  Result$summarized$Merged$Mean4<-rowMeans(Donor4Merged)*HumanBrainAtlas$Donor4$summary$BrainCompleteness
  Result$summarized$Merged$Mean5<-rowMeans(Donor5Merged)*HumanBrainAtlas$Donor5$summary$BrainCompleteness
  Result$summarized$Merged$Mean6<-rowMeans(Donor6Merged)*HumanBrainAtlas$Donor6$summary$BrainCompleteness
  
  Result$summarized$Merged$final2<-cbind(as.character(Result$summarized$Merged$structurename), rowMeans(cbind(Result$summarized$Merged$Mean1, Result$summarized$Merged$Mean2, Result$summarized$Merged$Mean3, Result$summarized$Merged$Mean4, Result$summarized$Merged$Mean5, Result$summarized$Merged$Mean6)))
  Result$summarized$Merged$final2<-as.data.frame(Result$summarized$Merged$final2)
  colnames(Result$summarized$Merged$final2)<- c('structurename','finalranking')
  Result$summarized$Merged$final2$finalranking<-as.numeric(as.character(Result$summarized$Merged$final2$finalranking))
  Result$summarized$Merged$final2<-Result$summarized$Merged$final2[order(-Result$summarized$Merged$final2$finalranking),c(1,2)]
  
  #############################################################################
  
  Result$summarized$Merged$Mean1<-rowMeans(Donor1Merged)*HumanBrainAtlas$Donor1$summary$BrainCompleteness*HumanBrainAtlas$Donor1$summary$BrainStructureDensity
  Result$summarized$Merged$Mean2<-rowMeans(Donor2Merged)*HumanBrainAtlas$Donor2$summary$BrainCompleteness*HumanBrainAtlas$Donor2$summary$BrainStructureDensity
  Result$summarized$Merged$Mean3<-rowMeans(Donor3Merged)*HumanBrainAtlas$Donor3$summary$BrainCompleteness*HumanBrainAtlas$Donor3$summary$BrainStructureDensity
  Result$summarized$Merged$Mean4<-rowMeans(Donor4Merged)*HumanBrainAtlas$Donor4$summary$BrainCompleteness*HumanBrainAtlas$Donor4$summary$BrainStructureDensity
  Result$summarized$Merged$Mean5<-rowMeans(Donor5Merged)*HumanBrainAtlas$Donor5$summary$BrainCompleteness*HumanBrainAtlas$Donor5$summary$BrainStructureDensity
  Result$summarized$Merged$Mean6<-rowMeans(Donor6Merged)*HumanBrainAtlas$Donor6$summary$BrainCompleteness*HumanBrainAtlas$Donor6$summary$BrainStructureDensity
  
  Result$summarized$Merged$final3<-cbind(as.character(Result$summarized$Merged$structurename), rowMeans(cbind(Result$summarized$Merged$Mean1, Result$summarized$Merged$Mean2, Result$summarized$Merged$Mean3, Result$summarized$Merged$Mean4, Result$summarized$Merged$Mean5, Result$summarized$Merged$Mean6)))
  Result$summarized$Merged$final3<-as.data.frame(Result$summarized$Merged$final3)
  colnames(Result$summarized$Merged$final3)<- c('structurename','finalranking')
  Result$summarized$Merged$final3$finalranking<-as.numeric(as.character(Result$summarized$Merged$final3$finalranking))
  Result$summarized$Merged$final3<-Result$summarized$Merged$final3[order(-Result$summarized$Merged$final3$finalranking),c(1,2)]
  
  return(Result)
}




findgenes<-function(x){
  tryCatch({
    
    GenesList<-as.data.frame(x)
    View(GenesList)
    #Initialize a variable where I will save the results and extra information
    #I use the SelectedProbes variable to store the result of filtering the Probe table,
    #according to the genes from the input available in the Atlas.
    if(!nrow(GenesList)==0){
      SelectedProbes<-probe[probe$probe %in% GenesList[,1],]
      `%ni%` <- Negate(`%in%`)
      
      FoundGenes<-as.data.frame(unique(SelectedProbes))
      FoundGenes<-as.data.frame(FoundGenes[order(FoundGenes),])
      colnames(FoundGenes)[1]<-paste(nrow(FoundGenes),'Genes Found on the Atlas')
      }
      
    else{
      FoundGenes<-as.data.frame(c('No valid genes in the input list'))
      colnames(FoundGenes)[1]<-'0 Genes Found on the Atlas'
    }
    if(nrow(FoundGenes)==0){
      FoundGenes<-as.data.frame(c())
      colnames(FoundGenes)[1]<-'0 Genes Found on the Atlas'
    }
    return(FoundGenes)
  },error=function(cond){
  
    FoundGenes<-as.data.frame(c('Error'))
    return(FoundGenes)
  })
  
  
}

