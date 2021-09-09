library(matrixStats)
library(randomForest)
library(penalizedSVM)
library(rknn)
library(caret)
library(mailR)
library(xlsx)
library(reshape2)
library(reshape)
library(sqldf)
library(tidyverse)
library(reshape2)
library(reshape)
library(gdata)



isrunning<<-FALSE

GetROIs<-function(x, FSFunction="ALL", test=FALSE){
  tryCatch( {
    print("GetROIs Started")
    set.seed(1)
    #Load the input list genes
    GenesList<-as.data.frame(x)
    #Initialize a variable where I will save the results and extra information
    Result<-list()
    #I use the SelectedProbes variable to store the result of filtering the Probe table,
    #according to the genes from the input available in the Atlas.
    
    
    SelectedProbes<-HumanBrainAtlas$Donor1$Probe[HumanBrainAtlas$Donor1$Probe$gene_symbol %in% GenesList[,1],]
    `%ni%` <- Negate(`%in%`)
    
    #Based on the ProbeIDs stored in SelectedProbes, I create a new column in Donor 1 Microarray Table,
    #marking as "1" those rows where ProbeIDs match with the available in SelectedProbes, and putting "0"
    #in the rest of the rows. Then, I save this new columns with the ProbeID columnd in Results$ProbeIDTarget
    
    Result$ProbeIDTarget<-as.data.frame(
      rbind(cbind(HumanBrainAtlas$Donor1$Microarray$Raw[HumanBrainAtlas$Donor1$Microarray$Raw[,1] %in% SelectedProbes[,1],]$ProbeID, 1)
            ,(cbind(HumanBrainAtlas$Donor1$Microarray$Raw[HumanBrainAtlas$Donor1$Microarray$Raw[,1] %ni% SelectedProbes[,1],]$ProbeID,0))))
    
    colnames(Result$ProbeIDTarget)<-c('ProbeID','Target')
    
    #Save the list of found genes, and the original input genes list.
    Result$FoundGenes<-as.list(unique(SelectedProbes$gene_symbol))
    print(paste0("Number of Found Genes:",length(Result$FoundGenes)))
    if(length(Result$FoundGenes)>300 | length(Result$FoundGenes)<50){
      FSFunction<-"W"
      print("FAST Wilcoxon")
    }
    
    Result$OriginalList<-as.list(GenesList)
    
    #IF test argument is FALSE, then the methods are applied over all the donors.
    #IF it is TRUE, then they are applied just over the donor 1 (test mode ON)
    if(test==FALSE){
      ndonors=6
    }else{
      ndonors=1
    }
    
    
    for(i in 1:ndonors){
      #Use Donor variable for looping through the Donors Microarrays in the Allen Human Brain Atlas
      isrunning<<-TRUE
      Donor<-paste('Donor',i, sep="")
      
      print(Donor)
      
      
      #En base al conjunto ProbeIDTarget, donde hab?a marcado con "1" los ProbeID que correspond?an a los genes con
      #los que me pidieron trabajar y con "0" el resto, traslado dichos valores binarios al conjunto de variables predictoras
      #de cada donante, utilizando la funci?n "merge" a partir de la variable ProbeID.
      data<-merge(Result$ProbeIDTarget, HumanBrainAtlas[[Donor]]$Microarray$Standarized, by='ProbeID')
      target<-data$Target
      #Borro las variables ProbeID y Target, las cuales no forman parte de la selecci?n de variables
      data$ProbeID<-NULL
      data$Target<-NULL
      #Renombro apropiadamente cada una de las variables, seg?n el SampleAnnot de cada donante.
      colnames(data)<-c(HumanBrainAtlas[[Donor]]$SampleAnnot$well_id)
      
      #De acuerdo al valor del par?metro FSFunction que hayamos usado al llamar a la funci?n GetTOIs, se invocar?
      #a algunos de los 4 m?todos de selecci?n de variables disponibles.
      #Luego de obtener el resultado por parte del m?todo invocado, se procede a postprocesar los mismos, renombrando
      #debidamente las columnas, ordenando los datos de acuerdo a su ?ndice de importancia, y agregando a cada regi?n del cerebro
      #aquellos datos de los que dispone el atlas (usando la funci?n merge(), entre el resultado y la informaci?n disponible en el SampleAnnot de cada donante)
      
      ################### WilconxonFS ##############################################
      if(FSFunction=="W" || FSFunction=="ALL"  ){
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
        #write.csv(Result$WilcoxonFS[[Donor]], file = paste(resultsdir,"/WilcoxonFS.",Donor,".csv", sep = ""))
        #write.xlsx(Result$WilcoxonFS[[Donor]], file = paste(resultsdir,"/WilcoxonFS.xlsx", sep=""), sheetName=Donor, row.names=FALSE, append = TRUE)
        #save(Result, file = paste(resultsdir,"/WilcoxonFS.Rda", sep=""))
      }
      
      
      
      #En los siguiente 3 m?todos wrappers, realizo el mismo procedimiento:
      #1) Llamo a la funci?n que realiza propiamente la selecci?n de variable.
      #2) Al dataframe varimp le adjunto la informaci?n disponible en el dataframe SampleAnnot de cada donante (merge()).
      #3) A la variable que contiene el ?ndice de importancia, lo convierto en num?rico y lo ordeno de mayor a menor.
      #4) Guardo en un archivo csv el dataframe ***varimp, con la denominaci?n correspondiente.
      
      ################### RandomForestFS ##############################################
      
      if(FSFunction=='RF' || FSFunction=="ALL" ){
        print('RandomForestFS')
        Result$RandomForestFS[[Donor]]<-RandomForestFS(data,target)
        Result$RandomForestFS[[Donor]]$rfvarimp<-merge(Result$RandomForestFS[[Donor]]$rfvarimp, HumanBrainAtlas[[Donor]]$SampleAnnot, by='well_id')
        Result$RandomForestFS[[Donor]]$rfvarimp$rfimportance<-as.numeric(Result$RandomForestFS[[Donor]]$rfvarimp$rfimportance)
        Result$RandomForestFS[[Donor]]$rfvarimp<-Result$RandomForestFS[[Donor]]$rfvarimp[order(-Result$RandomForestFS[[Donor]]$rfvarimp$rfimportance),]
        #write.csv(Result$RandomForestFS[[Donor]]$rfvarimp, file = paste(resultsdir,"/RandomForestFS.",Donor,".csv", sep = "") )
        #write.xlsx(Result$RandomForestFS[[Donor]]$rfvarimp, file = paste(resultsdir,"/RandomForestFS.xlsx", sep=""), sheetName=Donor, row.names=FALSE, append = TRUE)
        #save(Result, file = paste(resultsdir,"/RFFS.Rda", sep=""))
      }
      ################### SVMFS ##############################################
      
      if(FSFunction=='SVM' || FSFunction=="ALL" ){
        print('SVMFS')
        Result$SVMFS[[Donor]]<-SVMFS(data,target)
        Result$SVMFS[[Donor]]$svmvarimp<-merge(Result$SVMFS[[Donor]]$svmvarimp, HumanBrainAtlas[[Donor]]$SampleAnnot, by='well_id')
        Result$SVMFS[[Donor]]$svmvarimp$svmimportance<-as.numeric(Result$SVMFS[[Donor]]$svmvarimp$svmimportance)
        Result$SVMFS[[Donor]]$svmvarimp<-Result$SVMFS[[Donor]]$svmvarimp[order(-Result$SVMFS[[Donor]]$svmvarimp$svmimportance),]
        #write.csv(Result$SVMFS[[Donor]]$svmvarimp,file = paste(resultsdir,"/SVMFS.",Donor,".csv", sep = "") )
        #write.xlsx(Result$SVMFS[[Donor]]$svmvarimp, file = paste(resultsdir,"/SVMFS.xlsx", sep=""), sheetName=Donor, row.names=FALSE, append = TRUE)
        #save(Result, file = paste(resultsdir,"/SVMFS.Rda", sep=""))
      }
      
      ################### rKNN ##############################################
      
      
      if(FSFunction=='KNN' || FSFunction=="ALL" ){
        print('rKNNFS')
        Result$rKNNFS[[Donor]]<-rKNNFS(data,target)
        Result$rKNNFS[[Donor]]$knnvarimp<-merge(Result$rKNNFS[[Donor]]$knnvarimp, HumanBrainAtlas[[Donor]]$SampleAnnot, by='well_id')
        Result$rKNNFS[[Donor]]$knnvarimp$knnimportance<-as.numeric(Result$rKNNFS[[Donor]]$knnvarimp$knnimportance)
        Result$rKNNFS[[Donor]]$knnvarimp<-Result$rKNNFS[[Donor]]$knnvarimp[order(-Result$rKNNFS[[Donor]]$knnvarimp$knnimportance),]
        #write.csv(Result$rKNNFS[[Donor]]$knnvarimp, file = paste(resultsdir,"/KNN.",Donor,".csv", sep = "") )
        #write.xlsx(Result$rKNNFS[[Donor]]$knnvarimp, file = paste(resultsdir,"/KNN.xlsx", sep=""), sheetName=Donor, row.names=FALSE, append = TRUE)
        #save(Result, file = paste(resultsdir,"/KNN.Rda", sep=""))
      }
    }
    
    Result<-summarizePGL(Result)
    
    
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
  #Uno los predictores con la variable target
  data<-cbind(x, y)
  #Renombro la ?ltima variable como "target"
  colnames(data)[ncolX+1] <- "target"
  #Separo en dos conjuntos de datos distintos las filas con distinto Target
  data.class0<-subset(data, target==0)
  data.class1<-subset(data, target==1)
  
  #Borro la variable Target de ambos conjuntos, para excluirla del an?lisis
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
  
  
  #Por cada variable que posee el conjunto de datos recibido, realizo la prueba de Wilcoxon,
  #entre las filas con target=0 y las que poseen target=1
  for(i in 1:ncoldata){
    cat("\r",paste('Area:',i))
    wilcoxtest[[i]]<-try(wilcox.test(data.class1[,i], data.class0[,i], conf.int = confinter, conf.level = 0.95))
    #Si en alguna de las variables, durante la prueba de Wilcoxon se genera un error, coloco un cero en los resultados
    #que se hubiesen generado, y continuo con la siguiente variable.
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
###################################     RandomForestFS     ##################################################################
####################################################################################################################


RandomForestFS<-function(x,y){
  #library(randomForest)
  print("RanbdomForestFS Started")
  
  #Las l?neas a continuaci?n realizan un undersamplig sobre los genes No seleccionados
  #Uno las predictoras con la variable target en un la variable "x".
  x<-cbind(x,y)
  #Separo en una variable aparte las filas de "x" correspondientes a los SNPs seleccionados (Target=1)
  x.class1<-x[x[,ncol(x)]==1,]
  #Calculo cu?ntos SNPs seleccionados tengo, y multiplico ese n?mero por un factor correspondiente
  #a la cantidad de SNPs no seleccionados con los que quiero trabajar en los modelos
  nclass1<-nrow(x.class1*5)
  #Separo en una variable aparte los SNPs no seleccionados.
  x.class0<-x[x[,ncol(x)]==0,]
  #fijo un seed, para generar un undersamplig reproducible.
  set.seed(1)
  #Realizo un undersamplig sobre la variable donde guard? los SNPs no seleccionados, del tama?o calculado en
  #base a la cantidad de SNPs seleccionados multiplicado por un factor.
  class0.sample=sample(nrow(x.class0),nclass1)
  x.class0<-x.class0[class0.sample,]
  #Uno las filas de los SNPs seleccionados con los no seleccionados provenientes del undersamplig.
  x<-rbind(x.class1, x.class0)
  #Separo nuevamente el target y las predictoras en variables separadas "x" e "y".
  y<-x[,ncol(x)]
  x[,ncol(x)]<-NULL #Borra la variable target, para dejar en "x" solo las predictoras
  
  
  #Convierto la variable target en factor, para que las funciones correspondientes creen modelos
  #clasificatorios y no de regresi?n.
  y<-as.factor(y)
  start.time <- Sys.time()
  print("RF full model started")
  #Creo un modelo RandomForest con todas las variables (ROIs) disponibles.
  fullmodel.rf<-randomForest(x, y, cv.fold = 3, mtry = 50, ntree = 500)
  print("RF full model done")
  print("rfcv started")
  #Realizo la selecci?n de variables con RandomForest y Backward feature selection.
  rfcv.model<-rfcv(x, y , cv.fold = 3, scale="log", step = 0.95)
  print("rfcv ended")
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  print(time.taken)
  result<-list()
  
  
  x<-as.data.frame(x[,rfcv.model$n.var ])
  result$models$selectedmodel.rf<-randomForest(x, y, cv.fold = 3, mtry = 50, ntree = 500)
  
  result$models$fullmodel<-fullmodel.rf
  result$models$rfcv<-rfcv.model
  
  result$rfvarimp<-cbind(rownames(result$models$selectedmodel.rf$importance),as.numeric(result$models$selectedmodel.rf$importance))
  colnames(result$rfvarimp)<-c('well_id','rfimportance')
  
  
  #fullmodel.rf.importance<-as.data.frame(fullmodel.rf$importance)
  #fullmodel.rf.importance$varname<-rownames(fullmodel.rf.importance)
  #fullmodel.rf.importance<-fullmodel.rf.importance[order(fullmodel.rf.importance$MeanDecreaseGini),]
  #result<-list()
  #result$top10features<-fullmodel.rf.importance[1:10,]
  #result$rfcvmodel<-rfcv.model
  #result$rffullmodel<-fullmodel.rf
  
  return(result)
}


####################################################################################################################
###################################     SVMFS     ####################################################################
####################################################################################################################


SVMFS<-function(x,y){
  
  #library('penalizedSVM')
  
  
  result<-list()
  
  x<-cbind(x,y)
  x.class1<-x[x[,ncol(x)]==1,]
  nclass1<-nrow(x.class1*5)
  x.class0<-x[x[,ncol(x)]==0,]
  set.seed(1)
  class0.sample=sample(nrow(x.class0),nclass1)
  x.class0<-x.class0[class0.sample,]
  x<-rbind(x.class1, x.class0)
  y<-x[,ncol(x)]
  x[,ncol(x)]<-NULL #Borra la variable target
  
  
  y[y==0]<-(-1)
  y<-as.factor(y)
  x<-as.matrix(x)
  
  start.time <- Sys.time()
  print("svm penalized started")
  invisible(capture.output(svm.penalized<-svmfs(x, y, fs.method = "scad+L2", cross.inner = 3, maxIter=10, maxevals=10, seed = 1, verbose = FALSE)))
  print("svm penalized ended")
  
  
  
  x<-as.data.frame(x[,svm.penalized$model$xind])
  fitControl <- trainControl(method = "repeatedcv",number = 3, repeats = 3)
  print("svm L2 started")
  invisible(capture.output(svm.L2 <- train(x=x, y=y, method = "svmLinearWeights2", trControl = fitControl, verbose = FALSE)))
  print("svm L2 ended")
  #Calculo la importancia de cada varia ble seleccionada en el paso anterior, y creo un dataframe con
  #con el formato y nombre de columnas adecuado, ordenado de mayor a menor por el ?ndice de importancia ("svmimportance"),
  #convirti?ndolo previamente en num?rico.
  #Guardo los modelos creados y el dataframe con la informaci?n sobre la importancia de cada variable en una lista
  #de resultados (result), para finalmente retornarla.
  result$svmvarimp<-varImp(svm.L2, scale = FALSE)
  result$svmvarimp<-as.data.frame(cbind(rownames(result$svmvarimp$importance),result$svmvarimp$importance$X1))
  colnames(result$svmvarimp)<-c('well_id','svmimportance')
  result$svmvarimp$svmimportance<-as.numeric(result$svmvarimp$svmimportance)
  result$svmvarimp<-result$svmvarimp[order(-result$svmvarimp$svmimportance),]
  result$models$svm.penalized<-svm.penalized
  result$models$svm.L2<-svm.L2
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  print(time.taken)
  return(result)
  
}


####################################################################################################################
###################################     KNNFS     ##################################################################
####################################################################################################################


rKNNFS<-function(x,y){
  #library('rknn')
  result<-list()
  x<-cbind(x,y)
  x.class1<-x[x[,ncol(x)]==1,]
  nclass1<-nrow(x.class1*5)
  x.class0<-x[x[,ncol(x)]==0,]
  set.seed(1)
  class0.sample=sample(nrow(x.class0),nclass1)
  x.class0<-x.class0[class0.sample,]
  x<-rbind(x.class1, x.class0)
  y<-x[,ncol(x)]
  x[,ncol(x)]<-NULL #Borra la variable target
  
  y<-as.matrix(y)
  y<-as.factor(y)
  
  
  rknn.r<-r(ncol(x), m = floor(sqrt(ncol(x))))
  start.time <- Sys.time()
  invisible(capture.output(result$models$rknn.Beg<-rknnBeg(x, y, k = 40, r = rknn.r, mtry=trunc(sqrt(ncol(traindev.x))), cluster=NULL, seed = NULL)))
  result$models$rknn.bestset<-bestset(result$models$rknn.Beg, criterion=c("mean_accuracy", "mean_support"))
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  time.taken
  
  x<-as.data.frame(x[,result$models$rknn.bestset])
  
  fitControl <- trainControl(method = "repeatedcv",number = 3, repeats = 3)
  print("KNN started")
  
  
  invisible(capture.output(result$models$knn <- train(x=x, y=y, method = "kknn", trControl = fitControl)))
  print("KNN ended")
  result$knnvarimp<-varImp(result$models$knn, scale = FALSE)
  
  result$knnvarimp<-as.data.frame(cbind(rownames(result$knnvarimp$importance),result$knnvarimp$importance$X1))
  colnames(result$knnvarimp)<-c('well_id','knnimportance')
  result$knnvarimp$knnimportance<-as.numeric(result$knnvarimp$knnimportance)
  result$knnvarimp<-result$knnvarimp[order(-result$knnvarimp$knnimportance),]
  return(result)
}


####################################################################################################################
####################################################################################################################
###################################    ROIs To Genes     ##################################################################
####################################################################################################################
####################################################################################################################

#' @export
GetGenes<-function(genes,ROIs){
  tryCatch( {
    set.seed(1)
    Result<-list()
    GenesList<-as.character(unlist(genes))
    
    ROIsMatrix<-GetROIsinDonors(as.data.frame(ROIs))
    ProbeIDs<-as.character(HumanBrainAtlas$Donor1$Probe[HumanBrainAtlas$Donor1$Probe$gene_symbol %in% GenesList,]$probe_id)
    
    #Bucle para aplicar los m?todos de selecci?n de variables en los 6 donantes
    
    for(i in 1:6){
      #La variable Donor v? cambiando su valor a medida que el bucle avanza sobre cada uno de los 6 donantes
      #adquiriendo los valores Donor1, Donor2,...,Donor6.
      Donor<-paste('Donor',i, sep="")
      #Imprime en pantalla una referencia para indicar cu?l es el Donante que est? procesando
      print(Donor)
      
      
      ROIsList<-(unlist(ROIsMatrix[,i]))
      
      
      ROIsList<-as.character(na.omit(ROIsList))
      
      
      data<-as.data.frame(HumanBrainAtlas[[Donor]]$Microarray$Standarized)
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
      
      
      
      
    } 
    return(Result)
  }, 
  error = function(cond) {
    return(cond)
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


#La funci?n RanfomForestFS recibe el dataframe de Microarray transpuesto y la correspondiente variable target,
#de cada uno de los donantes del atlas.

GetGenes.RandomForestFS<-function(data){
  
  library(AUCRF)
  
  
  start.time <- Sys.time()
  set.seed(1)
  #La funci?n AUCRF realiza la selecci?n de variables, con eliminaci?n de variables hacia atr?s,
  #con la intenci?n de optimizar el ?rea najo la curva ROC.
  
  fit <- AUCRF(target~., data=data, pdel=0.2)
  
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  time.taken
  #Obtengo el ranking de variables (genes), seg?n la m?trica Mean Decrease Gini
  fit$ranking<-as.data.frame(fit$ranking)
  rows<-rownames(fit$ranking, prefix = "row")
  fit$ranking$probe_id<-gsub("X","",rows)
  return(fit$ranking)
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


