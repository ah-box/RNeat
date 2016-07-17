#' Create a new pool of neural networks trained using the NEAT algorithm using formula notation
#'
#' @param formula specifies the dependent and explantory varibles using a formula
#' @param trainingData Is the data used to train the networks
#' @param nTrainingGenerations Number of generations / breeding cycles to use in the genetic mating
#' @param maxNumberOfNodes The maximum number of neural network nodes
#' @param speciesPopulation The maximum bumber of species
#' @return rneatneuralnet class with pool of genomes and training data
#' @example R/Examples/squareRoot.R
#' @export
rneatneuralnet <- function(formula,trainingData,nTrainingGenerations=10,maxNumberOfNodes = 500, speciesPopulation = 200){
  result <- varify.variables(formula,trainingData,nTrainingGenerations,maxNumberOfNodes,speciesPopulation)

  config<-newConfigNEAT(length(result$model.list$variables),length(result$model.list$response),result$maxNumberOfNodes,result$speciesPopulation)
  neatSim <- newNEATFormulaSimulation(config)

  rneatsim <- list(inputs=result,neatConfig=config,simulation=neatSim)
  class(rneatsim) <- "rneatneuralnet"
  return (rneatneuralnetcontinuetraining(rneatsim,nTrainingGenerations))
}

#' Continues the training of the neural networks / runs more generations
#'
#' @param rneatsim Is the class rneatneuralnet created using rneatneuralnet function
#' @param nTrainingGenerations Number of new generations to train
#' @example R/Examples/squareRoot.R
#' @return rneatneuralnet class with pool of genomes and training data
#' @export
rneatneuralnetcontinuetraining <- function(rneatsim,nTrainingGenerations){
  assertTrueFunc(is(rneatsim,"rneatneuralnet"),"rneatsim must be a of class rneatneuralnet")
  for(i in seq(1,nTrainingGenerations)){
    rneatsim <- NEATFormulaSimulation.RunSingleGeneration(rneatsim)
  }
  return(rneatsim)
}

varify.variables <- function (formula,trainingData,nTrainingGenerations,maxNumberOfNodes,speciesPopulation){
  if (is.null(trainingData))
    stop("'trainingData' is missing", call. = FALSE)
  if (is.null(formula))
    stop("'formula' is missing", call. = FALSE)
  if(is.null(nTrainingGenerations) || !is.numeric(nTrainingGenerations)){
    nTrainingGenerations <- 10
  }
  if(is.null(maxNumberOfNodes) || !is.numeric(maxNumberOfNodes)){
    maxNumberOfNodes <- 500
  }
  if(is.null(speciesPopulation) || !is.numeric(speciesPopulation)){
    speciesPopulation <- 200
  }

  trainingData <- as.data.frame(trainingData)
  formula <- as.formula(formula)
  model.vars <- attr(terms(formula), "term.labels")
  formula.reverse <- formula
  formula.reverse[[3]] <- formula[[2]]
  model.resp <- attr(terms(formula.reverse), "term.labels")
  model.list <- list(response = model.resp, variables = model.vars)

  formulaNames <- c(model.list$response,model.list$variables)
  trainingDataInNames <- colnames(trainingData)
  missingtrainingData <- formulaNames[!formulaNames %in% trainingDataInNames]
  if(length(missingtrainingData)!=0){
    stop(paste("The formula requests variables '",toString(missingtrainingData),"' not present in the input trainingData", call. = FALSE))
  }
  return(list(trainingData = trainingData, formula = formula,model.list=model.list,nTrainingGenerations=nTrainingGenerations,maxNumberOfNodes=maxNumberOfNodes,speciesPopulation=speciesPopulation))
}


#' Takes a pool of trained neural networks, selects the most fit and uses it to predict the depend variable from the input data
#'
#' @param rneatsim Is the class rneatneuralnet created using rneatneuralnet function
#' @param data is input data
#' @example R/Examples/squareRoot.R
#' @return data.frame class with predicted dependent variable
#' @export
compute <- function(rneatsim,data){
  result <- varify.compute.variables(rneatsim,data)
  genome <- findMostFitGenome(result$rneatsim$simulation)
  genome <- generateNetwork(genome,rneatsim$simulation$Config)

  nnetInputData <- as.data.frame(result$data[,result$rneatsim$inputs$model.list$variables])
  colnames(nnetInputData) <- result$rneatsim$inputs$model.list$variables
  nnetOutput <- t(as.data.frame(apply(nnetInputData,1,function(x) { evaluateNetwork(genome$Network,x,result$rneatsim$simulation$Config) })))
  colnames(nnetOutput) <- paste(result$rneatsim$inputs$model.list$response,"Pred",sep="")
  res <- cbind(nnetInputData,nnetOutput)
  rownames(res) <- rownames(result$data)
  return(res)
}

varify.compute.variables <- function(rneatsim,data){
  assertTrueFunc(is(rneatsim,"rneatneuralnet"),"rneatsim must be a of class rneatneuralnet")

  formulaNames <- c(rneatsim$model.list$response,rneatsim$model.list$variables)
  dataInNames <- colnames(data)
  missingDataNames <- formulaNames[!formulaNames %in% dataInNames]
  if(length(missingDataNames)!=0){
    stop(paste("The formula requests variables '",toString(missingDataNames),"' not present in the input data", call. = FALSE))
  }
  return(list(rneatsim=rneatsim,data=as.data.frame(data)))
}

findMostFitGenome <- function(simulation){
  maxFitness <- -Inf
  maxFitnessGenome <- NA
  for(i in seq(1,length(simulation$Pool$species))){
    if(length(simulation$Pool$species[[i]]$genomes) > 0){
      for(j in seq(1,length(simulation$Pool$species[[i]]$genomes))){
         fitness <- simulation$Pool$species[[i]]$genomes[[j]]$Fitness
         if(fitness > maxFitness){
           maxFitness <- fitness
           maxFitnessGenome <-simulation$Pool$species[[i]]$genomes[[j]]
         }
      }
    }
  }
  return(maxFitnessGenome)
}

newNEATFormulaSimulation <- function(neatConfig){
  performanceTracker <-  data.frame(generation=numeric(),minFitness=numeric(),maxFitness=numeric(),meanFitness=numeric(),medianFitness=numeric(),stringsAsFactors=FALSE)
  res <- list(Config=neatConfig,Pool=initialisePool(neatConfig),PerformanceTracker=performanceTracker)
  class(res) <- "NEATFormulaSimulation"
  return(res)
}

NEATFormulaSimulation.RunSingleGeneration <- function(rneatsim){
  assertTrueFunc(is(rneatsim,"rneatneuralnet"),"rneatsim must be a of class rneatneuralnet")
  counter <- 1
  nTot <- calcTotalNumOfGenomes(rneatsim$simulation)
  for(i in seq(1,length(rneatsim$simulation$Pool$species))){
    for(j in seq(1,length(rneatsim$simulation$Pool$species[[i]]$genomes))){
      rneatsim$simulation <- simulationFormulaRunner(rneatsim$simulation,i,j,100*counter/nTot,rneatsim$inputs$trainingData,rneatsim$inputs$model.list)
      counter <- counter + 1
    }
  }
  rneatsim$simulation$Pool <- newGeneration(rneatsim$simulation$Pool,rneatsim$simulation$Config)
  #print(paste("MaxFitness:",simulation$Pool$maxFitness))

  rneatsim$simulation$PerformanceTracker[rneatsim$simulation$Pool$generation,] <- c(rneatsim$simulation$Pool$generation,rneatsim$simulation$Pool$minFitness,rneatsim$simulation$Pool$maxFitness,rneatsim$simulation$Pool$meanFitness,rneatsim$simulation$Pool$medianFitness)
  print(rneatsim$simulation$PerformanceTracker[rneatsim$simulation$Pool$generation,])
  return (rneatsim)
}


simulationFormulaRunner <- function(simulation,speciesNum,genomeNum,pctSimulated,data,model.list){
  i<-speciesNum
  j <-genomeNum
  if(length(simulation$Pool$species[[i]]$genomes[[j]]$ConnectionGenes)>0){

    simulation$Pool$species[[i]]$genomes[[j]] <- generateNetwork(simulation$Pool$species[[i]]$genomes[[j]],simulation$Config)
    simulation$Pool$species[[i]]$genomes[[j]]$Fitness <- 0

    nnetInputData <- as.data.frame(data[,model.list$variables])
    nnetTargetOutputData <- as.data.frame(data[,model.list$response])
    #print(model.list$variables)
    #print(nnetInputData)
    nnetOutput <- t(as.data.frame(apply(nnetInputData,1,function(x) { evaluateNetwork(simulation$Pool$species[[i]]$genomes[[j]]$Network,x,simulation$Config) })))
    outsum <- cbind(nnetTargetOutputData,nnetOutput)
    #print(outsum)
    mse <- sum((nnetTargetOutputData - nnetOutput)^2)
    updatedFitness<-1/mse
    #print(paste("mse",mse,"fitness",updatedFitness))

    simulation$Pool$species[[i]]$genomes[[j]]$Fitness <- updatedFitness
    print(paste(round(pctSimulated,2),"% Finished simulation of species",i,"/",length(simulation$Pool$species),"genome",j,"/",length(simulation$Pool$species[[i]]$genomes),"with fitness",updatedFitness,"mse",mse))
    #print(paste(round(pctSimulated,2),"% Finished simulation of species",i,"genome",j,"with fitness",updatedFitness))
  } else {
    print(paste("Skipped simulation of species",i,"genome",j,"due to no connection genes"))
  }

  return(simulation)
}

#' @export
plot.rneatneuralnet <- function(data){
  plotPerformanceTracker(data$simulation$PerformanceTracker)
}

#' @export
plot.NEATFormulaSimulation <- function(data){
  genome <- findMostFitGenome(data)
  drawGenotypeNEAT(genome,data$Config)
  drawPhenotypeNEAT(genome,data$Config)
}
