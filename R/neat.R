
#' @importFrom methods is
#' @importFrom stats median runif
#' @importFrom animation saveVideo ani.options

config.video.phenotypedurationseconds = 1
config.video.performancedurationseconds = 1

#' Configuration for setting the number of system inputs/outputs, the max number of nodes and the total number of genomes
#'
#' @param numInputs The number of inputs to the neural network
#' @param numOutputs The number of outputs from the neural network
#' @param maxNumOfNodes The maximum number of neural network nodes
#' @param speciesPopulation The number of genomes to simulate
#' @return configNEAT class
#' @export
newConfigNEAT <- function(numInputs,numOutputs,maxNumOfNodes, speciesPopulation=200){
  assertTrueFunc(is(numInputs,"numeric"),"NumInputs must be a number")
  assertTrueFunc(is(numOutputs,"numeric"),"NumOutputs must be a number")
  assertTrueFunc(is(maxNumOfNodes,"numeric"),"MaxNumOfNodes must be a number")
  assertTrueFunc(is(speciesPopulation,"numeric"),"SpeciesPopulation must be a number")

  assertTrueFunc(numInputs>0,"NumInputs must be greater than 0")
  assertTrueFunc(numOutputs>0,"NumOutputs must be greater than 0")
  assertTrueFunc(maxNumOfNodes>0,"MaxNumOfNodes must be greater than 0")
  assertTrueFunc(speciesPopulation>0,"SpeciesPopulation must be greater than 0")

  res <- list(MutateConnectionChance=0.25,
              MutateLinkChance=2,
              MutateBiasChance=0.4,
              MutateNodeChance=0.5,
              MutateEnableChance=0.2,
              MutateDisableChance=0.4,
              MutateStepSize=0.1,
              MutationIncOrDecRate=0.05,
              PerturbChance=0.9,
              StaleSpecies=25,
              CrossoverChance=0.75,
              SpeciesPopulation=speciesPopulation,
              SpeciesDeltaDisjoint=2,
              SpeciesDeltaWeight=0.4,
              SpeciesDeltaThreshold=1,
              Inputs=numInputs,
              Outputs=numOutputs,
              MaxNodes=maxNumOfNodes
  )
  class(res) <- "configNEAT"
  return(res)
}


print.newConfigNEAT <- function(configNEAT){
  dumpItemFunc(configNEAT)
}

neatseq <- function(from,to){
  if(to==0){
    return (list())
  } else {
    return (seq(from,to))
  }
}

assertTrueFunc <- function(testExpression, messageIfNotTrue){
  if(!testExpression){
    cat(paste("Assertion Fail:",messageIfNotTrue,"\n"))
    stop()
  }
}

dumpItemFunc <- function(item){
  for(key in ls(item)){
    print(paste(key,item[key]))
  }
}

sigmoid <- function(x,mu=4.9){
  return(2/(1+exp(-mu*x))-1)
}


newInnovation <- function(){
  UseMethod("newInnovation")
}

pkg.env <- new.env()
pkg.env$innovation <- 0
newInnovation.default <- function(){
  #stop("New innovation function not implemented")
  pkg.env$innovation <- pkg.env$innovation + 1
  return (pkg.env$innovation)
}

newPool <- function(neatConfig){
  UseMethod("newPool")
}

newPool.default <- function(neatConfig){
  res <- list(species=list(),generation=0,innovation=neatConfig$Outputs,currentSpecies=1,currentGenome=1,currentFrame=0,maxFitness=0,minFitness=0,meanFitness=0,medianFitness=0)
  class(res) <- "pool"
  return(res)
}

print.pool <- function(pool){
  dumpItemFunc(pool)
}

newSpecies <- function(){
  UseMethod("newSpecies")
}

newSpecies.default <- function(){
  res <- list(topFitness=0,staleness=0,genomes=list(),averageFitness=0)
  class(res) <- "species"
  return(res)
}

print.species <- function(species){
  dumpItemFunc(species)
}

newGenome <- function(neatConfig){
  UseMethod("newGenome")
}

newGenome.default <- function(neatConfig){
  #for(item in nodeGenes) { assertTrueFunc(is(item,"nodegene"),"Node genes list must only contain nodegene class") }
  #for(item in connectionGenes) { assertTrueFunc(is(item,"connectiongene"),"Connection genes list must only contain connectiongene class") }
  res <- list(ConnectionGenes=list(),
              Fitness=0,
              AdjustedFitness=0,
              Network=list(),
              MaxNeuron=neatConfig$Inputs,
              GlobalRank=0,
              MutationRate=list(Connections=neatConfig$MutateConnectionChance,
                                Link=neatConfig$MutateLinkChance,
                                Bias=neatConfig$MutateBiasChance,
                                Node=neatConfig$MutateNodeChance,
                                Enable=neatConfig$MutateEnableChance,
                                Disable=neatConfig$MutateDisableChance,
                                Step=neatConfig$MutateStepSize))
  class(res) <- "genome"
  return(res)
}

print.genome <- function(genome){
  dumpItemFunc(genome)
}

basicgenome <- function(neatConfig){
  g <- newGenome(neatConfig)
  g$MaxNeuron <- neatConfig$Inputs
  return (mutateGenome(g,neatConfig))
}

newConnectiongene <- function(){
  UseMethod("newConnectiongene")
}

newConnectiongene.default <- function(){
  #assertTrueFunc(is(inNode,"nodegene"),"InNode must be a node gene")
  #assertTrueFunc(is(outNode,"nodegene"),"OutNode must be a node gene")
  #assertTrueFunc(is.numeric(weight),"Weight must be a numeric value")
  #assertTrueFunc(is(enabled,"logical"),"Enabled must be a boolean")
  #assertTrueFunc(is.numeric(innovation),"Innovation must be a numeric value")
  res <- list(InNode=NA,OutNode=NA,Weight=0,Enabled=T,Innovation=0)
  class(res) <- "connectiongene"
  return(res)
}

print.connectiongene <- function(connectiongene){
  dumpItemFunc(connectiongene)
}

newNeuron <- function(){
  UseMethod("newNeuron")
}

newNeuron.default <- function(){
  res <- list(Incoming=list(),Value=0)
  class(res) <- "neuron"
  return(res)
}

print.neuron <- function(neuron){
  dumpItemFunc(neuron)
}

newNetwork <- function(){
  UseMethod("newNetwork")
}

newNetwork.default <- function(){
  res <- list(Neurons=list())
  class(res) <- "network"
  return (res)
}

print.network <- function(network){
  dumpItemFunc(network)
}

#NEED TO BE EXCEPTIONALLY CAREFUL DUE TO PASS BY VALUE AND NOT REF LIKE LUA
generateNetwork <- function(genome,neatConfig){
  assertTrueFunc(is(genome,"genome"),"Genome must be a of class genome")
  network <- newNetwork()

  for(i in seq(1,neatConfig$Inputs)){
    network$Neurons[[i]] <- newNeuron()
  }

  for(o in seq(1,neatConfig$Outputs)){
    network$Neurons[[o+neatConfig$MaxNodes]] <- newNeuron()
  }
  #print(network)
  sortOutNodeIndex <- order(sapply(genome$ConnectionGenes,"[[","OutNode"),na.last=F) #NA Treatment?
  #print(sortOutNodeIndex)
  for(idx in sortOutNodeIndex){
    gene <- genome$ConnectionGenes[[idx]]
    #if(gene$Enabled){
    if(is.null(network$Neurons[[gene$OutNode]])){
      network$Neurons[[gene$OutNode]] <- newNeuron()
    }
    #print("Trying to set the gene")
    #print(gene)
    network$Neurons[[gene$OutNode]]$Incoming[[length(network$Neurons[[gene$OutNode]]$Incoming)+1]] <- gene

    if(is.null(network$Neurons[gene$InNode])){
      network$Neurons[[gene$InNode]] <- newNeuron()
    }
    #}
  }

  genome$Network <- network
  return (genome)

}


evaluateNetwork <- function(network,inputs,neatConfig){
  assertTrueFunc(is(network,"network"),"Network must be a of class network")


  #print("EvaluateNetwork")
  #print(network)
  if(length(inputs) != neatConfig$Inputs){
    stop("Number of inputs does not match the number in the config")
  }

  for(i in seq(1,neatConfig$Inputs)){
    network$Neurons[[i]]$Value <- inputs[[i]]
  }

  #print(network$Neurons)
  for(j in seq(1,length(network$Neurons))){
    totalSum <- 0
    if(!is.null(network$Neurons[[j]])){

      for(i in neatseq(1,length(network$Neurons[[j]]$Incoming))){ #Do i need a null filter here?
        incoming <- network$Neurons[[j]]$Incoming[[i]]
        if(!is.null(incoming) && incoming$Enabled){
          other <- network$Neurons[[incoming$InNode]]
          #  print("incoming")
          #  print(incoming)
          #  print("other")
          #  print(other)
          #  print(paste("Incoming weight is:",incoming$Weight,"Incoming value",other$Value))
          totalSum <- totalSum + incoming$Weight * other$Value
          #  print(paste("Total sum",totalSum))
        }
      }
      if(!is.null(network$Neurons[[j]]$Incoming) && length(network$Neurons[[j]]$Incoming) >0){
        network$Neurons[[j]]$Value <- sigmoid(totalSum)
      }
    }
  }

  #  print(network$Neurons)

  outputs <- list()
  for(i in seq(1,neatConfig$Outputs)){
    #if(network$Neurons[[i+config.MaxNodes]]$Value > 0){
    #  outputs[i] <- T
    #} else {
    #  outputs[i] <- F
    #}
    #print("Seting output")
    #print(network$Neurons[[i+config.MaxNodes]])
    outputs[i] <- network$Neurons[[i+neatConfig$MaxNodes]]$Value
  }

  return (outputs)
}

crossover <- function(g1,g2,neatConfig){
  assertTrueFunc(is(g1,"genome"),"g1 must be a of class genome")
  assertTrueFunc(is(g2,"genome"),"g2 must be a of class genome")

  # print("crossover g1")
  # print(g1)
  # print("crossover g2")
  # print(g2)
  #Make sure g1 is higher fitness genome
  if(g2$Fitness >  g1$Fitness){
    tmp <- g1
    g1 <- g2
    g2 <- tmp
  }

  child <- newGenome(neatConfig)

  g1Innovations <- unlist(lapply(g1$ConnectionGenes,function(x) { x$Innovation }))
  g2Innovations <- unlist(lapply(g2$ConnectionGenes,function(x) { x$Innovation }))


  innovations1 <- list()
  for(i in seq(1,length(g1$ConnectionGenes))){
    gene <- g1$ConnectionGenes[[i]]
    innovations1[[gene$Innovation]]<-gene
  }

  innovations2 <- list()
  for(i in seq(1,length(g2$ConnectionGenes))){
    gene <- g2$ConnectionGenes[[i]]
    innovations2[[gene$Innovation]]<-gene
  }


  for(innovation in unique(sort(c(g1Innovations,g2Innovations)))){
    gene1 <- innovations1[innovation]
    gene2 <- innovations2[innovation]
    #print(innovation)
    #print(gene1)
    #print(gene2)
    if(!is.null(gene1[[1]]) & !is.null(gene2[[1]])){

      if(g1$Fitness == g2$Fitness){
        #Randomly select either
        if(sample(1:2,1)==1){
          child$ConnectionGenes[[length(child$ConnectionGenes)+1]] <- gene1[[1]]
        } else {
          child$ConnectionGenes[[length(child$ConnectionGenes)+1]] <- gene2[[1]]
        }
      } else {
        #Select the most fit gene
        child$ConnectionGenes[[length(child$ConnectionGenes)+1]] <- gene1[[1]]
      }
    } else {
      if(!is.null(gene1[[1]])){
        child$ConnectionGenes[[length(child$ConnectionGenes)+1]] <- gene1[[1]]
      } else {

        child$ConnectionGenes[[length(child$ConnectionGenes)+1]] <- gene2[[1]]
      }

    }
  }

  child$MaxNeuron = max(g1$MaxNeuron,g2$MaxNeuron)
  child$MutationRate <- g1$MutationRate
  return(child)
}

randomNeuron <- function(genes,nonInput,neatConfig){
  assertTrueFunc(is(nonInput,"logical"),"nonInput must be a of class logical")
  neurons <- list()
  #print("genes randomNeuron")
  #print(genes)
  if(!nonInput){
    for(i in seq(1,neatConfig$Inputs)){
      neurons[i] <- T
    }
  }

  for(i in seq(1,neatConfig$Outputs)){
    neurons[neatConfig$MaxNodes+i] <- T
  }

  for(gene in genes){
    #print("RandomNeuron")
    #print(gene)
    #print(class(gene))
    if(!nonInput | gene$InNode > neatConfig$Inputs){
      neurons[gene$InNode] <- T
    }

    if(!nonInput | gene$OutNode > neatConfig$Inputs){
      neurons[gene$OutNode] <- T
    }
  }

  matches <- unlist(lapply(neurons,function(x) { return(!is.null(x))}))
  candidates <- seq(1,length(matches))[matches] #Index of non-null neurons
  idx <- sample(candidates,1)
  return (idx)

}

containsLink <- function(genes,link){
  for(gene in genes){
    if(gene$InNode == link$InNode & gene$OutNode == link$OutNode){
      return (T)
    }
  }

  return (F)
}

pointMutate <- function(genome,neatConfig){
  assertTrueFunc(is(genome,"genome"),"Genome must be a of class genome")

  stepSize <- genome$MutationRate[["Step"]]

  for(i in neatseq(1,length(genome$ConnectionGenes))){
    if(runif(1) < neatConfig$PerturbChance){
      genome$ConnectionGenes[[i]]$Weight <- genome$ConnectionGenes[[i]]$Weight+runif(1)*stepSize*2-stepSize
    } else {
      genome$ConnectionGenes[[i]]$Weight <- runif(1)*4-2
    }
  }
  return(genome)
}



linkMutate <- function(genome,forceBias,neatConfig){
  assertTrueFunc(is(genome,"genome"),"Genome must be a of class genome")
  assertTrueFunc(is(forceBias,"logical"),"forceBias must be a of class logical")

  neuron1 <- randomNeuron(genome$ConnectionGenes,F,neatConfig)
  neuron2 <- randomNeuron(genome$ConnectionGenes,T,neatConfig)

  newLink <- newConnectiongene()

  if(neuron1 <= neatConfig$Inputs & neuron2 <= neatConfig$Inputs){
    #Both neurons are input nodes
    return(genome)
  }

  if(neuron2 <= neatConfig$Inputs){
    #Swap output and input
    tmp <- neuron1
    neuron1<-neuron2
    neuron2<-tmp
  }



  newLink$InNode <- neuron1
  newLink$OutNode <- neuron2
  if(newLink$InNode > neatConfig$MaxNodes){
    return(genome) #avoids wiring an outnode to another node
  }

  if(forceBias){
    newLink$InNode <- neatConfig$Inputs
  }

  if(newLink$InNode == newLink$OutNode){
    return(genome) #Avoids creating a link that points to itself
  }

  if(containsLink(genome$ConnectionGenes,newLink)){
    return(genome)
  }

  newLink$Innovation <- newInnovation()
  newLink$Weight <- runif(1)*4-2
  genome$ConnectionGenes[[length(genome$ConnectionGenes)+1]] <- newLink
  return (genome)
}

#Disables a connection and replaces it with two that mimic the original
#Is adding a node
nodeMutate <- function(genome){
  assertTrueFunc(is(genome,"genome"),"Genome must be a of class genome")

  if(length(genome$ConnectionGenes)==0){
    return(genome)
  }

  genome$MaxNeuron <- genome$MaxNeuron + 1

  geneIndex <- sample(1:length(genome$ConnectionGenes),1)
  gene <- genome$ConnectionGenes[[geneIndex]]

  if(!gene$Enabled){
    return (genome)
  }

  genome$ConnectionGenes[[geneIndex]]$Enabled = F

  gene1 <- gene
  gene1$OutNode <- genome$MaxNeuron
  gene1$Weight <- 1
  gene1$Innovation <- newInnovation()
  gene1$Enabled <- T
  genome$ConnectionGenes[[length(genome$ConnectionGenes)+1]] <- gene1

  gene2 <- gene
  gene2$InNode <- genome$MaxNeuron
  gene2$Innovation <- newInnovation()
  gene2$Enabled <- T
  genome$ConnectionGenes[[length(genome$ConnectionGenes)+1]] <- gene2


  return (genome)
}

enableDisableMutate <- function(genome,state){
  assertTrueFunc(is(genome,"genome"),"Genome must be a of class genome")
  assertTrueFunc(is(state,"logical"),"State must be a of class logical")
  matches <- unlist(lapply(genome$ConnectionGenes,function(x) { return(x$Enabled!=state)}))
  candidates <- seq(1,length(matches))[matches] #Index of matching candidates

  if(length(candidates)==0){
    return (genome)
  }

  idx <- sample(candidates,1)
  genome$ConnectionGenes[[idx]]$Enabled <- state
  return (genome)
}

mutateGenome <- function(genome,neatConfig){
  assertTrueFunc(is(genome,"genome"),"Genome must be a of class genome")

  #Update the rate at which things mutate
  for(key in names(genome$MutationRate)){
    rate <- genome$MutationRate[[key]]
    if(sample(1:2,1)==1){
      genome$MutationRate[[key]] <- rate * (1-neatConfig$MutationIncOrDecRate)
    } else {
      genome$MutationRate[[key]] <- rate * (1/(1-neatConfig$MutationIncOrDecRate))
    }
  }

  if(runif(1) < genome$MutationRate[["Connections"]]){
    genome <- pointMutate(genome,neatConfig) #Done
  }

  p <- genome$MutationRate[["Link"]]
  while(p>0){
    if(runif(1)<p){
      genome <- linkMutate(genome,F,neatConfig) #done
    }
    p <- p - 1
  }

  p <- genome$MutationRate[["Bias"]]
  while(p>0){
    if(runif(1)<p){
      genome <- linkMutate(genome,T,neatConfig) #done
    }
    p <- p - 1
  }

  p <- genome$MutationRate[["Node"]]
  while(p>0){
    if(runif(1)<p){
      genome <- nodeMutate(genome) #done
    }
    p <- p - 1
  }

  p <- genome$MutationRate[["Enable"]]
  while(p>0){
    if(runif(1)<p){
      genome <- enableDisableMutate(genome,T) #done
    }
    p <- p - 1
  }

  p <- genome$MutationRate[["Disable"]]
  while(p>0){
    if(runif(1)<p){
      genome <- enableDisableMutate(genome,F) #done
    }
    p <- p - 1
  }
  return(genome)
}

#Finds the opposite of an intersection
outersect <- function(x, y) {
  sort(c(setdiff(x, y),
         setdiff(y, x)))
}

disjoint <-function(connectionGenesA,connectionGenesB){
  #print(connectionGenesA)
  #print(connectionGenesB)
  innovationsA <- unlist(lapply(connectionGenesA,function(x){x$Innovation}))
  innovationsB <- unlist(lapply(connectionGenesB,function(x){x$Innovation}))

  if(is.null(innovationsA)){
    innovationsA <- numeric()
  }
  if(is.null(innovationsB)){
    innovationsB <- numeric()
  }
  disjointInnovations <- outersect(innovationsA,innovationsB)
  nGenes <- max(length(connectionGenesA),length(connectionGenesB))
  if(nGenes == 0){
    return(0)
  }
  return (length(disjointInnovations)/nGenes)
}

weightsFunc <- function(connectionGenesA,connectionGenesB){

  innovations2 <- list()

  for(i in neatseq(1,length(connectionGenesB))){
    gene <- connectionGenesB[[i]]
    innovations2[[gene$Innovation]]<-gene
  }

  total <- 0
  coincident <- 0


  for(i in neatseq(1,length(connectionGenesA))){
    gene <- connectionGenesA[[i]]
    gene2 <- innovations2[gene$Innovation]

    if(!is.null(gene2[[1]])){
      total <- total + abs(gene$Weight+gene2[[1]]$Weight)
      coincident <- coincident + 1
    }

  }

  if(coincident==0){
    return(0)
  }
  return (total/coincident)
}

sameSpecies <-function(genomeA, genomeB,neatConfig){
  #print("sameSpeciesGenomeA")
  #print(genomeA)
  #print("sameSpeciesGenomeB")
  #print(genomeB)
  dd <- neatConfig$SpeciesDeltaDisjoint*disjoint(genomeA$ConnectionGenes, genomeB$ConnectionGenes)
  dw <- neatConfig$SpeciesDeltaWeight*weightsFunc(genomeA$ConnectionGenes, genomeB$ConnectionGenes)
  return ((dd + dw) < neatConfig$SpeciesDeltaThreshold)
}

rankGlobally <- function(pool){
  assertTrueFunc(is(pool,"pool"),"pool must be a of class pool")
  allFitnessScores <- sort(unlist(lapply(pool$species,function(y) { unlist(lapply(y$genomes,function(x){x$Fitness}))})))
  #print(paste("allfitnessscores",allFitnessScores))
  for(i in seq(1,length(pool$species))){

    for(j in neatseq(1,length(pool$species[[i]]$genomes))){
      pool$species[[i]]$genomes[[j]]$GlobalRank<-match(pool$species[[i]]$genomes[[j]]$Fitness,allFitnessScores)
    }
  }
  pool$maxFitness <- max(allFitnessScores)
  pool$minFitness <- min(allFitnessScores)
  pool$meanFitness <- mean(allFitnessScores)
  pool$medianFitness  <- median(allFitnessScores)
  return (pool)
}

calculateAverageFitness <- function(pool){
  assertTrueFunc(is(pool,"pool"),"pool must be a of class pool")
  for(i in seq(1,length(pool$species))){
    total <- sum(unlist(lapply(pool$species[[i]]$genomes,function(x){x$GlobalRank})))
    pool$species[[i]]$averageFitness <- (total/length(pool$species[[i]]$genomes))
    if(!is.numeric(pool$species[[i]]$averageFitness) | is.nan(pool$species[[i]]$averageFitness)){
      pool$species[[i]]$averageFitness <- 0
    }
  }
  return(pool)
}

totalAverageFitness <- function(pool){
  assertTrueFunc(is(pool,"pool"),"pool must be a of class pool")
  return (sum(unlist(lapply(pool$species,function(x){x$averageFitness}))))
}

cullSpecies <- function(pool,cutToOne){
  assertTrueFunc(is(pool,"pool"),"pool must be a of class pool")
  if(length(pool$species)== 0){
    print("cullSpecies for some reason there are no species in the pool???")
    return (pool)
  }
  for(i in seq(1,length(pool$species))){
    if(length(pool$species[[i]]$genomes)>0){
      fitness <- unlist(lapply(pool$species[[i]]$genomes,function(x){x$Fitness}))
      if(!is.numeric(fitness)){
        fitness <- 0
      }
      fitnessCutoff <- 0
      if(cutToOne){
        fitnessCutoff <- max(fitness)
      } else {
        fitnessCutoff <- median(fitness)
      }

      survivedGenomeIdx <- fitness >= fitnessCutoff
      pool$species[[i]]$genomes<-pool$species[[i]]$genomes[survivedGenomeIdx]
    }
  }
  return (pool)
}

breedChild <- function(species,neatConfig){
  #assertTrueFunc(is(species,"species"),"species must be a of class species")
  if(runif(1)<neatConfig$CrossoverChance){
    g1 <- species$genomes[[sample(1:length(species$genomes),1)]]
    g2 <- species$genomes[[sample(1:length(species$genomes),1)]]
    crossOverChild <- crossover(g1,g2,neatConfig)
    return (mutateGenome(crossOverChild,neatConfig))
  } else {
    child <- species$genomes[[sample(1:length(species$genomes),1)]]
    return (mutateGenome(child,neatConfig))
  }
}

removeStaleSpecies <- function(pool,neatConfig){
  assertTrueFunc(is(pool,"pool"),"pool must be a of class pool")
  for(i in seq(1,length(pool$species))){
    fitness <- unlist(lapply(pool$species[[i]]$genomes,function(x){x$Fitness}))
    if(!is.numeric(fitness)){
      fitness <- 0
    }
    if(max(fitness) > pool$species[[i]]$topFitness){
      pool$species[[i]]$topFitness <- max(fitness)
      pool$species[[i]]$staleness <- 0
    } else {
      pool$species[[i]]$staleness <- pool$species[[i]]$staleness + 1
    }
  }
  survivedIdx <- unlist(lapply(pool$species, function(x) { return (x$staleness < neatConfig$StaleSpecies ||x$topFitness >= pool$maxFitness) }))
  #print("survivedIdx")
  #print(survivedIdx)
  pool$species <- pool$species[survivedIdx]
  return (pool)
}

removeWeakSpecies <- function(pool,neatConfig){
  assertTrueFunc(is(pool,"pool"),"pool must be a of class pool")
  tavf <- totalAverageFitness(pool)
  survivedIdx <- unlist(lapply(pool$species,function(x){ return (floor(x$averageFitness/tavf*neatConfig$SpeciesPopulation) >= 1)}))
  pool$species <- pool$species[survivedIdx]
  return (pool)
}

addToSpecies <- function(pool,child,neatConfig){
  assertTrueFunc(is(pool,"pool"),"pool must be a of class pool")
  assertTrueFunc(is(child,"genome"),"child must be a of class genome")
  assertTrueFunc(is(neatConfig,"configNEAT"),"neatConfig must be a of class configNEAT")
  # print("child")
  # print(child)
  foundSpecies <- F

  for(i in neatseq(1,length(pool$species))){
    #print(pool$species[[i]])
    if(sameSpecies(child,pool$species[[i]]$genomes[[1]],neatConfig)){
      pool$species[[i]]$genomes[[length(pool$species[[i]]$genomes)+1]] <- child
      foundSpecies <- T
      break
    }
  }

  if(!foundSpecies){
    childspecies <- newSpecies()
    childspecies$genomes[[1]]<- child
    pool$species[[length(pool$species)+1]]<-childspecies
  }
  return (pool)
}

newGeneration <- function(pool,neatConfig){
  pool <- cullSpecies(pool,F)
  pool <- rankGlobally(pool)
  pool <- removeStaleSpecies(pool,neatConfig)
  pool <- rankGlobally(pool)
  pool <- calculateAverageFitness(pool)
  pool <- removeWeakSpecies(pool,neatConfig)
  tavf <- totalAverageFitness(pool)

  children <- list()

  for(i in neatseq(1,length(pool$species))){
    if(length(pool$species[[i]]$genomes) > 0){
      nBreed <- floor(pool$species[[i]]$averageFitness/tavf*neatConfig$SpeciesPopulation)-1
      if(!is.numeric(nBreed)){
        nBreed <- 0
      }
      #print(paste("nBreed",nBreed))
      #print(paste("averageFitness",pool$species[[i]]$averageFitness))
      #print(paste("tavf",tavf))
      if(nBreed >0){
        for(j in seq(1,nBreed)){
          children[[length(children)+1]]<-breedChild(pool$species[[i]],neatConfig)
        }
      }
    }
  }

  pool <- cullSpecies(pool,T)
  while((length(children)+length(pool$species))<neatConfig$SpeciesPopulation){
    children[[length(children)+1]] <- breedChild(pool$species[[sample(1:length(pool$species),1)]],neatConfig)
  }


  for(i in neatseq(1,length(children))){
    pool <- addToSpecies(pool,children[[i]],neatConfig)
  }


  pool$generation <- pool$generation+1

  return(pool)
}

initialisePool <- function(neatConfig){
  pool <- newPool(neatConfig)
  for(i in seq(1,neatConfig$SpeciesPopulation)){
    basic <- basicgenome(neatConfig)
    pool <- addToSpecies(pool,basic,neatConfig)
  }

  return(pool)
}


#' Create a new NEAT Simulation which contains the pool of species and genomes
#'
#' @param neatConfig Takes a NEATConfig class
#' @param processInitialStateFunc A function that specifies the initial condition (state) of the system
#' @param processUpdateStateFunc A function that takes the current system state and updates the state
#' @param processStateToNeuralInputFunc A function that takes the current state and converts it to neural net input
#' @param fitnessUpdateFunc A function that takes the current fitness level, old state and new state and returns the new fitness
#' @param terminationCheckFunc A function that returns TRUE if the simulation should be terminated
#' @param plotStateFunc A function that will plot what the current state is
#' @return NEATSimulation class with new pool of genomes
#' @example R/Examples/polebalance.R
#' @export
newNEATSimulation <- function(neatConfig,processInitialStateFunc,
                                      processUpdateStateFunc,
                                      processStateToNeuralInputFunc,
                                      fitnessUpdateFunc,
                                      terminationCheckFunc,
                                      plotStateFunc){
  performanceTracker <-  data.frame(generation=numeric(),minFitness=numeric(),maxFitness=numeric(),meanFitness=numeric(),medianFitness=numeric(),stringsAsFactors=FALSE)
  res <- list(Config=neatConfig,
              ProcessInitialStateFunc=processInitialStateFunc,
              ProcessUpdateStateFunc=processUpdateStateFunc,
              ProcessStateToNeuralInputFunc=processStateToNeuralInputFunc,
              FitnessUpdateFunc=fitnessUpdateFunc,
              TerminationCheckFunc=terminationCheckFunc,
              PlotStateFunc=plotStateFunc,
              Pool=initialisePool(neatConfig),PerformanceTracker=performanceTracker)
  class(res) <- "NEATSimulation"
  return(res)
}

#state <- ProcessInitialStateFunc()
#neuralNetInputs <- ProcessStateToNeuralInputFunc(currentState)
#updatedState <- ProcessUpdateStateFunc(currentState,neuralNetOutputs)
#newFitness <- FitnessUpdateFunc(oldState,updatedState,oldFitness)
#bool <- TerminationCheckFunc(oldState,updatedState,oldFitness,newFitness)

calcTotalNumOfGenomes <- function(simulation){
  counter <- 0
  for(i in seq(1,length(simulation$Pool$species))){
    counter <- counter + length(simulation$Pool$species[[i]]$genomes)
  }
  return(counter)
}

simulationRunner <- function(simulation,speciesNum,genomeNum,plotScene, pctSimulated,framesPerSecond=1){
  i<-speciesNum
  j <-genomeNum
  if(length(simulation$Pool$species[[i]]$genomes[[j]]$ConnectionGenes)>0){
    if(plotScene){
      tCount <- 0
      while(tCount < config.video.phenotypedurationseconds){
        tryCatch({drawPhenotypeNEAT(simulation$Pool$species[[i]]$genomes[[j]],simulation$Config)},error=function(e) print(paste("Error during draw:",e)))
        tCount <- tCount + 1/framesPerSecond
      }
    }
    #print(paste("Started simulation of species",i,"/",length(simulation$Pool$species),"genome",j,"/",length(simulation$Pool$species[[i]]$genomes)))
    state <- simulation$ProcessInitialStateFunc()
    fitness <- 0
    simulation$Pool$species[[i]]$genomes[[j]] <- generateNetwork(simulation$Pool$species[[i]]$genomes[[j]],simulation$Config)
    simulation$Pool$species[[i]]$genomes[[j]]$Fitness <- 0
    #Repeat acts like a do-while loop
    frameNum <- 0
    repeat{
      if(plotScene){
        simulation$PlotState(state)
      }
      neuralNetInputs <- simulation$ProcessStateToNeuralInputFunc(state)

      #  print(simulation$Pool$species[[i]]$genomes[[j]])
      pkg.env$debugGenome <-simulation$Pool$species[[i]]$genomes[[j]]
      neuralNetOutputs <- evaluateNetwork(simulation$Pool$species[[i]]$genomes[[j]]$Network,neuralNetInputs,simulation$Config)
      updatedState<-simulation$ProcessUpdateStateFunc(state,neuralNetOutputs)
      updatedFitness <- simulation$FitnessUpdateFunc(state,updatedState,simulation$Pool$species[[i]]$genomes[[j]]$Fitness)
      simulation$Pool$species[[i]]$genomes[[j]]$Fitness <- updatedFitness


      if(simulation$TerminationCheckFunc(frameNum,state,updatedState,fitness,updatedFitness)){
        break
      }

      fitness <- updatedFitness
      state <- updatedState
      frameNum <- frameNum + 1
    }
    simulation$Pool$species[[i]]$genomes[[j]]$Fitness <- updatedFitness
    print(paste(round(pctSimulated,2),"% Finished simulation of species",i,"/",length(simulation$Pool$species),"genome",j,"/",length(simulation$Pool$species[[i]]$genomes),"with fitness",updatedFitness))
    #print(paste(round(pctSimulated,2),"% Finished simulation of species",i,"genome",j,"with fitness",updatedFitness))
  } else {
    print(paste("Skipped simulation of species",i,"genome",j,"due to no connection genes"))
  }

  if(plotScene){
    tCount <- 0
    while(tCount < config.video.performancedurationseconds){
      plotPerformanceTracker(simulation$PerformanceTracker)
      tCount <- tCount + 1/framesPerSecond
    }
  }
  return(simulation)
}

pkg.env$debugGenome <- 0

#' Runs a single generation
#'
#' Takes in a simulation, runs all the genomes, evaluates fitness and breeds the new generation
#' @param simulation Takes a NEATSimulation class
#' @param createVideo True/False to save a video of the highest fitness simulation
#' @param videoPath Path to where to save the video
#' @param videoName Name of the video
#' @param framesPerSecond The frames per second of the video
#' @return NEATSimulation class with new generation of genomes
#' @export
NEATSimulation.RunSingleGeneration <- function(simulation, createVideo=F, videoPath="videos",videoName="", framesPerSecond=1){
  assertTrueFunc(is(simulation,"NEATSimulation"),"simulation must be a of class NEATSimulation")
  oldMaxFitness <- simulation$Pool$maxFitness

  print("Starting simulations...")
  counter <- 1
  nTot <- calcTotalNumOfGenomes(simulation)
  for(i in seq(1,length(simulation$Pool$species))){
    for(j in seq(1,length(simulation$Pool$species[[i]]$genomes))){
      simulation <- simulationRunner(simulation,i,j,F,100*counter/nTot)
      counter <- counter + 1
    }
  }
  simulation$Pool <- newGeneration(simulation$Pool,simulation$Config)
  #print(paste("MaxFitness:",simulation$Pool$maxFitness))

  simulation$PerformanceTracker[simulation$Pool$generation,] <- c(simulation$Pool$generation,simulation$Pool$minFitness,simulation$Pool$maxFitness,simulation$Pool$meanFitness,simulation$Pool$medianFitness)
  print(simulation$PerformanceTracker[simulation$Pool$generation,])
  if(createVideo){
    createdVid <- F
    for(i in seq(1,length(simulation$Pool$species))){
      for(j in neatseq(1,length(simulation$Pool$species[[i]]$genomes))){
        if(simulation$Pool$species[[i]]$genomes[[j]]$Fitness==simulation$Pool$maxFitness & !createdVid){
          if(oldMaxFitness != simulation$Pool$maxFitness){
            videoName = paste(videoPath,"/",videoName,"generation",simulation$Pool$generation,"fitness",simulation$Pool$maxFitness,"species",i,"genome",j,".mpeg",sep="")
            print(paste("Creating video",videoName,"..."))
            oopt = ani.options(ani.width = 1200, ani.height = 800, other.opts = "-define png:color-type=2")
            saveVideo(simulationRunner(simulation,i,j,T,100,framesPerSecond),interval=1/framesPerSecond,ani.options=oopt,video.name=videoName)
            ani.options(oopt)
          } else {
            print("Max Fitness did not increase on this generation so skipping creating a video")
          }
          createdVid <- T
        }
      }
    }
  }

  return (simulation)
}

#' Runs a genome and tracks the state history (will run the most fit by default)
#'
#' Runs a genome and tracks the state history (will run the most fit by default)
#' @param simulation Takes a NEATSimulation class
#' @param genomeNum the genome number to run
#' @param speciesNum the species number to run
#' @return State history
#' @export
NEATSimulation.GetStateHistoryForGenomeAndSpecies <- function(simulation, genomeNum=NA, speciesNum=NA){
    assertTrueFunc(is(simulation,"NEATSimulation"),"simulation must be a of class NEATSimulation")
 if(is.na(genomeNum) || is.na(speciesNum)){
    print("No genome or species specified, using the most fit genome")
    genome <- findMostFitGenome(simulation)
 } else {
    assertTrueFunc(is(genomeNum,"numeric"),"genomeNum must be numeric")
    assertTrueFunc(is(speciesNum,"numeric"),"speciesNum must be numeric")
    genome <- simulation$Pool$species[[i]]$genomes[[j]]
    assertTrueFunc(is(genome,"genome"),"Invalid species or genome number")
 }


    state <- simulation$ProcessInitialStateFunc()
    stateHist <- as.data.frame(t(unlist(state)))
    fitness <- 0
    genome$Fitness <- fitness
    genome <- generateNetwork(genome,simulation$Config)

    #Repeat acts like a do-while loop
    frameNum <- 0
    repeat{

      neuralNetInputs <- simulation$ProcessStateToNeuralInputFunc(state)
      neuralNetOutputs <- evaluateNetwork(genome$Network,neuralNetInputs,simulation$Config)
      updatedState<-simulation$ProcessUpdateStateFunc(state,neuralNetOutputs)
      updatedFitness <- simulation$FitnessUpdateFunc(state,updatedState,genome$Fitness)
      genome$Fitness <- updatedFitness

      stateHist<-rbind(stateHist,unlist(updatedState))

      if(simulation$TerminationCheckFunc(frameNum,state,updatedState,fitness,updatedFitness)){
        break
      }

      fitness <- updatedFitness
      state <- updatedState
      frameNum <- frameNum + 1
    }
    genome <- updatedFitness

  return((stateHist))

}

#' @export
plot.NEATSimulation <- function(simulation){
  plotPerformanceTracker(simulation$PerformanceTracker)
}

plotPerformanceTracker <- function(data){
  plot(x=data[,"generation"],y=data[,"maxFitness"],col="blue",main="Fitness",xlab="Generation",ylab="Fitness",type="o",ylim=c(0,max(data[,"maxFitness"])),lwd=2)
  lines(x=data[,"generation"],y=data[,"minFitness"],col="red",type="o",lwd=2)
  lines(x=data[,"generation"],y=data[,"meanFitness"],col="green",type="o",lwd=2)
  lines(x=data[,"generation"],y=data[,"medianFitness"],col="purple",type="o",lwd=2)
  legend(x='bottomright', c("Min","Max","Mean","Median"),  fill=c("red","blue","green","purple"), bty='n')

}















