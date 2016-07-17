
#' @importFrom grDevices dev.new
#' @importFrom graphics box layout legend lines mtext par plot polygon symbols text
#' @importFrom igraph add_edges add_vertices layout.circle make_empty_graph



config.nodegene.drawwidth <- 60
config.nodegene.drawheight <- 35
config.nodegene.textmargin <- 3
config.nodegene.radius <- 5
config.nodegene.nodemargin <- 5*config.nodegene.radius

config.connectiongene.drawwidth <- 60
config.connectiongene.drawheight <- 60
config.connectiongene.textmargin <- 3

# Function to create a blank canvas / scene for drawing objects onto later
#' @export
# @param bottomLeftX Bottom left x co-ordinate of the canvas
# @param bottomLeftY Bottom left y co-ordinate of the canvas
# @param width Canvas width
# @param height Canvas height
# @param main Text title of the plot
# @param xlab Label for x-axis
# @param ylab Label for y-axis
# @param ann See plot ann
# @param xaxt See plot xaxt
# @param yaxt See plot yaxt
# @param xlim See plot xlim
# @param ylim See plot ylim
# @param frame.plot Bool to enable or disable drawing of a frame around the canvas
createSceneFunc <- function(bottomLeftX, bottomLeftY, width,height,main="",xlab="",ylab="",ann=T,xaxt=NULL,yaxt=NULL,xlim=NULL,ylim=NULL,frame.plot=T){
  plot(c(bottomLeftX, bottomLeftX+width), c(bottomLeftY,bottomLeftY+height), type = "n",ann=ann, xaxt=xaxt, yaxt=yaxt,xlim=xlim,ylim=ylim,main=main,xlab=xlab,ylab=ylab,frame.plot=frame.plot )
}

# Function to draw a box on the scene
# @param topLeftX x co-ordinate of top left corner of box
# @param topLeftY y co-ordinate of top left corner of box
# @param width Width of the box
# @param height Height of the box
# @param fillColour Colour to fill the box in with
# @param borderColour Colour of the box edge
#' @export
createBoxFunc <- function(topLeftX, topLeftY, width, height, fillColour=NA, borderColour="black"){
  polygon(c(topLeftX,topLeftX+width,topLeftX+width,topLeftX),
          c(topLeftY,topLeftY,topLeftY-height,topLeftY-height),
          col = fillColour, border=borderColour)
}

# Function to draw a circle on the scene
# @param centerX x co-ordinate of the center of the circle
# @param centerY y co-ordinate of the center of the circle
# @param radius Radius of the circle
# @param fillColour Colour to fill the circle in with
# @param borderColour Colour of the circle edge
#' @export
createCircleFunc <- function(centerX,centerY,radius,fillColour=NA, borderColour="black"){
  symbols(centerX,centerY,circles=radius,inches=F,add=T,fg=borderColour,bg=fillColour)
}

#Function to write text scene
createTextFunc <- function(topLeftX, topLeftY, width, height, message,adjustment=c(0,1),textcolour="black"){
  text(topLeftX,topLeftY,message,adj=adjustment,col=textcolour)
}


drawGenotypeNEAT <- function(object,neatConfig,x,y){
  UseMethod("drawGenotypeNEAT")
}

drawGenotypeNEAT.connectiongene <- function(connectiongene,neatConfig,topLeftX,topLeftY){
  #print("drawing connectiongene")
  createBoxFunc(topLeftX,topLeftY,config.connectiongene.drawwidth,config.connectiongene.drawheight)

  txtSpace <- config.connectiongene.drawheight / 6
  txtColour <- "black"
  if(!connectiongene$Enabled){
    txtColour <- "red"
  }
  createTextFunc(topLeftX+config.connectiongene.textmargin,topLeftY-config.connectiongene.textmargin-txtSpace*0,message=paste("InNode:",connectiongene$InNode),textcolour=txtColour)
  createTextFunc(topLeftX+config.connectiongene.textmargin,topLeftY-config.connectiongene.textmargin-txtSpace*1,message=paste("OutNode:",connectiongene$OutNode),textcolour=txtColour)
  createTextFunc(topLeftX+config.connectiongene.textmargin,topLeftY-config.connectiongene.textmargin-txtSpace*2,message=paste("Weight:",round(connectiongene$Weight,5)),textcolour=txtColour)
  createTextFunc(topLeftX+config.connectiongene.textmargin,topLeftY-config.connectiongene.textmargin-txtSpace*3,message=paste("Enabled:",connectiongene$Enabled),textcolour=txtColour)
  createTextFunc(topLeftX+config.connectiongene.textmargin,topLeftY-config.connectiongene.textmargin-txtSpace*4,message=paste("Innovation:",connectiongene$Innovation),textcolour=txtColour)

}

drawGenotypeNEAT.genome <- function(genome,neatConfig,topLeftX=0,topLeftY=0){
  #print("drawing genome")

  connectionOffset <- 5

  nodeGenesHeight <- config.nodegene.drawheight
  nodeGenesWidth <- (sum(unlist(lapply(generateNetwork(genome,neatConfig)$Network$Neurons,function(x) { return(!is.null(x))})))
  )*config.connectiongene.drawwidth
  connectionGenesWidth <- max(length(genome$ConnectionGenes)*config.connectiongene.drawwidth)
  connectionGenesHeight <- config.nodegene.drawheight+connectionOffset+config.connectiongene.drawheight

  width <- max(nodeGenesWidth,connectionGenesWidth)
  height <- max(nodeGenesHeight,connectionGenesHeight)

  par(mar=c(1,1,1,1) + 0.0)
  createSceneFunc(0,-height,width,height,xaxt="n",yaxt="n")


  networkNeurons <- generateNetwork(genome,neatConfig)$Network$Neurons
  drawCount <- 1
  for(i in seq(1,length(networkNeurons))){
    networkNeuron <- networkNeurons[[i]]
    if(!is.null(networkNeuron)){
      nodeType <- "Hidden"
      if(i <= neatConfig$Inputs){
        nodeType <- "Input"
      }
      if(i >= neatConfig$MaxNodes-neatConfig$Outputs){
        nodeType <- "Output"
      }

      createBoxFunc((drawCount-1)*config.connectiongene.drawwidth,topLeftY,config.nodegene.drawwidth,config.nodegene.drawheight)
      txtSpace <- config.nodegene.drawheight / 3
      createTextFunc((drawCount-1)*config.connectiongene.drawwidth+config.nodegene.textmargin,topLeftY-config.nodegene.textmargin-txtSpace*0,message=paste("NodeId:",i),textcolour="black")
      createTextFunc((drawCount-1)*config.connectiongene.drawwidth+config.nodegene.textmargin,topLeftY-config.nodegene.textmargin-txtSpace*1,message=paste("NodeType:",nodeType),textcolour="black")
      drawCount <- drawCount +1
    }

  }

  for(i in seq(1,length(genome$ConnectionGenes))){
    drawGenotypeNEAT(genome$ConnectionGenes[[i]],neatConfig,(i-1)*config.connectiongene.drawwidth,-config.nodegene.drawheight-connectionOffset)
  }

}


drawPhenotypeNEAT <- function(object,neatConfig,x,y){
  UseMethod("drawPhenotypeNEAT")
}

drawPhenotypeNEAT.nodegene <- function(nodegene,neatConfig,topLeftX,topLeftY){
  centerOffset <- config.nodegene.nodemargin+config.nodegene.radius
  centerX <- topLeftX + centerOffset
  centerY <- topLeftY - centerOffset
  fillColour <- "white"
  if(nodegene$NodeType=="input"){ fillColour <- "cadetblue1" }
  if(nodegene$NodeType=="output"){ fillColour <- "coral" }
  createCircleFunc(centerX,centerY,radius=config.nodegene.radius,fillColour)
  createTextFunc(centerX,centerY,message=nodegene$NodeId,adjustment=c(0.5,0.5))
}

drawPhenotypeNEAT.connectiongene <- function(connectiongene,neatConfig,topLeftX,topLeftY){
  stop("Not implemented, connections drawn in the drawPhenotypeNEAT.genome function")
}


drawPhenotypeNEAT.genome <- function(genome,neatConfig,topLeftX,topLeftY){

  #inputNodeGenes <- Filter(function(x) { return(x$NodeType=="input")},genomeOne$NodeGenes)
  #hiddenNodeGenes <- Filter(function(x) { return(x$NodeType=="hidden")},genomeOne$NodeGenes)
  #outputGenes <- Filter(function(x) { return(x$NodeType=="output")},genomeOne$NodeGenes)

  #http://www.shizukalab.com/toolkits/sna/plotting-directed-networks
  #set.seed(1)


  gg <- createGraph(genome,neatConfig)

  #cleanLayout <- layout.fruchterman.reingold(gg)
  cleanLayout <- layout.circle(gg)
  #  cleanLayout[c(inputsIdx,outputsIdx),] <- ll[c(inputsIdx,outputsIdx),]
  #    plot(gg,layout=ll)
  plot(gg,layout=cleanLayout)
  #print(gg)


}

createGraph <- function(genome,neatConfig){
  g <- make_empty_graph()
  nodeIdToVertexIdMap <- list()

  drawnNodes <- list()
  counter <- 1


  neurons <- generateNetwork(genome,neatConfig)$Network$Neurons
  for(i in seq(1,length(neurons))){
    vertexColour <- "white"
    if(i <= neatConfig$Inputs){  vertexColour<-"cadetblue1"}
    if(i > neatConfig$MaxNodes){  vertexColour<-"coral"}
    if(!is.null(neurons[[i]])){
      g<-add_vertices(g,1, color = vertexColour,name=i)
      nodeIdToVertexIdMap[i] <- counter
      counter <- counter+1
    }
  }

  for(i in seq(1,length(neurons))){
    neuron <- neurons[[i]]
    if(!is.null(neuron) && !is.null(neuron$Incoming)){
      for(connectionGene in neuron$Incoming){
        edgeColour <- "black"
        if(!connectionGene$Enabled){
          edgeColour <- "red"
        }
        g<-add_edges(g,color=edgeColour,c(nodeIdToVertexIdMap[connectionGene$InNode],nodeIdToVertexIdMap[connectionGene$OutNode]),attr=list(weight=c(connectionGene$Weight)))
      }
    }
  }
  return (g)
}


drawNEAT <- function(object,neatConfig){
  UseMethod("drawNEAT")
}

drawNEAT.genome <- function(genome,neatConfig){
  dev.new()
  layout(matrix(c(1,2), 2, 1, byrow = TRUE),heights=c(1,2))
  drawGenotypeNEAT(genome,neatConfig)
  drawPhenotypeNEAT(genome,neatConfig)
}

drawSideBySideNEAT <- function(genomeOne,genomeTwo,neatConfig,chartDescription){
  dev.new()
  layout(matrix(c(5,5,1,3,2,4), 3, 2, byrow = TRUE),heights=c(1,4,8))

  drawGenotypeNEAT(genomeOne,neatConfig)
  drawPhenotypeNEAT(genomeOne,neatConfig)

  box(lty = '1373', col = 'red')
  drawGenotypeNEAT(genomeTwo,neatConfig)
  drawPhenotypeNEAT(genomeTwo,neatConfig)
  box(lty = '1373', col = 'red')
  createSceneFunc(0,-0,0,0,ann=F,xaxt="n",yaxt="n",frame.plot=F)
  mtext(chartDescription, outer = F, cex = 1.5,side=1)

}

drawCrossoverNEAT <- function(genomeA,genomeB,neatConfig,genomeCross){
  dev.new()
  layout(matrix(c(5,5,5,4,4,4,1,2,3), 3, 3, byrow = TRUE),heights=c(1,6,10))

  par(mar=c(1,1,1,1))
  drawPhenotypeNEAT(genomeA,neatConfig)
  box(lty = '1373', col = 'red')
  par(mar=c(1,1,1,1))
  drawPhenotypeNEAT(genomeB,neatConfig)
  box(lty = '1373', col = 'red')
  par(mar=c(1,1,1,1))
  drawPhenotypeNEAT(genomeCross,neatConfig)
  box(lty = '1373', col = 'blue')


  smalloffset <- 5
  largeoffset <-10
  innovationsA <- unlist(lapply(genomeA$ConnectionGenes,function(x){x$Innovation}))
  innovationsB <- unlist(lapply(genomeB$ConnectionGenes,function(x){x$Innovation}))
  innovationsCross <- unlist(lapply(genomeCross$ConnectionGenes,function(x){x$Innovation}))
  uniqueInnovations <- unique(sort(c(innovationsA,innovationsB)))
  #print(innovationsA)
  #print(innovationsB)
  #print(uniqueInnovations)
  width <- max(length(uniqueInnovations),length(innovationsCross))*config.connectiongene.drawwidth
  height <- 3*config.connectiongene.drawheight+smalloffset+largeoffset

  #print(width)
  #print(height)
  par(mar=c(1,1,1,1) + 0.0)
  createSceneFunc(0,-height,width,height,xaxt="n",yaxt="n")


  for(i in seq(1,length(genomeA$ConnectionGenes))){
    drawGenotypeNEAT(genomeA$ConnectionGenes[[i]],neatConfig,(match(genomeA$ConnectionGenes[[i]]$Innovation,uniqueInnovations)-1)*config.connectiongene.drawwidth,0)
  }
  for(i in seq(1,length(genomeB$ConnectionGenes))){
    drawGenotypeNEAT(genomeB$ConnectionGenes[[i]],neatConfig,(match(genomeB$ConnectionGenes[[i]]$Innovation,uniqueInnovations)-1)*config.connectiongene.drawwidth,-config.connectiongene.drawheight-smalloffset)
  }
  for(i in seq(1,length(genomeCross$ConnectionGenes))){
    drawGenotypeNEAT(genomeCross$ConnectionGenes[[i]],neatConfig,(i-1)*config.connectiongene.drawwidth,-2*config.connectiongene.drawheight-smalloffset-largeoffset)
  }

  createSceneFunc(0,-0,0,0,ann=F,xaxt="n",yaxt="n",frame.plot=F)
  mtext("Genome Crossover - Combining the topological features", outer = F, cex = 1.5,side=1)
}

#inNodeOne <- newConnectiongene()
#inNodeOne$InNode <- 1
#inNodeOne$OutNode <- 3
#inNodeOne$Enabled = T
#inNodeOne$Weight=7
#inNodeOne$Innovation=1
#
#hiddenNodeOne <- newConnectiongene()
#hiddenNodeOne$InNode <- 2
#hiddenNodeOne$OutNode <- 3
#hiddenNodeOne$Enabled = T
#hiddenNodeOne$Weight=5
#hiddenNodeOne$Innovation=2
#
#hiddenNodeTwo <- newConnectiongene()
#hiddenNodeTwo$InNode <- 2
#hiddenNodeTwo$OutNode <- config.MaxNodes+1
#hiddenNodeTwo$Enabled = T
#hiddenNodeTwo$Weight=3
#hiddenNodeTwo$Innovation=4
#
#outNodeOne <- newConnectiongene()
#outNodeOne$InNode <- 3
#outNodeOne$OutNode <- config.MaxNodes+1
#outNodeOne$Enabled = T
#outNodeOne$Weight=3
#outNodeOne$Innovation=3
#
#genomeOne = newGenome()
#genomeOne$ConnectionGenes <-list(inNodeOne,outNodeOne,hiddenNodeOne,hiddenNodeTwo)
#genomeOne$MaxNeuron <- 5
#genomeOne
#drawNEAT(genomeOne)
#
#genomeOne <- generateNetwork(genomeOne)
#
#netInputs <- list()
#netInputs[1] <- 1.2
#netInputs[2] <- 2.5
#evaluateNetwork(genomeOne$Network,netInputs)
#
##expect to see  (1.2*7+2.5*5)*3 = 62.7 (if override sigmoid with sum)
##drawNEAT(mutateGenome(genomeOne))
#innovationGlob <- 4; t <-nodeMutate(genomeOne); drawNEAT(t)
#innovationGlob <- 4; t<- enableDisableMutate(genomeOne,F); drawNEAT(t)
#innovationGlob <- 4; t<-pointMutate(genomeOne); drawNEAT(t)
#
#
#
#demoConnection <- newConnectiongene()
#demoConnection$InNode <- 2
#demoConnection$OutNode <- config.MaxNodes+1
#demoConnection$Enabled = T
#demoConnection$Weight=5
#demoConnection$Innovation=1
#innovationGlob <- 1
#
#demoGenome = newGenome()
#demoGenome$ConnectionGenes <-list(demoConnection)
#demoGenome$MaxNeuron <- 2
#demoGenome
#
#drawNEAT(demoGenome)
#innovationGlob <- 2; drawNEAT(linkMutate(demoGenome,T)) #Adds a new connection
#innovationGlob <- 2; drawNEAT(pointMutate(demoGenome)) #Changes the weight in the network
#innovationGlob <- 2; drawNEAT(nodeMutate(demoGenome)) #Adds a new node
#innovationGlob <- 2; drawNEAT(enableDisableMutate(demoGenome,F))
#
#
#set.seed(2)
#imgWidth <- 1600
#imgHeight <- 1200
##png(filename="linkMutate.png",width=imgWidth,height=imgHeight)
#innovationGlob <- 1; mutatedGenome<-linkMutate(demoGenome,T)
#drawSideBySideNEAT(demoGenome,mutatedGenome,"Link Mutate - Add a new connection")
##dev.off()
#
##png(filename="pointMutate.png",width=imgWidth,height=imgHeight)
#innovationGlob <- 1; mutatedGenome<-pointMutate(demoGenome)
#drawSideBySideNEAT(demoGenome,mutatedGenome,"Point Mutate - Mutate the weights")
##dev.off()
#
##png(filename="nodeMutate.png",width=imgWidth,height=imgHeight)
#innovationGlob <- 1; mutatedGenome<-nodeMutate(demoGenome)
#drawSideBySideNEAT(demoGenome,mutatedGenome,"Node Mutate - Add a new node, by replacing a connection with equivalent connections")
##dev.off()
#
##png(filename="enableDisableMutate.png",width=imgWidth,height=imgHeight)
#innovationGlob <- 1; mutatedGenome<-enableDisableMutate(demoGenome,F)
#drawSideBySideNEAT(demoGenome,mutatedGenome,"Enable/Disable Mutate - Enables/Disables a connection")
##dev.off()





