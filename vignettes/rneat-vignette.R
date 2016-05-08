## ----eval=FALSE----------------------------------------------------------
#  
#  drawPoleFunc <- function(fixedEnd.x,fixedEnd.y,poleLength, theta,
#                           fillColour=NA, borderColour="black"){
#    floatingEnd.x <- fixedEnd.x-poleLength * sin(theta)
#    floatingEnd.y <- fixedEnd.y+poleLength * cos(theta)
#  
#    polygon(c(fixedEnd.x,floatingEnd.x,floatingEnd.x,fixedEnd.x),
#            c(fixedEnd.y,floatingEnd.y,floatingEnd.y,fixedEnd.y),
#            col = fillColour, border=borderColour)
#  }
#  
#  drawPendulum <- function(fixedEnd.x,fixedEnd.y,poleLength, theta,
#                           radius,fillColour=NA, borderColour="black"){
#    floatingEnd.x <- fixedEnd.x-poleLength * sin(theta)
#    floatingEnd.y <- fixedEnd.y+poleLength * cos(theta)
#    createCircleFunc(floatingEnd.x,floatingEnd.y,radius,fillColour,borderColour)
#  }
#  
#  #Parameters to control the simulation
#  simulation.timestep = 0.005
#  simulation.gravity = 9.8 #meters per second^2
#  simulation.numoftimesteps = 2000
#  
#  pole.length = 1 #meters, total pole length
#  pole.width = 0.2
#  pole.theta = pi
#  pole.thetaDot = 0
#  pole.thetaDotDot = 0
#  pole.colour = "purple"
#  
#  
#  pendulum.centerX = NA
#  pendulum.centerY = NA
#  pendulum.radius = 0.1
#  pendulum.mass = 0.1
#  pendulum.colour = "purple"
#  
#  cart.width=0.5
#  cart.centerX = 0
#  cart.centerY = 0
#  cart.height=0.2
#  cart.colour="red"
#  cart.centerXDot = 0
#  cart.centerXDotDot = 0
#  cart.mass = 0.4
#  cart.force = 0
#  cart.mu=2
#  
#  
#  track.limit= 10 #meters from center
#  track.x = -track.limit
#  track.height=0.01
#  track.y = 0.5*track.height
#  track.colour = "blue"
#  
#  leftBuffer.width=0.1
#  leftBuffer.height=0.2
#  leftBuffer.x=-track.limit-0.5*cart.width-leftBuffer.width
#  leftBuffer.y=0.5*leftBuffer.height
#  leftBuffer.colour = "blue"
#  
#  rightBuffer.width=0.1
#  rightBuffer.height=0.2
#  rightBuffer.x=track.limit+0.5*cart.width
#  rightBuffer.y=0.5*rightBuffer.height
#  rightBuffer.colour = "blue"
#  
#  #Define the size of the scene (used to visualise what is happening in the simulation)
#  scene.width = 2*max(rightBuffer.x+rightBuffer.width,track.limit+pole.length+pendulum.radius)
#  scene.bottomLeftX = -0.5*scene.width
#  scene.height=max(pole.length+pendulum.radius,scene.width)
#  scene.bottomLeftY = -0.5*scene.height
#  
#  poleBalance.InitialState <- function(){
#    state <- list()
#    state[1] <- cart.centerX
#    state[2] <- cart.centerXDot
#    state[3] <- cart.centerXDotDot
#    state[4] <- cart.force
#    state[5] <- pole.theta
#    state[6] <- pole.thetaDot
#    state[7] <- pole.thetaDotDot
#    return(state)
#  }
#  
#  poleBalance.ConvertStateToNeuralNetInputs <- function(currentState){
#    return (currentState)
#  }
#  
#  poleBalance.UpdatePoleState <- function(currentState,neuralNetOutputs){
#    #print("Updating pole state")
#    #print(neuralNetOutputs)
#    cart.centerX <- currentState[[1]]
#    cart.centerXDot <- currentState[[2]]
#    cart.centerXDotDot <- currentState[[3]]
#    cart.force <- currentState[[4]]+neuralNetOutputs[[1]]
#    pole.theta <- currentState[[5]]
#    pole.thetaDot <- currentState[[6]]
#    pole.thetaDotDot <- currentState[[7]]
#  
#    costheta = cos(pole.theta)
#    sintheta = sin(pole.theta)
#    totalmass = cart.mass+pendulum.mass
#    masslength = pendulum.mass*pole.length
#  
#    pole.thetaDotDot = (simulation.gravity*totalmass*sintheta+costheta*
#                      (cart.force-masslength*pole.thetaDot^2*sintheta-cart.mu*cart.centerXDot))/
#                      (pole.length*(totalmass-pendulum.mass*costheta^2))
#  
#    cart.centerXDotDot =(cart.force+masslength*(pole.thetaDotDot*costheta-pole.thetaDot^2*sintheta)-
#                        cart.mu*cart.centerXDot)/totalmass
#  
#    cart.centerX = cart.centerX+simulation.timestep*cart.centerXDot
#    cart.centerXDot = cart.centerXDot+simulation.timestep*cart.centerXDotDot
#    pole.theta = (pole.theta +simulation.timestep*pole.thetaDot )
#    pole.thetaDot = pole.thetaDot+simulation.timestep*pole.thetaDotDot
#  
#    currentState[1] <- cart.centerX
#    currentState[2] <- cart.centerXDot
#    currentState[3] <- cart.centerXDotDot
#    currentState[4] <- cart.force
#    currentState[5] <- pole.theta
#    currentState[6] <- pole.thetaDot
#    currentState[7] <- pole.thetaDotDot
#    return (currentState)
#  }
#  
#  
#  
#  poleBalance.UpdateFitness <- function(oldState,updatedState,oldFitness){
#    #return (oldFitness+1) #fitness is just how long we've ran for
#    #return (oldFitness+((track.limit-abs(updatedState[[1]]))/track.limit)^2)
#    #More reward for staying near middle of track
#  
#    height <- cos(updatedState[[5]]) #is -ve if below track
#    heightFitness <- max(height,0)
#    centerFitness <- (track.limit-abs(updatedState[[1]]))/track.limit
#    return (oldFitness+(heightFitness + heightFitness*centerFitness))
#  }
#  
#  poleBalance.CheckForTermination <- function(frameNum,oldState,updatedState,oldFitness,newFitness){
#    cart.centerX <- updatedState[[1]]
#    cart.centerXDot <- updatedState[[2]]
#    cart.centerXDotDot <- updatedState[[3]]
#    cart.force <- updatedState[[4]]
#    pole.theta <- updatedState[[5]]
#    pole.thetaDot <- updatedState[[6]]
#    pole.thetaDotDot <- updatedState[[7]]
#  
#    oldpole.theta <- oldState[[5]]
#    if(frameNum > 20000){
#      print("Max Frame Num Exceeded , stopping simulation")
#      return (TRUE)
#    }
#  
#    height <- cos(pole.theta)
#    oldHeight <- cos(oldpole.theta)
#    if(height==-1 & cart.force==0){
#      return(TRUE)
#    }
#  
#    if(oldHeight >= 0 & height < 0){
#      #print("Pole fell over")
#      return (TRUE)
#    }
#    if(cart.centerX < track.x | cart.centerX > (track.x+2*track.limit)){
#      #print("Exceeded track length")
#      return (TRUE)
#    } else {
#      return (FALSE)
#    }
#  }
#  
#  poleBalance.PlotState <-function(updatedState){
#    cart.centerX <- updatedState[[1]]
#    cart.centerXDot <- updatedState[[2]]
#    cart.centerXDotDot <- updatedState[[3]]
#    cart.force <- updatedState[[4]]
#    pole.theta <- updatedState[[5]]
#    pole.thetaDot <- updatedState[[6]]
#    pole.thetaDotDot <- updatedState[[7]]
#  
#    createSceneFunc(scene.bottomLeftX,scene.bottomLeftY,scene.width,scene.height,
#                    main="Simulation of Inverted Pendulum - www.gekkoquant.com",xlab="",
#                    ylab="",xlim=c(-0.5*scene.width,0.5*scene.width),
#                    ylim=c(-0.5*scene.height,0.5*scene.height))
#  
#    createBoxFunc(track.x,track.y,track.limit*2,track.height,track.colour)
#    createBoxFunc(leftBuffer.x,leftBuffer.y,leftBuffer.width,leftBuffer.height,leftBuffer.colour)
#    createBoxFunc(rightBuffer.x,rightBuffer.y,rightBuffer.width,
#                  rightBuffer.height,rightBuffer.colour)
#    createBoxFunc(cart.centerX-0.5*cart.width,cart.centerY+0.5*cart.height,cart.width,cart.height,
#                  cart.colour)
#    drawPoleFunc(cart.centerX,cart.centerY,2*pole.length,pole.theta,pole.colour)
#    drawPendulum(cart.centerX,cart.centerY,2*pole.length,pole.theta,pendulum.radius,pendulum.colour)
#  
#  }
#  
#  config <- newConfigNEAT(7,1,500,50)
#  poleSimulation <- newNEATSimulation(config, poleBalance.InitialState,
#                                      poleBalance.UpdatePoleState,
#                                      poleBalance.ConvertStateToNeuralNetInputs,
#                                      poleBalance.UpdateFitness,
#                                      poleBalance.CheckForTermination,
#                                      poleBalance.PlotState)
#  
#  nMax <- 1 #Number of generations to run
#  for(i in seq(1,nMax)){
#    poleSimulation <- NEATSimulation.RunSingleGeneration(poleSimulation)
#    #poleSimulation <- NEATSimulation.RunSingleGeneration(poleSimulation,T,"videos",
#    #                                            "poleBalance",1/simulation.timestep)
#  }
#  
#  
#  
#  
#  

