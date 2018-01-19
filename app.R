install.packages("DiGGer_1.0.1_R_x86_64-pc-linux-gnu.tar.gz", repos = NULL)#, type = "win.binary")

require(shiny)
require(ggplot2)
require(DiGGer)
require(magrittr)

ui <- fluidPage(
      sidebarLayout(
            
            sidebarPanel(
                   numericInput("nChecks", "Number of Checks", 3, 2, 10, 1),
                   numericInput("rowNumber", "Number of Rows", 7, 3, 30, 1),
                   numericInput("colNumber", "Number of Columns", 7, 3, 30, 1),  
            
                   conditionalPanel(
                   condition = "input.rowNumber != input.colNumber",
                         selectInput("bDim", "Balanced block factor", list("Rows", "Columns"))
                   ),       
                   actionButton("Go", "Create design")
            ),
        
        
            mainPanel(
                plotOutput("lay")     
            )  
            
            
        
        
        
      )
)


server <- function(input, output){
  
### Optimization function  
Dig_R <- function(nRows,nCols,nChecks,bDim="BOTH",Zvar=1,
                  target=0,maxE_1=10000,maxE_2=10000, Intensity1=100,Intensity2=100){
    
    require(DiGGer)
    
    #Define balanced block factor 
    if(bDim=="BOTH"|bDim=="COLUMN")
    { 
      bBlock <- nCols
      ibDim  <- "ROW"
    } else if(bDim=="ROW"){
      bBlock <- nRows
      ibDim  <- "COLUMN"
    }
    
    #Calculate treatments numbers  
    entries <- nRows*nCols-nChecks*bBlock
    nTrts   <- entries+nChecks 
    
    repeat{                                           #Repeat loop to reach target value
      
      #Standard Initial Design 
      D <- DiGGer(nTrts,nRows,nCols,nRows,nCols, 
                  treatRepPerRep = rep(c(bBlock,1),c(nChecks,entries)),
                  treatGroup = rep(c(1,2),c(nChecks,entries)))
      
      #Two objectives, one for each phase
      O1 <- Objective(1,"NONE",0,corr = Correlation("ID",0,"ID",0))
      O2 <- Objective(1,"NONE",0,corr = Correlation("ID",0,"ID",0))
      
      #Add blocks to the objectives, random 1st phase and fixed 2nd
      if(bDim=="ROW"|bDim=="COLUMN"){
        addBlock.Objective(O1,newblock = Block(1,nCols,Zvar))
        addBlock.Objective(O1,newblock = Block(nRows,1,Zvar))
        addBlock.Objective(O2,newblock = Block(1,nCols,0))
        addBlock.Objective(O2,newblock = Block(nRows,1,0))
      }
      else{
        addBlock.Objective(O1,newblock = Block(1,nCols,0))
        addBlock.Objective(O1,newblock = Block(nRows,1,0))
        addBlock.Objective(O2,newblock = Block(1,nCols,0))
        addBlock.Objective(O2,newblock = Block(nRows,1,0))
      }
      
      #Define two phases, 1st Agg then A22,
      P1 <-Phase(nRows,nCols,nRows,nCols,aType="Agg",
                 objectives = list(O1),maxInterchanges = maxE_1,     #Search Intensity and 
                 searchIntensity = Intensity1, targetAValue = target)#max Exchanges can be
      P2 <-Phase(nRows,nCols,nRows,nCols,aType="A22",          #adjusted;
                 objectives = list(O2), maxInterchanges = maxE_2,
                 searchIntensity = Intensity2)
      
      #Add the two phases and run search
      addPhase(D,P1)
      addPhase(D,P2)
      R <- run(D)
      
      #Break the loop only if design does not pass over target value
      ifelse(R$ddphase[[1]]$aMeasures[1]<(target-0.00001), NA,break) 
    }
    
    #Adjust Design File
    DesF <- R$dlist
    DesF$COLUMN                 <- DesF$RANGE %>% as.factor
    DesF$ROW                    <- DesF$ROW %>% as.factor
    DesF$TRT                    <- DesF$TRT %>% as.factor
    keeps                       <- c("ROW","COLUMN","TRT")
    DesF                        <- DesF[keeps]
    L <- R$ddphase[[1]]$design
    Agg <- R$ddphase[[1]]$aMeasures[1]
    
    #Return layout, Design File and Agg
    output <- list(L,DesF,Agg)
    return(output)
  } 

## Event reactive plot creation 
CreateLayout <-   eventReactive(input$Go, {ggplot(data=Dig_R(input$rowNumber, input$colNumber,input$nChecks)[[2]], aes(ROW, COLUMN)) +
                    theme(axis.line = element_line(size=0.5), panel.grid = element_blank(),
                    axis.text=element_text(size=10), 
                    plot.title=element_text(hjust = 0.5,size=15,face = "bold"),
                    panel.background = element_rect(fill = 'white', colour = 'white'),
                    axis.title = element_text(size=10,face="bold"))+
                    geom_tile(colour="black",size=0.7) + 
                    geom_text(aes(ROW, COLUMN, label = TRT), size = 5) +
                    theme(legend.position="none")+ coord_flip()+
                    scale_x_discrete(name="Row",expand = c(0, 0))+
                    scale_y_discrete(name="Column",expand = c(0, 0))
                    })

output$lay <- renderPlot(CreateLayout()) ##call reactive function to render plot
}

shinyApp(ui = ui, server = server)


### Deploy app ###
rsconnect::setAccountInfo(name='buzzwinkel', token='D2B6A38F5FBE4C51B5BACEC30E2DA66F', secret='ggsXQnt15mHBUfGefD9F8vd701FBhjqdGAkzlDGn')
rsconnect::deployApp("F:/GitHub/ARCD", appTitle = "ARCD Creator")
