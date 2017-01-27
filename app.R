library(shiny)
library(shinydashboard)
library(shinyFiles)
library(shinyDND)
library(boot)
library(splines)
library(DiagrammeR)
library(DT)


g.function<-function(dat, interv, outcome,outcome.model,model.t, bootn){
  dat$int<-rep(-1,nrow(dat))
  dat.notreat<-dat
  dat.notreat$int<-rep(0,nrow(dat.notreat))
  dat.notreat[,interv]<-rep(0,nrow(dat.notreat))
  dat.notreat[,outcome]<-rep(NA, nrow(dat.notreat))
  dat.treat<-dat
  dat.treat$int<-rep(1,nrow(dat.treat))
  dat.treat[,interv]<-rep(1, nrow(dat.treat))
  dat.treat[,outcome]<-rep(NA, nrow(dat.treat))
  full.dat<-as.data.frame(rbind(dat, dat.notreat, dat.treat))
  outcomemodel1<-glm(formula= outcome.model,family = model.t, data=full.dat)
  full.dat$meanOutcome<-predict(outcomemodel1, full.dat, type="response")
  results1<-with(full.dat, tapply(meanOutcome, list(int),mean))
  
  boot.function<-function(dat2,index,type){
    boot.dat<-dat2[index,]
    fit<-glm(formula= outcome.model , data = boot.dat)
    boot.dat$meanOutcome<-predict(fit, boot.dat, type="response")
    sd.mean<-with(boot.dat, tapply(meanOutcome, list(int), mean))
    results2<-as.numeric(unlist(c(sd.mean, diff(sd.mean[c(2,3)]))))
    return(results2)
  }
  
  goes<-bootn
  set.seed(8675309)
  results3<-boot(full.dat, boot.function, R=goes)
  results4<-matrix(NA,4,4)
  for(i in 1:4){
    results4[i,]<-c(results3$t0[i],
                    sd(results3$t[,i]),
                    results3$t0[i]-1.96 * sd(results3$t[,i]),
                    results3$t0[i]+1.96 * sd(results3$t[,i]))
  }
  
  dimnames(results4)<-list(c("Observed","No Intervention","Intervention","Difference (intervened vs not)"),
                           c("Mean","SD","Lower CI","Upper CI"))
  print(results4)
  
}

ui <- dashboardPage(
  dashboardHeader(title = "Parametric G-function"),
  dashboardSidebar(
    br(),
    sidebarMenu(id="tabs",
                menuItem("Inputs", tabName = "tab1")
    )
  ),
  dashboardBody(
    tabItems(
      tabItem(tabName = "tab1",
              box(title = "Upload Data",
                  fileInput("file1","Choose a file",accept=c('text/csv', 
                                                             'text/comma-separated-values,text/plain', 
                                                             '.csv')),
                  tags$hr(),
                  h5(helpText("Select file parameters:")),
                  checkboxInput(inputId = 'header1', label= 'Header', value= TRUE),
                  br(),
                  radioButtons(inputId = 'sep1', label = 'Seperator', choices = c(Comma=',', Semicolon=';', Tab='\t', Space= ' '), selected= ','),
                  actionButton("upload1", "Upload File")
              ),
              tabBox(title="Variable Selection",
                     tabPanel("Outcome",
                              uiOutput("g.out1"),
                              textInput("g.out2","Outcome Label:"),
                              selectInput("modelt","Family", choices=c("binomial","gaussian"), selected=2)
                     ),
                     tabPanel("Intervention",
                              uiOutput("g.expo1"),
                              textInput("g.expo2","Intervention Label:"),
                              textInput("g.cut1","Cut point for continuous interventions:")),
                     tabPanel("Confounders (Factors)",
                              uiOutput("g.cov1"),
                              uiOutput("g.cov2")),
                     tabPanel("Confounders (Numeric)",
                              uiOutput("g.cov3"),
                              uiOutput("g.cov4"),
                              numericInput("trans1","Include additional terms:", value=1, min = -10, max = 10),
                              uiOutput("trans2")
                              
                     )),
              box(title="DAG",
                  grVizOutput("og.dag1")),
              box(title="g-formula",
                  verbatimTextOutput("model1"),
                  numericInput("boot1","Number of Bootstrap resamples", value=50),
                  actionButton("go2","Calculate!"),
                  verbatimTextOutput("g.results1")))
      
      
    )
  ))

server <- shinyServer(function(input, output) {
  values<-reactiveValues()
  
  g.dat<-eventReactive(input$upload1,{
    file2<-input$file1
    if(is.null(file2)){return()}
    dat1<-read.table(file=file2$datapath, sep= input$sep1, header= input$header1, stringsAsFactors= FALSE)
  })
  
  output$g.cov1<-renderUI({
    df<-g.dat()
    if(is.null(df))return(NULL)
    items=colnames(df)
    names(items)=items
    selectInput("g.cov1.1","Select ALL Factor Confounders:",items, multiple=TRUE)
  })
  
  output$g.cov3<-renderUI({
    df<-g.dat()
    if(is.null(df))return(NULL)
    items=colnames(df)
    names(items)=items
    selectInput("g.cov3.1","Select ALL Numeric Confounders:",items, multiple=TRUE)
  })
  
  output$g.expo1<-renderUI({
    selectInput("g.expo1.1","Choose Exposure",choices = colnames(g.dat()))
  })
  
  output$g.out1<-renderUI({
    selectInput("g.out1.1","Choose Outcome",choices = colnames(g.dat()))
  })
  
  output$g.cov2<-renderUI({
    num.cov<-as.integer(length(input$g.cov1.1))
    lapply(1:num.cov, function(i){
      textInput(inputId =paste("g.cov2.1",i, sep=""),label=paste("Confounder",i,sep=" ") )
    })
  })
  
  output$g.cov4<-renderUI({
    num.cov<-as.integer(length(input$g.cov3.1))
    lapply(1:num.cov, function(i){
      textInput(inputId =paste("g.cov4.1",i, sep=""),label=paste("Confounder",i,sep=" ") )
    })
  })
  
  output$trans2<-renderUI({
    checkboxGroupInput("trans2.1","Additional Units for:",choices =c(unlist(input$g.cov3.1)),selected=NULL)
  })
  
  
  #######making the dag##########
  output$og.dag1 <- renderGrViz({
    cov.num<-length(input$g.cov1.1) +length(input$g.cov3.1)
    
    nodes_1<-create_nodes(
      nodes<-c(sapply(1:(length(input$g.cov1.1)), function(i){(input[[paste("g.cov2.1",i,sep="")]])}),
               sapply(1:(length(input$g.cov3.1)), function(i){(input[[paste("g.cov4.1",i,sep="")]])}))
    )
    
    edges_1<-create_edges(
      from = c(sapply(1:(length(input$g.cov1.1)), function(i){(input[[paste("g.cov2.1",i,sep="")]])}),
               sapply(1:(length(input$g.cov3.1)), function(i){(input[[paste("g.cov4.1",i,sep="")]])})),
      to = c(rep(input$g.expo2, cov.num))
    )
    
    
    g1 <- create_graph(
      nodes_df = nodes_1,
      edges_df = edges_1,
      graph_attrs = "rankdir = LR",
      node_attrs = "shape = plaintext")
    
    grViz(g1$dot_code)
    
    nodes_2<-create_nodes(
      nodes<-input$g.out2
    )
    
    edges_2.2<-create_edges(
      from = c(sapply(1:(length(input$g.cov1.1)), function(i){(input[[paste("g.cov2.1",i,sep="")]])}),
               sapply(1:(length(input$g.cov3.1)), function(i){(input[[paste("g.cov4.1",i,sep="")]])})),
      to = rep(input$g.out2, cov.num)
    )
    
    g2<-add_node_df(g1, nodes_2)
    g2.2<-add_edge(g2, input$g.expo2, input$g.out2)
    g2.3<-add_edge_df(g2.2, edges_2.2)
    grViz(g2.3$dot_code)
    
  })
  #######################################
  
  # g.form<-eventReactive(input$go2,{
  #   factor.pred<-c(unlist(input$g.cov1.1))
  #   numeric.pred<-c(unlist(input$g.cov3.1))
  #   if(is.null(input$trans2.1)){
  #     predictors<-c(input$g.expo1.1,factor.pred, numeric.pred)
  #   }
  #   if(!is.null(input$trans2.1)){
  #     other.numeric<-c(sapply(1:(length(input$trans2.1)), function(i){paste("I(",input$trans2.1[[i]],"^",input$trans1,")",sep="")}))
  #     predictors<-c(input$g.expo1.1,factor.pred, numeric.pred, other.numeric)
  #   }
  #   g.model<-(as.formula(paste(input$g.out1.1,"~",paste(predictors,collapse="+"),sep="")))
  #   g.function(g.dat(), input$g.expo1.1, input$g.out1.1, g.model , input$modelt, input$boot1)
  # })
  
  g.form<-eventReactive(input$go2,{
    if(input$g.cut1==""){
      factor.pred<-c(unlist(input$g.cov1.1))
      numeric.pred<-c(unlist(input$g.cov3.1))
      if(is.null(input$trans2.1)){
        predictors<-c(input$g.expo1.1,factor.pred, numeric.pred)
      }
      if(!is.null(input$trans2.1)){
        other.numeric<-c(sapply(1:(length(input$trans2.1)), function(i){paste("I(",input$trans2.1[[i]],"^",input$trans1,")",sep="")}))
        predictors<-c(input$g.expo1.1,factor.pred, numeric.pred, other.numeric)
      }
      g.model<-(as.formula(paste(input$g.out1.1,"~",paste(predictors,collapse="+"),sep="")))
      g.function(g.dat(), input$g.expo1.1, input$g.out1.1, g.model , input$modelt, input$boot1)
    }
    if(input$g.cut1!=""){
      dat<-g.dat()
      dat$intervention4gformula<-ifelse(dat[,input$g.expo1.1]>=as.numeric(input$g.cut1),1,0)
      factor.pred<-c(unlist(input$g.cov1.1))
      numeric.pred<-c(unlist(input$g.cov3.1))
      if(is.null(input$trans2.1)){
        predictors<-c(dat$intervention4gformula,factor.pred, numeric.pred)
      }
      if(!is.null(input$trans2.1)){
        other.numeric<-c(sapply(1:(length(input$trans2.1)), function(i){paste("I(",input$trans2.1[[i]],"^",input$trans1,")",sep="")}))
        predictors<-c(dat$intervention4gformula,factor.pred, numeric.pred, other.numeric)
      }
      g.model<-(as.formula(paste(input$g.out1.1,"~",paste(predictors,collapse="+"),sep="")))
      g.function(dat, dat$intervention4gformula, input$g.out1.1, g.model , input$modelt, input$boot1)
    }
  })
  
  
  
  output$g.results1<-renderPrint({
    g.form()
  })
  
  output$model1<-renderPrint({
    factor.pred<-c(unlist(input$g.cov1.1))
    numeric.pred<-c(unlist(input$g.cov3.1))
    if(is.null(input$trans2.1)){
      predictors<-c(input$g.expo1.1,factor.pred, numeric.pred)
    }
    if(!is.null(input$trans2.1)){
      other.numeric<-c(sapply(1:(length(input$trans2.1)), function(i){paste("I(",input$trans2.1[[i]],"^",input$trans1,")",sep="")}))
      predictors<-c(input$g.expo1.1,factor.pred, numeric.pred, other.numeric)
    }
    print(paste(input$g.out1.1,"~",paste(predictors,collapse="+"),sep=""))
  })
  
})

# Run the application 
shinyApp(ui = ui, server = server)
