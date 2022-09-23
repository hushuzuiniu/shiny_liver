## app.R ##
# library the packages we need
library('rsconnect')
library('shinydashboard')
library('shinydashboardPlus')
library('edgeR')
library('pheatmap')
library('ggplot2')
library('PCAtools')
library('GSEABase')
library('DT')
library('stringr')
library('dashboardthemes')
library("markdown")

# read the data that we need
counts<-read.csv("cpm.csv",header = TRUE,row.names = 1)
sample_detail <- read.csv("sample_info.csv",header = TRUE,row.names = 1)
pca_output <- readRDS("pca.rds")

# load the .gmt files for DE set analysis
# kegg_data <- getGmt('c2.cp.kegg.v7.5.1.symbols.gmt')
gomf_data <- getGmt('go.mf.18022019.symbols.gmt')
gobp_data <- getGmt('go.bp.18022019.symbols.gmt')
gocc_data <- getGmt('go.cc.18022019.symbols.gmt')
reactome_data <- getGmt('reactome.05092018.symbols.gmt')

# KEGG_Pathways<-geneIds(kegg_data)
GO_Molecular_Functions<-geneIds(gomf_data)
GO_Biological_Processes<-geneIds(gobp_data)
GO_Cellular_Components<-geneIds(gocc_data)
Reactome<-geneIds(reactome_data)
pathway_name <- c('GO_Biological_Processes',
                  'GO_Cellular_Components',
                  'GO_Molecular_Functions',
                  'Reactome')

# read the DE results and DE set results
file_list <-list.files()
file_name<-str_match(file_list,'(.*?)\\.(GOBP|GOCC|GOMF|KEGG|Reactome|results)\\.')[,2]
unique_name<-unique(file_name)
unique_name<-na.omit(unique_name)
# results <- GOBPs <- GOCCs <- GOMFs <- KEGGs <- Reactomes <- list()
results <- GOBPs <- GOCCs <- GOMFs <- Reactomes <- list()
for (i in 1:length(unique_name)){
  results[[i]] <- read.csv(paste0(unique_name[i],'.results.csv'),header=TRUE)
  GOBPs[[i]] <- read.csv(paste0(unique_name[i],'.GOBP.csv'),header=TRUE,row.names = 1)
  GOCCs[[i]] <- read.csv(paste0(unique_name[i],'.GOCC.csv'),header=TRUE,row.names = 1)
  GOMFs[[i]] <- read.csv(paste0(unique_name[i],'.GOMF.csv'),header=TRUE,row.names = 1)
  Reactomes[[i]] <- read.csv(paste0(unique_name[i],'.Reactome.csv'),header=TRUE,row.names = 1)
  # KEGGs[[i]] <- read.csv(paste0(unique_name[i],'.KEGG.csv'),header=TRUE,row.names = 1)
  
  results[[i]]<-results[[i]][!duplicated(results[[i]]$gene_name), ]
  rownames(results[[i]])<- results[[i]]$gene_name
  results[[i]]<-results[[i]][,-1]
  results[[i]]<-results[[i]][order(results[[i]]$FDR),]
}
# unique_name<-c()
# names(results) <- names(GOBPs) <- names(GOCCs) <- names(GOMFs) <- names(KEGGs) <- names(Reactomes) <- unique_name
names(results) <- names(GOBPs) <- names(GOCCs) <- names(GOMFs) <- names(Reactomes) <- unique_name
select_choices <- names(sample_detail)[c(1:length(unique_name),length(sample_detail))]
#########################################################################################################################


####################################
#         User interface           #
####################################

ui <- dashboardPage(                        
  dashboardHeader(                   #header of the app
    title = span(tagList(icon("react"),
                         "Differential gene & gene set expression analysis")),
    titleWidth = 470),
  dashboardSidebar( 
    width = 300,
    # side menu and its items
    sidebarMenu(
      
      menuItem("Gene Differential Expression Analysis" , 
               tabname = "menu1", 
               icon = icon("table"),
               startExpanded = TRUE,
               menuSubItem("Gene Differential Expression Analysis Intro",
                           tabName ="DE_Analysis",
                           icon = icon("th")),
               menuSubItem("PCA Plot", 
                           tabName = "PCA", 
                           icon = icon("eye")),
               menuSubItem("Differential Expression Gene Table",
                           tabName = "Gene-table", 
                           icon = icon("map-signs")),
               menuSubItem("Single Gene Histogram",
                           tabName = "Barplot", 
                           icon = icon("share-alt")),
               menuSubItem("Single Gene Boxplot",
                           tabName = "box-Plot",
                           icon = icon("chart-bar")),
               menuSubItem("Heatmap", 
                           tabName = "Heatmap", 
                           icon = icon("dashboard")),
               menuSubItem("MA Plot",
                           tabName = "MA-Plot",
                           icon = icon("sync")),
               menuSubItem("Volcano Plot",
                           tabName = "Volcano-Plot",
                           icon = icon("filter"))
      ),
      menuItem("Gene Set Enrichment Analysis" , 
               tabname = "menu2", 
               icon = icon("table"),
               startExpanded = TRUE,
               menuSubItem("Gene Set Enrichment Analysis Intro",
                           tabName = "DE_set_analysis",
                           icon = icon("th")),
               menuSubItem("Enrichment Table", 
                           tabName = "pathway", 
                           icon = icon("map-signs")),
               menuSubItem("Gene Set Table", 
                           tabName = "gene", 
                           icon = icon("share-alt")),
               menuSubItem("DE Set Heatmap",
                           tabName = "Heatmap-new", 
                           icon = icon("chart-line"))
      ),
      selectInput("sel.group", label="Sample Contrast Group", 
                  choices = unique_name[1:length(unique_name)],
                  selected = NULL)
    )
  ),
  dashboardBody(    # content of the app 
    shinyDashboardThemes(theme='onenote'),
    # customTheme,    # theme of the app
    tags$style(
      '
        @media (min-width: 768px){
          .sidebar-mini.sidebar-collapse .main-header .logo {
              width: 470px; 
          }
          .sidebar-mini.sidebar-collapse .main-header .navbar {
              margin-left: 470px;
          }
        }
        '
    ),
    tabItems(
      # DE intro content
      tabItem(tabName = "DE_Analysis",
              fluidRow(
                includeMarkdown("DE_analysis.Rmd"),
                width = 12),
      ),
      # First tab content
      tabItem(tabName = "PCA",
              fluidRow(
                box(title = "Principal Component Analysis Plot",
                    plotOutput("plot.PCA"),width = 12,height = 300),
              ),
              fluidRow(
                box(
                  title = "X axis",
                  width = 4,
                  selectInput("sel1", label=NULL, 
                              choices = paste("PC",1:10,sep = ""),
                              selected = "PC1")
                ),
                box(
                  title = "Y axis",
                  width = 4,
                  selectInput("sel2", label=NULL,
                              choices = paste("PC",1:10,sep = ""),
                              selected = "PC2")
                ),
                box(
                  title = "Point Size",
                  width = 4,
                  sliderInput("slider1", label=NULL, 0.5, 10, 5)
                )),
              fluidRow(
                box(
                  title = "Legend Group",
                  width = 6,
                  selectInput("sel15",label = NULL,
                              choices = select_choices)
                  
                ),
                box(
                  title = "Add labels to points",
                  width = 4,
                  radioButtons("radio1", label=NULL, 
                               choices = c("YES", "NO"),
                               selected = "YES")
                ),
                
                downloadButton('downloadPlot1', 'Download Plot')
              )
      ),
      
      # second single gene expression barplot tab content
      tabItem(tabName = "Barplot",
              fluidRow(
                box(
                  title = "First please enter a gene name you want to view the differential expression histogram of the gene",
                  width = 8,
                  textInput("sel5", 
                            label=NULL,
                            placeholder = "please enter a gene name(eg:SPATA21)")
                ),
                box(
                  title = "Show the grouped histogram",
                  width = 4,
                  radioButtons("radio12", label=NULL, 
                               choices = c("TRUE", "FALSE"),
                               selected = "FALSE")
                  
                ),
                fluidRow(
                  box(title = "Individual Gene Expression Histogram",
                      plotOutput("plot.bar"),width = 12,height = 300)
                  
                ),
                downloadButton('downloadPlot2', 'Download Plot'),
                
                fluidRow(
                  box(title = "Differential Expression Top Genes",
                      DT::dataTableOutput("DE_top_gene"),width = 12,height = 300)
                ),
                
              )
      ),
      
      # third top gene table
      tabItem(tabName = "Gene-table",
              fluidRow(
                box(title = "Gene Differential Expression Results Table",
                    DT::dataTableOutput("table_gene"),width = 12,height = 300)
              ),
              downloadButton('downloadData1', 'Download Data')
      ),
      
      # fourth  heatmap tab content
      tabItem(tabName = "Heatmap",
              
              fluidRow(
                box(title = "DE Top Genes Heatmap",
                    plotOutput("plot.heatmap"),
                    width = 12,height = 300),
              ),
              
              fluidRow(
                box(
                  title = "Sample Annotation", width = 8,
                  selectInput("sel6", label=NULL, 
                              choices = select_choices, 
                              selected = NULL,
                              multiple = T)
                ),
                box(
                  title = "Annotation for column", width = 4,
                  radioButtons("radio2", label=NULL, 
                               choices = c("TRUE", "FALSE"),
                               selected = "FALSE")
                )
              ),
              fluidRow( 
                box(
                  title = "Number of top DE genes", width = 4,
                  sliderInput("slider4", label=NULL, 10, 100, 15)
                ),
                box(
                  title = "Show sample names", width = 4,
                  radioButtons("radio3", label=NULL, 
                               choices = c("TRUE", "FALSE"),
                               selected = "FALSE")
                ),
                box(
                  title = "Show gene names", width = 4,
                  radioButtons("radio10", label=NULL, 
                               choices = c("TRUE", "FALSE"),
                               selected = "TRUE")
                ),
                
              ),
              fluidRow(
                
                
                box(
                  title = "Cluster genes", width = 4,
                  radioButtons("radio4", label=NULL, 
                               choices = c("TRUE", "FALSE"),
                               selected = "TRUE")
                ),
                box(
                  title = "Cluster samples", width = 4,
                  radioButtons("radio5", label=NULL, 
                               choices = c("TRUE", "FALSE"),
                               selected = "TRUE")
                ),
                box(
                  title = "Font size of gene names",width = 4,
                  sliderInput("slider13",label = NULL,2,20,10)
                  
                ),
                downloadButton('downloadPlot3', 'Download Plot')
              )
      ),
      
      # fifth box plot tab content
      tabItem(tabName = "box-Plot",
              fluidRow(
                box(
                  title = "First please enter a gene name you want to view the differential expression boxplot of the gene",
                  width = 12,
                  textInput("sel16", 
                            label=NULL,
                            placeholder = "please enter a gene name(eg:SPATA21)")
                ),
                
                
                fluidRow(
                  box(title = "Individual Gene Expression Boxplot",
                      plotOutput("plot.box"),width = 12,height = 300)
                ),
                downloadButton('downloadPlot4', 'Download Plot'),
                
                fluidRow(
                  box(title = "Differential Expression Top Genes",
                      DT::dataTableOutput("DE_top_gene1"),width = 12,height = 300)
                ),
              )
      ),
      # fifth MA plot tab content
      tabItem(tabName = "MA-Plot",
              fluidRow(
                box(title = "M-versus-A Plot",
                    plotOutput("plot.MA"),width = 12,height = 300),
              ),
              fluidRow(
                box(
                  title = "FDR threshold", width = 6,
                  sliderInput("slider5", label=NULL, 10^-10, 0.05, 0.05)
                ),
                box(
                  title = "FC threshold", width = 6,
                  sliderInput("slider6", label=NULL, 0.5, 5, 2)
                ),
                downloadButton('downloadPlot7', 'Download Plot')
              )
      ),
      
      
      # sixth volcano plot tab content
      tabItem(tabName = "Volcano-Plot",
              fluidRow(
                box(title = "Volcano Plot",
                    plotOutput("plot.vol"),width = 12,height = 300),
              ),
              fluidRow(
                box(
                  title = "log2 fold change", width = 4,
                  sliderInput("slider7", label=NULL,  0.1, 4, 0.6)
                ),
                box(
                  title = "-log10(FDR)", width = 4,
                  sliderInput("slider8", label=NULL, 10^-10, 0.1, 0.05)
                ),
                box(
                  title = "labeled points abs(logFC)", width = 4,
                  sliderInput("slider9", label=NULL, 2, 12, 4)
                ),
                downloadButton('downloadPlot5', 'Download Plot')
              )
      ),
      # DE set intro tab content
      tabItem(tabName = "DE_set_analysis",
              fluidRow(
                includeMarkdown("DE_set_analysis.Rmd"),
                width = 12),
      ),
      # seventh DE set heatmap tab content
      tabItem(tabName = "Heatmap-new",
              fluidRow(
                box(
                  title = "Gene Set Heatmap",
                  plotOutput("plot.heatmapnew"),width = 12,height = 300),
              ),
              fluidRow(
                box(
                  title = "Choose your enrichment", width = 4,
                  selectInput("sel7", label=NULL, 
                              choices = pathway_name)
                ),
                box(
                  title = "Gene set selection", width = 8,
                  selectInput("sel8", label=NULL, choices = NULL)
                ),
                
                fluidRow(
                  box(
                    title = "Number of gene set ", width = 4,
                    sliderInput("slider10", label=NULL, 2, 10, 2)
                  ),
                  box(
                    title = "Sample Annotation", width = 8,
                    selectInput("sel13", label=NULL, 
                                choices = select_choices, 
                                selected = NULL,
                                multiple = T)
                  )
                )
              ),
              fluidRow(
                box(
                  title = "Annotation for column", width = 4,
                  radioButtons("radio9", label=NULL, 
                               choices = c("TRUE", "FALSE"),
                               selected = "FALSE")
                ),
                box(
                  title = "Show sample names", width = 4,
                  radioButtons("radio6", label=NULL, 
                               choices = c("TRUE", "FALSE"),
                               selected = "FALSE")
                ),
                box(
                  title = "Show gene names", width = 4,
                  radioButtons("radio11", label=NULL, 
                               choices = c("TRUE", "FALSE"),
                               selected = "TRUE")
                ),
              ),
              fluidRow(
                box(
                  title = "Cluster genes", width = 4,
                  radioButtons("radio7", label=NULL, 
                               choices = c("TRUE", "FALSE"),
                               selected = "TRUE")
                ),
                
                box(
                  title = "Cluster samples", width = 4,
                  radioButtons("radio8", label=NULL, 
                               choices = c("TRUE", "FALSE"),
                               selected = "TRUE")
                ),
                box(
                  title = "Font size of gene names",width = 4,
                  sliderInput("slider14",label = NULL,2,20,10)
                  
                ),
                
                downloadButton('downloadPlot6', 'Download Plot')
              ),
      ),
      
      # eighth pathway table tab content
      tabItem(tabName = "pathway",
              fluidRow(
                box(
                  title = "Choose your enrichment", width = 12,
                  selectInput("sel10", label=NULL, 
                              choices = pathway_name)
                )
              ),
              fluidRow(
                box(title = "Pathway Table", width =12,
                    DT::dataTableOutput("pathwaytable")
                ),
                downloadButton('downloadData2', 'Download Data')
              )
              
      ),
      
      # ninth set genes table tab content
      tabItem(tabName = "gene",
              fluidRow(
                box(
                  title = "Choose your enrichment", width = 4,
                  selectInput("sel11", label=NULL, 
                              choices = pathway_name)
                ),
                box(
                  title = "Gene set selection", width = 8,
                  selectInput("sel12", label=NULL, choices = NULL)
                )
              ),
              fluidRow(
                box(title = "Gene Set Table", width =12,
                    DT::dataTableOutput("genetable")
                ),
                downloadButton('downloadData3', 'Download Data')
              )
      )
      
    )
  )
)
####################################
#         Server interface         #
####################################
server <- function(input, output,session) {
  
  # plot the PCA server
  plotInput1 <- reactive({
    observe({
      updateSelectInput(session, "sel1", 
                        choices = pca_output$components, 
                        selected = input$sel1)
    })
    observe({
      updateSelectInput(session, "sel2", 
                        choices = pca_output$components, 
                        selected = input$sel2)
    })
    
    
    # temp<-pca_output$metadata[,input$sel15]
    
    if(input$radio1 == "YES"){
      biplot(pca_output,
             pointSize = input$slider1,
             x=input$sel1,
             y=input$sel2,
             lab= rownames(pca_output$metadata),
             colby = input$sel15,
             legendPosition = 'right',
             legendLabSize = 8,
             legendIconSize = 3,
             legendTitleSize  =  10)
    }
    else{
      biplot(pca_output,
             pointSize = input$slider1,
             x=input$sel1,
             y=input$sel2, 
             lab = NULL,
             colby = input$sel15,
             legendPosition = 'right',
             legendLabSize = 8,
             legendIconSize = 3,
             legendTitleSize  =  10)
    }
  })
  
  output$plot.PCA <- renderPlot({
    plotInput1()
  })
  output$downloadPlot1 <- downloadHandler(
    filename = function() { paste('PCA plot', '.png', sep='') },
    content = function(file) {
      ggsave(width=8,height=4,dpi = 800,file,plotInput1())
    }
  )
  
  # get the top genes table
  tableInput1<- reactive({
    top_100_genes<- results[[input$sel.group]]
    
  }
  )
  output$table_gene <- DT::renderDataTable({
    tableInput1()
    
  })
  output$downloadData1 <- downloadHandler(
    filename = function() { paste(input$sel.group,' top genes table', '.csv', sep='') },
    content = function(file) {
      write.csv(tableInput1(), file)
    }
  )
  
  # plot individual gene differential expression histogram
  plotInput2 <- reactive({
    output$DE_top_gene <- DT::renderDataTable({
      DE_top_gene<- results[[input$sel.group]][1:100,c(1,4,5)]
      DE_top_gene
    })
    gene_name<-rownames(results[[input$sel.group]])
    if(input$radio12 == FALSE){
      if(input$sel5 %in% gene_name){
        if(input$sel.group== "Histological diagnosis: yes V no"){
          expr <- cpm(counts)
          plot_data <- cbind(sample_detail, expression = expr[input$sel5,])
          index <- which(unique_name == input$sel.group)
          plot_data <- plot_data[order(plot_data[,index]),] 
          plot_data$histological_diagnosis_of_nash<-as.factor(plot_data$histological_diagnosis_of_nash)
          plot_data$test <- row.names(plot_data)
          ggplot(plot_data, aes(x=test, y=expression))+ 
            geom_col(aes(fill=histological_diagnosis_of_nash))+
            scale_x_discrete(limit = plot_data$test)+
            theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
            theme(legend.key.size = unit(0.1,"inches"))+
            labs(x = input$sel.group, y = 'CPM',fill = names(plot_data)[index])+
            ggtitle(input$sel5)
          
        }else{
          
          expr <- cpm(counts)
          plot_data <- cbind(sample_detail, expression = expr[input$sel5,])
          index <- which(unique_name == input$sel.group)
          plot_data <- plot_data[order(plot_data[,index]),] 
          ggplot(plot_data, aes(x=reorder(rownames(plot_data),plot_data[,index]), 
                                y=expression, fill=factor(plot_data[,index]))) + 
            geom_bar(stat = 'identity') +  theme_bw() +
            theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
            theme(legend.key.size = unit(0.1,"inches"))+
            labs(x = input$sel.group, y = 'CPM', fill = names(plot_data)[index])+
            ggtitle(input$sel5)
        }
      }else{
        return("please enter a gene name")
      }
    }else{
      if(input$sel5 %in% gene_name){
        expr <- cpm(counts)
        plot_data <- cbind(sample_detail, expression = expr[input$sel5,])
        index <- which(unique_name == input$sel.group)
        plot_data <- plot_data[order(plot_data[,index]),] 
        ggplot(plot_data, aes(x=factor(plot_data[,index]), y=expression, fill=factor(plot_data[,index]))) + 
          geom_bar(stat = 'identity') +  theme_bw() +
          theme(axis.text.x = element_text(angle = 0, hjust = 0)) + 
          labs(x = input$sel.group, y = 'CPM', fill = names(plot_data)[index])+
          ggtitle(input$sel5)
      }else
        return("please enter a gene name")
    }
    
  })
  output$plot.bar <- renderPlot({
    plotInput2()
  })
  output$downloadPlot2 <- downloadHandler(
    filename = function() { paste(input$sel5,' gene expression','.png', sep='') },
    content = function(file) {
      ggsave(width=8,height=4,dpi = 800,file,plotInput2())
    }
  )
  
  
  # plot the boxplot server
  plotInput4 <- reactive({
    output$DE_top_gene1 <- DT::renderDataTable({
      DE_top_gene1<- results[[input$sel.group]][1:100,c(1,4,5)]
      DE_top_gene1
    })
    gene_name<-rownames(results[[input$sel.group]])
    if(input$sel16 %in% gene_name){
      expr <- cpm(counts)
      box_data <- cbind(sample_detail, expression = expr[input$sel16,])
      index <- which(unique_name == input$sel.group)
      box_data <- box_data[order(box_data[,index]),] 
      ggplot(box_data, aes(x=as.factor(box_data[,index]), y=expression))+
        geom_boxplot(aes(fill=as.character(box_data[,index]))) +
        labs(x = input$sel.group, y = 'CPM', fill = names(box_data)[index])+
        ggtitle(input$sel16)+
        theme_bw() 
    }else{
      return("please enter a gene name")
    }
    
  })
  
  output$plot.box <- renderPlot({
    plotInput4()
  })
  
  output$downloadPlot4 <- downloadHandler(
    filename = function() { paste(input$sel.group,' box plot','.png', sep='') },
    content = function(file) {
      ggsave(width=8,height=4,dpi = 800,file,plotInput4())
    }
  )
  
  
  # plot heatmap server
  plotInput3 <- reactive({
    top_genes<- results[[input$sel.group]][1:input$slider4,]
    plotmatrix <- cpm(counts,log=TRUE)[rownames(top_genes),]
    
    annotation_col <- dplyr::select(sample_detail, input$sel6)
    
    if(input$radio2 == "TRUE"){
      annotation_col<- annotation_col[order(annotation_col[,1]),,drop=FALSE]
      plotmatrix<-plotmatrix[,rownames(annotation_col)]
      pheatmap(
        plotmatrix, 
        show_rownames = as.logical(input$radio10),
        show_colnames = as.logical(input$radio3),
        border_color = NA,
        legend = TRUE, 
        cluster_cols = as.logical(input$radio5),
        cluster_rows = as.logical(input$radio4),
        scale = 'row',
        fontsize_row = input$slider13,
        annotation_names_col = FALSE,
        annotation_col = annotation_col,
        
      )
    }else{
      pheatmap(
        plotmatrix, 
        show_rownames = as.logical(input$radio10),
        show_colnames = as.logical(input$radio3),
        border_color = NA,
        legend = TRUE, 
        cluster_cols = as.logical(input$radio5),
        cluster_rows = as.logical(input$radio4),
        scale = 'row',
        fontsize_row = input$slider13,
        annotation_names_col = FALSE,
        annotation_col = NULL
      )
    }
  })
  output$plot.heatmap <- renderPlot({
    plotInput3()
  })
  
  output$downloadPlot3 <- downloadHandler(
    filename = function() { paste(input$sel.group,' heatmap','.png', sep='') },
    content = function(file) {
      ggsave(width=8,height=4,dpi = 800,file,plotInput3())
    }
  )
  # plot the MA server
  plotInput7 <- reactive({
    
    result_file <- results[[input$sel.group]]
    result_file$significant <- 'no'
    result_file$significant[result_file$FDR <= input$slider5  ] <- 'yes'
    ggplot(result_file, aes(logCPM, logFC,  color=significant )) + 
      geom_point(alpha = 0.2) +   
      scale_colour_manual(name = 'significant', 
                          values = setNames(c('red','grey'),c('yes', 'no'))) +
      geom_hline(yintercept=log2(input$slider6), linetype= "dashed") +
      geom_hline(yintercept=-1*log2(input$slider6), linetype= "dashed") +
      theme_bw()
  })
  
  output$plot.MA <- renderPlot({
    plotInput7()
  })
  
  output$downloadPlot7 <- downloadHandler(
    filename = function() { paste(input$sel.group,' MA plot','.png', sep='') },
    content = function(file) {
      ggsave(file,plotInput7())
    }
  )
  
  
  # plot the volcano server
  plotInput5 <- reactive({
    
    genes<-results[[input$sel.group]]
    genes$diffexpressed <- "NO"
    genes$diffexpressed[genes$logFC > input$slider7 & genes$FDR < input$slider8] <- "UP"
    genes$diffexpressed[genes$logFC < -input$slider7 & genes$FDR < input$slider8] <- "DOWN"
    ggplot(data=genes, aes(logFC, -log10(FDR))) +
      geom_point(aes(col=diffexpressed),)+
      geom_text_repel(aes(label=ifelse(abs(logFC)>=input$slider9, 
                                       as.character(row.names(genes)),'')))+
      theme_classic()+
      geom_vline(xintercept=c(input$slider7,-input$slider7), linetype = 'dashed')
    
  })
  
  output$plot.vol <- renderPlot({
    plotInput5()
  })
  
  output$downloadPlot5 <- downloadHandler(
    filename = function() { paste(input$sel.group,' volcano plot','.png', sep='') },
    content = function(file) {
      ggsave(width=8,height=4,dpi = 800,file,plotInput5())
    }
  )
  
  # plot DE set heatmap server
  plotInput6 <- reactive({
    
    annotation_col <- dplyr::select(sample_detail, input$sel13)
    if(input$sel7 == "GO_Biological_Processes"){
      
      observe({
        updateSelectInput(session, "sel8", choices = row.names(GOBPs[[input$sel.group]]), selected = input$sel8)
      })
      
      genes<-GO_Biological_Processes[[input$sel8]]
      genes<-genes[which(genes %in% rownames(results[[input$sel.group]]))]
      
    }else if(input$sel7 == "GO_Cellular_Components"){
      
      observe({
        updateSelectInput(session, "sel8", choices = row.names(GOCCs[[input$sel.group]]), selected = input$sel8)
      })
      
      genes<-GO_Cellular_Components[[input$sel8]]
      genes<-genes[which(genes %in% rownames(results[[input$sel.group]]))]
      
    }else if(input$sel7 == "GO_Molecular_Functions"){
      
      
      observe({
        updateSelectInput(session, "sel8", choices = row.names(GOMFs[[input$sel.group]]), selected = input$sel8)
      })
      
      genes<-GO_Molecular_Functions[[input$sel8]]
      genes<-genes[which(genes %in% rownames(results[[input$sel.group]]))]
      
    }else if(input$sel7 == "KEGG_Pathways"){
      
      observe({
        updateSelectInput(session, "sel8", choices = row.names(KEGGs[[input$sel.group]]), selected = input$sel8)
      })
      
      genes<-KEGG_Pathways[[input$sel8]]
      genes<-genes[which(genes %in% rownames(results[[input$sel.group]]))]
      
    }else{
      
      observe({
        updateSelectInput(session, "sel8", choices = row.names(Reactomes[[input$sel.group]]), selected = input$sel8)
      })
      
      genes<-Reactome[[input$sel8]]
      genes<-genes[which(genes %in% rownames(results[[input$sel.group]]))]
    }
    
    plotmatrix.original <- cpm(counts,log=TRUE)[genes,]
    
    
    if(nrow(plotmatrix.original) < 2){
      plot.new()
    }
    else{
      if(input$radio9 == "TRUE"){
        observe({
          updateSliderInput(session, "slider10", min = 2, max = nrow(plotmatrix.original))
        })
        
        plotmatrix <- plotmatrix.original[1:input$slider10,]
        annotation_col<-annotation_col[order(annotation_col[,1]),,drop=FALSE]
        plotmatrix<-plotmatrix[,rownames(annotation_col)]
        pheatmap(
          plotmatrix, 
          show_rownames = as.logical(input$radio11),
          show_colnames = as.logical(input$radio6),
          border_color = NA,
          legend = TRUE, 
          cluster_cols = as.logical(input$radio8),
          cluster_rows = as.logical(input$radio7),
          scale = 'row',
          fontsize_row = input$slider14,
          annotation_names_col = FALSE,
          annotation_col = annotation_col
        ) }else{
          observe({
            updateSliderInput(session, "slider10", min = 2, max = nrow(plotmatrix.original))
          })
          
          plotmatrix <- plotmatrix.original[1:input$slider10,]
          
          pheatmap(
            plotmatrix, 
            show_rownames = as.logical(input$radio11),
            show_colnames = as.logical(input$radio6),
            border_color = NA,
            legend = TRUE, 
            cluster_cols = as.logical(input$radio8),
            cluster_rows = as.logical(input$radio7),
            fontsize_row = input$slider14,
            scale = 'row',
            annotation_names_col = FALSE,
            annotation_col = NULL)
          
        }
      
    }
    
  })
  output$plot.heatmapnew <- renderPlot({
    plotInput6()
  })
  
  output$downloadPlot6 <- downloadHandler(
    filename = function() { paste('DE set heatmap','.png', sep='') },
    content = function(file) {
      ggsave(width=8,height=4,dpi = 800,file,plotInput6())
    }
  )
  
  
  # pathway table
  tableInput2<-reactive({
    if(input$sel10 == "GO_Biological_Processes"){
      pathway_table <- GOBPs[[input$sel.group]]
      
    }else if(input$sel10 == "GO_Cellular_Components"){
      pathway_table <- GOCCs[[input$sel.group]]
      
    }else if(input$sel10 == "GO_Molecular_Functions"){
      pathway_table <- GOMFs[[input$sel.group]]
      
    }else if(input$sel10 == "KEGG_Pathways"){
      pathway_table <- KEGGs[[input$sel.group]]
      
    }else{
      pathway_table <- Reactomes[[input$sel.group]]
      
    }
    
    pathway_table
    
  })
  output$pathwaytable <- DT::renderDataTable({
    tableInput2()
    
  })
  
  output$downloadData2 <- downloadHandler(
    filename = function() { paste(input$sel10,' table', '.csv', sep='') },
    content = function(file) {
      write.csv(tableInput2(), file)
    }
  )
  
  # gene set table 
  tableInput3<-reactive({
    if(input$sel11 == "GO_Biological_Processes"){
      
      observe({
        updateSelectInput(session, "sel12", 
                          choices = row.names(GOBPs[[input$sel.group]]), 
                          selected = input$sel12)
      })
      
      genes<-GO_Biological_Processes[[input$sel12]]
      genes<-genes[which(genes %in% rownames(results[[input$sel.group]]))]
      
    }else if(input$sel11 == "GO_Cellular_Components"){
      
      observe({
        updateSelectInput(session, "sel12", 
                          choices = row.names(GOCCs[[input$sel.group]]), 
                          selected = input$sel12)
      })
      
      genes<-GO_Cellular_Components[[input$sel12]]
      genes<-genes[which(genes %in% rownames(results[[input$sel.group]]))]
      
    }else if(input$sel11 == "GO_Molecular_Functions"){
      
      
      observe({
        updateSelectInput(session, "sel12", 
                          choices = row.names(GOMFs[[input$sel.group]]), 
                          selected = input$sel12)
      })
      
      genes<-GO_Molecular_Functions[[input$sel12]]
      genes<-genes[which(genes %in% rownames(results[[input$sel.group]]))]
      
    }else if(input$sel11 == "KEGG_Pathways"){
      
      observe({
        updateSelectInput(session, "sel12", 
                          choices = row.names(KEGGs[[input$sel.group]]), 
                          selected = input$sel12)
      })
      
      genes<-KEGG_Pathways[[input$sel12]]
      genes<-genes[which(genes %in% rownames(results[[input$sel.group]]))]
      
    }else{
      
      observe({
        updateSelectInput(session, "sel12", 
                          choices = row.names(Reactomes[[input$sel.group]]), 
                          selected = input$sel12)
      })
      
      genes<-Reactomes[[input$sel12]]
      genes<-genes[which(genes %in% rownames(results[[input$sel.group]]))]
    }
    
    results[[input$sel.group]][genes,]
    
  })
  output$genetable <- DT::renderDataTable({
    
    tableInput3()
  })
  
  output$downloadData3 <- downloadHandler(
    filename = function() { paste('DE set genes table', '.csv', sep='') },
    content = function(file) {
      write.csv(tableInput3(), file)
    }
  )
}

####################################
# Create the shiny app             #
####################################

shinyApp(ui, server)