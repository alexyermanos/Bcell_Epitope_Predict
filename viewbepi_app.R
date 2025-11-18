required_packages <- c(
  "bio3d", "NGLVieweR", "shiny", "dplyr", "ggplot2", "readr", "scales", "stringr", "tools", "tidyr"
)

for (pkg in required_packages) {
  if(!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg, repos="https://cloud.r-project.org")
  library(pkg, character.only = TRUE)
}

get_residueindex_colors <- function(residues) {
  n <- length(residues)
  cols <- colorRampPalette(c('red','orange','yellow','green','blue'))(n)
  data.frame(resno = residues, color = cols)
}

num_bins <- 21

# --- UI ---
ui <- fluidPage(
  titlePanel("Protein Structure Viewer"),
  sidebarLayout(
    sidebarPanel(
      fileInput("pdb_file", "Upload PDB File", accept = c(".pdb")),
      fileInput("csv_file", "Upload CSV File", accept = c(".csv")),
      hr(),
      selectInput('struc_representation', 'Structure representation',
                  c('cartoon','ribbon','ball+stick'), selected='cartoon'),
      selectInput('struc_color_scheme', 'Base color scheme',
                  c('uniform','residueindex','bfactor'), selected='residueindex'),
      uiOutput("viewmode_ui"),  # <-- dynamic
      selectInput('type', 'Selection representation',
                  c('ball+stick','cartoon','ribbon','surface','spacefill'), selected='surface'),
      selectInput('viewmethod', 'Selection method', c('score','topn','threshold')),
      conditionalPanel(condition="input.viewmethod=='topn'",
                       sliderInput('topn_value','Top n%',0,100,25,5)),
      conditionalPanel(condition="input.viewmethod=='threshold'",
                       uiOutput('threshold_slider_ui')),
      checkboxInput('normalize','Normalize values', value=FALSE),
      actionButton('add','Add selection'),
      actionButton('remove','Remove selection'),
      actionButton('download_png','Download Structure as PNG')
    ),
    mainPanel(
      verbatimTextOutput('selection'),
      plotOutput('score_plot', height='250px'),
      NGLVieweROutput('structure', height='600px')
    )
  )
)

# --- SERVER ---
server <- function(input, output, session) {
  
  # Reactively read CSV
  data_reactive_file <- reactive({
    req(input$csv_file)
    read_csv(input$csv_file$datapath, show_col_types=FALSE) %>%
      drop_na(chain, resno)
  })
  
  # Reactive PDB file
  pdb_file_reactive <- reactive({
    req(input$pdb_file)
    input$pdb_file$datapath
  })
  
  # Dynamically generate numeric column options
  output$viewmode_ui <- renderUI({
    req(data_reactive_file())
    numeric_cols <- names(data_reactive_file())[sapply(data_reactive_file(), is.numeric)]
    selectInput('viewmode', 'Selection type', choices = numeric_cols)
  })
  
  # Reactive version of your 'data'
  data_reactive <- reactive({
    req(input$viewmode)
    d <- data_reactive_file()
    d$view <- d[[input$viewmode]]
    if (isTRUE(input$normalize))
      d$view <- scales::rescale(d$view, to=c(0,1))
    d
  })
  
  output$threshold_slider_ui <- renderUI({
    req(input$viewmode)
    vals <- data_reactive()$view
    sliderInput("threshold_value", "Threshold",
                min = min(vals, na.rm = TRUE),
                max = max(vals, na.rm = TRUE),
                value = mean(vals, na.rm = TRUE),
                step = 0.01)
  })
  
  output$score_plot <- renderPlot({
    req(input$viewmode)
    d <- data_reactive()
    
    if (input$viewmethod == "score") {
      d$bin <- cut(d$view, breaks = num_bins, labels = FALSE)
      palette <- colorRampPalette(c("blue", "white", "red"))(num_bins)
      d$color <- palette[d$bin]
    } else {
      d <- left_join(d, get_residueindex_colors(na.omit(d$resno)), by='resno')
    }
    
    p <- ggplot(d, aes(x=resno, y=view)) +
      geom_col(aes(fill=color), width=1) +
      scale_fill_identity() +
      labs(x="Residue Number", y=input$viewmode,
           title=paste("Per-residue", input$viewmode)) +
      theme_bw()
    
    if (input$viewmethod == "threshold")
      p <- p + geom_hline(yintercept = input$threshold_value, color="red", linetype="dashed")
    
    if (input$viewmethod == "topn") {
      cutoff <- quantile(d$view, probs = 1 - input$topn_value/100, na.rm = TRUE)
      p <- p + geom_hline(yintercept = cutoff, color="red", linetype="dashed")
    }
    p
  })
  
  output$structure <- renderNGLVieweR({
    req(pdb_file_reactive())
    d <- data_reactive()
    NGLVieweR(pdb_file_reactive()) %>%
      addRepresentation(
        type = input$struc_representation,
        param = list(
          name='base',
          colorScheme=input$struc_color_scheme,
          sele = paste0("/0 AND :", unique(na.omit(d$chain)))
        )
      ) %>%
      stageParameters(backgroundColor='white') %>%
      setQuality('medium') %>%
      setFocus(0)
  })
  
  proxy <- NGLVieweR_proxy("structure")
  
  observeEvent(input$remove, {
    for (i in 1:num_bins) proxy %>% removeSelection(paste0("sel",i))
    proxy %>% removeSelection("threshold") %>% removeSelection("topn")
  })
  
  observeEvent(input$add, {
    d <- data_reactive()
    if (input$viewmethod == "score") {
      d$bin <- cut(d$view, breaks = num_bins, labels = FALSE)
      d$color <- colorRampPalette(c("blue","white","red"))(num_bins)[d$bin]
      d$sele <- paste0(d$resno)
      color_groups <- d %>% group_by(color) %>%
        summarise(sele = paste0(sele, collapse=" "), .groups="drop")
      
      for (i in seq_len(nrow(color_groups))) {
        proxy <- proxy %>%
          addSelection(input$type,
                       param = list(
                         name = paste0("sel",i),
                         sele = paste0("/0 AND :", unique(na.omit(d$chain)), " AND (", color_groups$sele[i], ")"),
                         colorValue = color_groups$color[i]))
      }
    } else if (input$viewmethod == "threshold") {
      top_res <- d %>% filter(view >= input$threshold_value)
      proxy <- proxy %>%
        addSelection(input$type,
                     param = list(
                       name="threshold",
                       sele = paste0("/0 AND :", unique(na.omit(d$chain)),
                                     " AND (", paste0(top_res$resno, collapse=" "), ")"),
                       colorScheme="residueindex"))
    } else if (input$viewmethod == "topn") {
      cutoff <- quantile(d$view, probs = 1 - input$topn_value/100)
      top_res <- d %>% filter(view >= cutoff)
      proxy <- proxy %>%
        addSelection(input$type,
                     param = list(
                       name="topn",
                       sele = paste0("/0 AND :", unique(na.omit(d$chain)),
                                     " AND (", paste0(top_res$resno, collapse=" "), ")"),
                       colorScheme="residueindex"))
    }
  })
  
  observeEvent({
    input$topn_value
    input$threshold_value
    input$normalize
  }, {
    req(input$viewmode)
    data <- data_reactive()
    proxy <- NGLVieweR_proxy("structure")
    
    if (input$viewmethod == "topn") {
      req(!is.na(input$topn_value))
      cutoff <- quantile(data$view, probs = 1 - input$topn_value/100, na.rm = TRUE)
      top_res <- data %>% filter(view >= cutoff)
      sele_string <- paste0(top_res$resno, ":", top_res$chain, collapse = " or ")
      proxy %>% updateSelection("topn", sele = paste0("/0 AND :", unique(na.omit(data$chain)), " AND (", sele_string, ")"))
    } else if (input$viewmethod == "threshold") {
      req(!is.na(input$threshold_value))
      top_res <- data %>% filter(view >= input$threshold_value)
      sele_string <- paste0(top_res$resno, ":", top_res$chain, collapse = " or ")
      proxy %>% updateSelection("threshold", sele = paste0("/0 AND :", unique(na.omit(data$chain)), " AND (", sele_string, ")"))
    }
  })
  
  observeEvent({
    input$download_png
  },{
    req(pdb_file_reactive())
    
    filename <- paste0(
      tools::file_path_sans_ext(basename(input$pdb_file$name)),
      "_snapshot.png"
    )
    
    proxy <- NGLVieweR_proxy("structure")
    
    proxy %>%
      snapShot(fileName = filename,
               param = list(
                 antialias = TRUE,
                 trim = TRUE,
                 transparent = TRUE,
                 scale = 1
               )
      )
  })
}

# --- Run App ---
shinyApp(ui, server, options=list(launch.browser=TRUE))