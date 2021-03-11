library(shiny)
library(DT)
library(plotly)
library(shinydashboard)
library(shinycssloaders)

shinyUI( dashboardPage(

  dashboardHeader(title="stressome"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Select Gene", tabName="tab_gene"),
      menuItem("Single-nuclei", tabName="tab_snrna"),
      menuItem("Proteome", tabName="tab_prot"),
      menuItem("Phosphoproteome", tabName="tab_phos"),
      menuItem("Phosphopeptides", tabName="tab_phos_pep")
    )
  ),
  dashboardBody(
    tagList(
      tags$head(
        tags$link(rel="stylesheet", type="text/css", href=system.file("extdata/style.css", package="stressome")),
        tags$style(".inlineInputs div { display:inline-block; }")
      )
    ),
    tabItems(
      tabItem("tab_gene", column(8,tags$h3(textOutput("gene_name"))),
        column(4, selectizeInput("gene_input", "Select Gene", choices=c(), multiple=FALSE)),
        column(3, selectizeInput("select_plottype", "Type of Plot", choices=c("violin plot","box plot"), multiple=FALSE)),
        column(2, checkboxInput('select_logaxis','Logarithmic Axis', value=T)),
        column(2, checkboxInput('select_plotpoints','Plot Points', value=T)),
        column(5, htmlOutput("availability")),
        box(title = "Transcriptome data (press + to view)", width=12, withSpinner(plotOutput("gene_plot", height=1200, width = 600)),collapsible = T, collapsed = T)
    ),
      tabItem("tab_snrna",box(title="single-nucleus RNA-seq", width=12,
                              withSpinner(plotOutput("snrna_plot", width=600, height=500)))
    ),
      tabItem("tab_prot",box(title = "Proteomic datasets", width=12, withSpinner(plotOutput("protplot_plot", width = 600, height = 800)))
    ),
      tabItem("tab_phos",box(title = "Phosphoproteomic datasets", width=12,
                             selectizeInput("phos_assay", "Assay", choices=c("scaledSVA","log2FC"), multiple=FALSE),
                             withSpinner(plotOutput("phospho_plot", width = 1000)))
    ),
      tabItem("tab_phos_pep",box(title = "Peptide Explorer (warning: may be slow if protein has many peptides)", width=12,withSpinner(plotOutput("phospho_pep_plot", width = 1000)))
    ))
  )
))
