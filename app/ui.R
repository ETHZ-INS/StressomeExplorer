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
      menuItem("Phosphopeptides", tabName="tab_phos_pep"),
      tags$div(id="lablink", tags$a(href="https://bohaceklab.ethz.ch/", "bohaceklab.ethz.ch"))
    )
  ),
  dashboardBody(
    tagList(
      tags$head(
        tags$style(".inlineInputs div { display:inline-block; }"),
	tags$style(HTML('
#lablink {
  position: absolute;
  bottom:2px;
  left:8px;
  font-size: 110%;
  font-weight: bold;
  text-decoration: none;
}
#lablink a:hover {
  text-decoration: underline;
}
    '))
      )
    ),
    tabItems(
      tabItem("tab_gene", column(8,tags$h3(textOutput("gene_name"))),
        column(4, selectizeInput("gene_input", "Select Gene", choices=c(), multiple=FALSE)),
        column(3, selectizeInput("select_plottype", "Type of Plot", choices=c("violin plot","box plot"), multiple=FALSE)),
        column(2, checkboxInput('select_logaxis','Logarithmic Axis', value=T)),
        column(2, checkboxInput('select_plotpoints','Plot Points', value=T)),
        column(5, htmlOutput("availability")),
        box(title = "Transcriptome data (press + to view)", width=8, withSpinner(plotOutput("gene_plot", height=1200, width = 600)),collapsible = T, collapsed = T),
        box(title = "Experimental designs", width = 4,
            plotOutput("EDTS", height = "250px"),
            plotOutput("FvM", height = "250px"),
            plotOutput("LvR", height = "200px"),
            p("adult (~2-3 months) female mice"),
            plotOutput("CMV", height = "200px"),
            plotOutput("TRAP", height = "300px"),
            p("adult male (2-3 months) C57BL/6 mice were used unless otherwise stated"),
            p("vHC = ventral hippocampus, dHC = dorsal hippocampus, AS = acute stress"),
            collapsible = T, collapsed = T)
    ),
      tabItem("tab_snrna",
              column(12,tags$h3(textOutput("gene_name_sn"))),
              box(title = "Experimental design",width = 12,
                  plotOutput("EDsnRNA", height = "110px"),
                  p("adult male (2-3 months) C57BL/6 mice"),
                  p("AS = acute stress"),
                  collapsible = T, collapsed = T),
              box(title="single-nucleus RNA-seq", width=12,
                              withSpinner(plotOutput("snrna_plot", width=600, height=800)))
    ),
    tabItem("tab_prot",
            column(12,tags$h3(textOutput("gene_name_prot"))),
            box(title = "Experimental design",width = 12,
                plotOutput("EDprot", height = "90px"),
                p("adult male (2-3 months) C57BL/6 mice"),
                p("vHC = ventral hippocampus, dHC = dorsal hippocampus, AS = acute stress, CA = Cornu Ammonis 1, DG = dentate gyrus"),
                collapsible = T, collapsed = T),
            box(title = "Proteomic datasets", width=12, withSpinner(plotOutput("protplot_plot", width = 600, height = 800)))
    ),
      tabItem("tab_phos",
              column(12,tags$h3(textOutput("gene_name_phos"))),
              box(title = "Experimental design",width = 12,
                  plotOutput("EDphos", height = "200px"),
                  p("adult male (2-3 months) C57BL/6 mice"),
                  p("vHC = ventral hippocampus, dHC = dorsal hippocampus, AS = acute stress"),
                  collapsible = T, collapsed = T),
              box(title = "Phosphoproteomic datasets", width=12,
                             selectizeInput("phos_assay", "Assay", choices=c("scaledSVA","log2FC"), multiple=FALSE),
                             withSpinner(plotOutput("phospho_plot", width = 1000)))
    ),
      tabItem("tab_phos_pep",
              column(12,tags$h3(textOutput("gene_name_phospep"))),
              box(title = "Experimental design",width = 12, 
                  plotOutput("EDphos2", height = "200px"),
                  p("adult male (2-3 months) C57BL/6 mice"),
                  p("vHC = ventral hippocampus, dHC = dorsal hippocampus, AS = acute stress"),
                  collapsible = T, collapsed = T),
              box(title = "Peptide Explorer (warning: may be slow if protein has many peptides)", width=12,
                  withSpinner(plotOutput("phospho_pep_plot", width = 1000)),
                  p("missing data points indicate that peptide was not detected in a sample"))
    ))
  )
))
