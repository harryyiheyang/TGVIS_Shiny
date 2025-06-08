library(shiny)
library(data.table)
library(dplyr)
library(shinyjs)
library(glue)

# Fix function conflicts
filter <- dplyr::filter
select <- dplyr::select
mutate <- dplyr::mutate

options(shiny.timeout = 300)  

# --- Section: Load data ---
trait_info <- readRDS("Trait.rds") %>% mutate(LinkName = tolower(ShortName))
w_list <- vector("list", nrow(trait_info))
for(i in 1:nrow(trait_info)){
  file_path <- glue("{trait_info$ShortName[i]}/FineMap.rds")
  if (file.exists(file_path)) {
    w_list[[i]] <- readRDS(file_path) %>%
      mutate(Outcome = trait_info$FullName[i])
  } else {
    w_list[[i]] <- NULL
  }
}
w_all <- bind_rows(w_list[!sapply(w_list, is.null)])%>%dplyr::select(Ensembl=Gene,Symbol=GeneSymbol,Tissue,Outcome,Identifier,
                                                                     `Pratt index of credible set`=CS.Pratt,
                                                                     `PIP of credible set`=CS.PIP,xQTL,`Splicing information`=Splicing)

# --- Section: UI Layout ---
ui <- fluidPage(
  useShinyjs(),
  
  tags$head(
    tags$style(HTML("
      .magnify-container {
        position: relative;
      }
      .magnify-container img {
        width: 100%;
        border: 1px solid #ccc;
      }
      .magnify-lens {
        position: absolute;
        border: 1px solid #000;
        width: 100px;
        height: 100px;
        visibility: hidden;
        pointer-events: none;
        background-repeat: no-repeat;
        background-size: 200% 200%;
        z-index: 1000;
      }
    ")),
    tags$script(HTML("
      Shiny.addCustomMessageHandler('enableMagnify', function(x) {
        function enableMagnify(el) {
          if (el.querySelector('.magnify-lens')) return;
          const lens = document.createElement('div');
          lens.classList.add('magnify-lens');
          el.appendChild(lens);

          const img = el.querySelector('img');
          img.addEventListener('mousemove', moveLens);
          lens.addEventListener('mousemove', moveLens);
          img.addEventListener('mouseenter', () => lens.style.visibility = 'visible');
          img.addEventListener('mouseleave', () => lens.style.visibility = 'hidden');

          function moveLens(e) {
            const posX = e.offsetX;
            const posY = e.offsetY;
            lens.style.left = (posX - lens.offsetWidth / 2) + 'px';
            lens.style.top = (posY - lens.offsetHeight / 2) + 'px';
            lens.style.backgroundImage = 'url(' + img.src + ')';
            lens.style.backgroundPosition = `-${posX}px -${posY}px`;
          }
        }

        document.querySelectorAll('.magnify-container').forEach(enableMagnify);
      });
    "))
  ),
  
  titlePanel("Trait Viewer"),
  
  tabsetPanel(
    tabPanel("Trait View",
             sidebarLayout(
               sidebarPanel(
                 selectizeInput(
                   inputId = "trait",
                   label = "Select or type a trait:",
                   choices = setNames(trait_info$ShortName, trait_info$FullName),
                   selected = NULL,
                   multiple = FALSE,
                   options = list(
                     placeholder = 'e.g. LDL or LDL Cholesterol',
                     create = FALSE
                   )
                 ),
                 conditionalPanel(
                   condition = "input.trait != ''",
                   br(),
                   downloadButton("download_finemap", "Download Trait Table")
                 )
               ),
               mainPanel(
                 uiOutput("main_display"),
                 tags$div(style = "margin-top: 10px;",
                          helpText("• You can locate a locus using Chromosome + BP, or search by gene name.",
                                   "• Clear the search box and press Enter to return to the summary table.")
                 )
               )
             )
    ),
    
    tabPanel("Gene View",
             sidebarLayout(
               sidebarPanel(
                 selectizeInput(
                   inputId = "gene_search",
                   label = "Select or type a gene:",
                   choices = NULL,
                   selected = NULL,
                   multiple = FALSE,
                   options = list(
                     placeholder = 'e.g. PCSK9 or ENSG00000169174',
                     create = FALSE,
                     maxOptions = 100
                   )
                 ),
                 helpText("• This view shows all traits where the queried gene is causal."),
                 br(),
                 downloadButton("download_gene_table", "Download Gene Table")
               ),
               mainPanel(
                 dataTableOutput("gene_table"),
                 textOutput("gene_msg")
               )
             )
    )
  ),
  
  # Fixed navigation buttons
  tags$div(
    id = "doc_button_area",
    actionButton("goto_doc", "Document", icon = icon("book")),
    style = "position: fixed; bottom: 20px; left: 20px; z-index: 1000;"
  ),
  
  tags$div(
    actionButton("goto_home", "Return to Home", icon = icon("home")),
    style = "position: fixed; bottom: 20px; left: 130px; z-index: 1000;"
  )
)

# --- Section: Server Logic ---
server <- function(input, output, session) {
  current_page <- reactiveVal("home")
  confirmed_trait <- reactiveVal(NULL)
  
  observe({
    if (current_page() == "home") {
      shinyjs::show("doc_button_area")
    } else {
      shinyjs::hide("doc_button_area")
    }
  })
  
  observeEvent(input$trait, {
    val <- input$trait
    if (val %in% trait_info$ShortName) {
      confirmed_trait(val)
    } else if (val %in% trait_info$FullName) {
      sn <- trait_info$ShortName[trait_info$FullName == val]
      updateSelectizeInput(session, "trait", selected = sn)
      confirmed_trait(sn)
    }
  }, ignoreInit = TRUE)
  
  observeEvent(input$goto_doc, {
    current_page("doc")
  })
  
  observeEvent(input$goto_home, {
    current_page("home")
    confirmed_trait(NULL)
    updateSelectizeInput(session, "trait", selected = "")
  })
  
  # Initialize gene choices
  observe({
    gene_symbols <- unique(w_all$Symbol)
    gene_ensembl <- unique(w_all$Ensembl)
    
    gene_symbols <- gene_symbols[!is.na(gene_symbols)]
    gene_ensembl <- gene_ensembl[!is.na(gene_ensembl)]
    
    gene_choices <- c(
      setNames(gene_symbols, gene_symbols),
      setNames(gene_ensembl, gene_ensembl)
    )
    
    gene_choices <- gene_choices[order(names(gene_choices))]
    
    updateSelectizeInput(
      session,
      "gene_search",
      choices = gene_choices,
      server = TRUE
    )
  })
  
  output$main_display <- renderUI({
    page <- current_page()
    
    if (page == "doc") {
      ui_doc <- tagList(
        h3("Document"),
        tags$ul(
          tags$li(HTML('eQTL and sQTL data are from GTEx v10: <a href="https://www.gtexportal.org/home/downloads/adult-gtex/qtl" target="_blank">GTEx download</a>')),
          tags$li("Tissues include 31 common body and brain tissues, such as fat, artery, liver, cortex, cerebellum, heart, etc."),
          tags$li(HTML('45 cardiometabolic traits are analyzed. All GWAS summary data are imputed via SBayesRC using a 7.3M SNP LD reference panel: <a href="https://github.com/zhilizheng/SBayesRC" target="_blank">SBayesRC GitHub</a>')),
          tags$li("Loci are defined as ±1MB windows centered on lead SNPs with P < 5×10⁻⁸ in trait GWAS."),
          tags$li("Instrument variants are selected with P < 1×10⁻⁵ in at least one gene-tissue pair, followed by Clumping+Thresholding (LD < 0.64) using UKBB (9,680 EUR) reference panel.")
        )
      )
      session$sendCustomMessage("enableMagnify", list())
      return(ui_doc)
    }
    
    trait <- confirmed_trait()
    if (is.null(trait) || trait == "") {
      ui_home <- tagList(
        h1("Tissue-Gene Pairs, Direct Causal Variants, and Infinitesimal Effects Selector"),
        div(class = "magnify-container",
            tags$img(src = "https://raw.githubusercontent.com/harryyiheyang/TGVIS/main/Overview.png",
                     width = "90%", style = "border: 1px solid #ccc; margin-bottom: 10px;"))
      )
      session$sendCustomMessage("enableMagnify", list())
      return(ui_home)
    }
    
    path <- file.path(trait, "FineMap.rds")
    if (!file.exists(path)) {
      return(tags$p("No fine-mapping results of causal gene-tissue pairs.", style = "color:red;"))
    }
    
    ui_search <- tagList(
      h4(HTML('Search by <em>Chromosome + Base Pair</em>, or <em>Gene</em>')),
      fluidRow(
        column(4, textInput("search_chr", "Chromosome:", placeholder = "e.g. 1–22")),
        column(4, textInput("search_pos", "Base Pair (GRCh38):", placeholder = "e.g. 53.8e6")),
        column(4, textInput("gene_query", "Search by Gene Symbol or Ensembl ID:", placeholder = "e.g. SORT1 or ENSG00000134243"))
      ),
      tags$p(
        style = "font-style: italic; color: gray;",
        "Enter chromosome + position or a gene to locate a locus. ",
        "Clear all fields and press Enter to return to the summary table."
      ),
      uiOutput("show_leading_summary"),
      uiOutput("locus_match_output")
    )
    
    session$sendCustomMessage("enableMagnify", list())
    return(ui_search)
  })
  
  output$locus_match_output <- renderUI({
    trait <- confirmed_trait()
    req(trait)
    
    path <- file.path(trait, "FineMap.rds")
    df <- readRDS(path)
    df$CHR <- as.integer(df$CHR)
    
    gene <- input$gene_query
    chr <- as.integer(input$search_chr)
    pos <- as.numeric(input$search_pos)
    
    if (!is.null(gene) && gene != "") {
      match_idx <- which(tolower(df$GeneSymbol) == tolower(gene) | tolower(df$Gene) == tolower(gene))
      if (length(match_idx) == 0) {
        return(tags$p("Target gene is not causal in current fine-mapping results.", style = "color: red; font-style: italic;"))
      }
      matched_row <- match_idx[1]
    } else {
      if (is.na(chr) || is.na(pos)) {
        return(tags$p("Please enter valid numeric values for both chromosome and position.",
                      style = "color: red; font-style: italic;"))
      }
      in_block <- df$CHR == chr & abs(df$BP - pos) <= 1e6
      if (!any(in_block)) {
        return(tags$p("No causal gene-tissue pairs found in the specified locus.",
                      style = "color: red; font-style: italic;"))
      }
      matched_row <- which(in_block)[1]
    }
    
    chr0 <- df$CHR[matched_row]
    bp0 <- df$BP[matched_row]
    block_locus <- df$Locus[matched_row]
    
    link_name <- trait_info$LinkName[trait_info$ShortName == trait]
    img_prefix <- sprintf("https://tgvis-result.netlify.app/%s", link_name)
    
    fluidRow(
      column(6,
             tags$a(href = sprintf("%s/locuszoom_%s_%s.png", img_prefix, chr0, bp0), target = "_blank",
                    div(class = "magnify-container",
                        tags$img(src = sprintf("%s/locuszoom_%s_%s.png", img_prefix, chr0, bp0),
                                 width = "100%", style = "border: 1px solid #ccc;")))),
      column(6,
             tags$a(href = sprintf("%s/finemap_%s_%s.png", img_prefix, chr0, bp0), target = "_blank",
                    div(class = "magnify-container",
                        tags$img(src = sprintf("%s/finemap_%s_%s.png", img_prefix, chr0, bp0),
                                 width = "100%", style = "border: 1px solid #ccc;")))),
      column(12, br(), h5(sprintf("Gene-tissue pairs in Locus %s (chr%s:%s ±1MB)", block_locus, chr0, bp0)),
             dataTableOutput("df_locus_subset"))
    )
  })
  
  output$show_leading_summary <- renderUI({
    gene <- input$gene_query
    chr <- input$search_chr
    pos <- input$search_pos
    
    if ((is.null(gene) || gene == "") &&
        (is.null(chr) || chr == "") &&
        (is.null(pos) || pos == "")) {
      return(tagList(
        br(), h5("Leading Genes Across All Loci"),
        dataTableOutput("leading_summary")
      ))
    } else {
      return(NULL)
    }
  })
  
  output$leading_summary <- renderDataTable({
    trait <- confirmed_trait()
    req(trait)
    
    path <- file.path(trait, "FineMap.rds")
    df <- readRDS(path)
    
    df %>%
      group_by(Locus) %>%
      slice_max(order_by = CS.Pratt, n = 1, with_ties = FALSE) %>%
      ungroup() %>%
      transmute(
        Locus,
        `Leading Gene` = GeneSymbol,
        `Pratt index of Credible Set` = signif(CS.Pratt, 3)
      ) %>%
      dplyr::arrange(desc(`Pratt index of Credible Set`))
  })
  
  output$df_locus_subset <- renderDataTable({
    trait <- confirmed_trait()
    req(trait)
    
    path <- file.path(trait, "FineMap.rds")
    df <- readRDS(path)
    df$CHR <- as.integer(df$CHR)
    
    gene <- input$gene_query
    chr <- as.integer(input$search_chr)
    pos <- as.numeric(input$search_pos)
    
    if (!is.null(gene) && gene != "") {
      match_idx <- which(tolower(df$GeneSymbol) == tolower(gene) | tolower(df$Gene) == tolower(gene))
      if (length(match_idx) == 0) return(NULL)
      block_locus <- df$Locus[match_idx[1]]
    } else {
      in_block <- df$CHR == chr & abs(df$BP - pos) <= 1e6
      if (!any(in_block)) return(NULL)
      block_locus <- df$Locus[which(in_block)[1]]
    }
    
    df %>%
      dplyr::filter(Locus == block_locus) %>%
      transmute(
        Identifier,
        `Pratt index of Credible set` = signif(CS.Pratt, 3),
        Ensembl = Gene,
        Symbol = GeneSymbol,
        Tissue,
        xQTL,
        Type,
        `Splicing Information`=Splicing
      )
  })
  
  output$download_finemap <- downloadHandler(
    filename = function() {
      trait <- confirmed_trait()
      paste0(trait, "_FineMap.csv")
    },
    content = function(file) {
      trait <- confirmed_trait()
      req(trait)
      path <- file.path(trait, "FineMap.rds")
      df <- readRDS(path)
      full_name <- trait_info$FullName[trait_info$ShortName == trait]
      df_out <- df %>%      mutate(Outcome = full_name) %>%
        transmute(
          Outcome,
          Identifier,
          `Pratt index of Credible set` = signif(CS.Pratt, 3),
          `PIP of Credible set` = signif(CS.PIP, 3),
          Ensembl = Gene,
          Symbol = GeneSymbol,
          Tissue,
          xQTL,
          Type,
          `Splicing Information` = Splicing
        )
      fwrite(df_out, file)
    }
  )
  
  # Gene View reactive data with improved matching
  gene_matched_data <- reactive({
    req(input$gene_search)
    gene_input <- trimws(input$gene_search)
    
    if (gene_input == "") {
      return(data.frame())
    }
    
    # Exact match first
    matched <- w_all %>%
      dplyr::filter(
        Symbol == gene_input | Ensembl == gene_input
      )
    
    # If no exact match, try case-insensitive
    if (nrow(matched) == 0) {
      matched <- w_all %>%
        dplyr::filter(
          tolower(Symbol) == tolower(gene_input) |
            tolower(Ensembl) == tolower(gene_input)
        )
    }
    
    # If still no match, try partial matching
    if (nrow(matched) == 0) {
      matched <- w_all %>%
        dplyr::filter(
          grepl(gene_input, Symbol, ignore.case = TRUE) |
            grepl(gene_input, Ensembl, ignore.case = TRUE)
        )
    }
    
    return(matched)
  })
  
  # Gene table output
  output$gene_table <- renderDataTable({
    matched <- gene_matched_data()
    
    if (nrow(matched) > 0) {
      matched %>%
        dplyr::mutate(
          `Pratt index of credible set` = round(`Pratt index of credible set`, 3),
          `PIP of credible set` = round(`PIP of credible set`, 3)
        )
    } else {
      data.frame(Message = "No results found")
    }
  }, options = list(pageLength = 10, scrollX = TRUE))
  
  # Gene message output
  output$gene_msg <- renderText({
    matched <- gene_matched_data()
    
    if (is.null(input$gene_search) || input$gene_search == "") {
      return("")
    }
    
    if (nrow(matched) == 0) {
      "Target gene is currently not causal in our database."
    } else {
      paste("Found", nrow(matched), "records for", input$gene_search)
    }
  })
  
  # Gene table download
  output$download_gene_table <- downloadHandler(
    filename = function() {
      paste0("GeneView_", gsub("[^a-zA-Z0-9]", "_", input$gene_search), ".csv")
    },
    content = function(file) {
      matched <- gene_matched_data()
      
      if (nrow(matched) > 0) {
        write.csv(matched, file, row.names = FALSE)
      } else {
        write.csv(data.frame(Message = "Gene not found in database."), file, row.names = FALSE)
      }
    }
  )
  
}

shinyApp(ui, server)