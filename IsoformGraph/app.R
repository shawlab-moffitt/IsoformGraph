# app.R

library(shiny)
library(EnsDb.Hsapiens.v86)
library(ensembldb)
library(AnnotationFilter)
library(dplyr)
library(tibble)
library(readr)
library(ggplot2)
library(ggtranscript)
library(patchwork)

edb <- EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86

#appris_url <- "https://apprisws.bioinfo.cnio.es/pub/current_release/datafiles/homo_sapiens/e110v48/appris_data.appris.txt"
#https://github.com/shawlab-moffitt/appris_primary_transcripts/blob/main/Gencode48_Ensembl114_appris_data.principal_score.txt.zip

url <- "https://raw.githubusercontent.com/shawlab-moffitt/appris_primary_transcripts/main/Gencode48_Ensembl114_appris_data.principal_score.txt.zip"

download.file(
  "https://raw.githubusercontent.com/shawlab-moffitt/appris_primary_transcripts/main/Gencode48_Ensembl114_appris_data.principal_score.txt.zip",
  destfile = "appris.zip",
  mode = "wb"
)

appris <- readr::read_tsv(
  unz("appris.zip",
      "Gencode48_Ensembl114_appris_data.principal_score.txt"),
  col_names = TRUE,
  show_col_types = FALSE
)

ui <- shiny::fluidPage(
  shiny::titlePanel("Exon Splicing Graph"),
  
  shiny::sidebarLayout(
    shiny::sidebarPanel(
      shiny::textInput("gene_symbol", "Gene symbol", value = "PTPRC"),
      shiny::numericInput(
        "tsl_cutoff",
        "Minimum transcript support level",
        value = 1,
        min = 0
      ),
      shiny::actionButton("run", "Generate plot")
    ),
    
    shiny::mainPanel(
      shiny::plotOutput("combined_plot", height = "900px"),
      shiny::h4("APPRIS transcripts used"),
      shiny::tableOutput("tx_table")
    )
  )
)

server <- function(input, output, session) {
  
  gene_data <- shiny::eventReactive(input$run, {
    
    gene_symbol <- toupper(input$gene_symbol)
    
    appris_pc_gene <- appris %>%
      dplyr::mutate(
        `Transcript support level` =
          suppressWarnings(as.numeric(`Transcript support level`))
      ) %>%
      dplyr::filter(
        `Transcript type` == "protein_coding",
        `Gene name (HGNC)` == gene_symbol,
        `Transcript support level` >= input$tsl_cutoff
      )
    
    shiny::validate(
      shiny::need(
        nrow(appris_pc_gene) > 0,
        "No APPRIS protein-coding transcripts found for this gene."
      )
    )
    
    ex <- ensembldb::exons(
      edb,
      filter = AnnotationFilter::AnnotationFilterList(
        AnnotationFilter::GeneNameFilter(gene_symbol),
        AnnotationFilter::TxBiotypeFilter("protein_coding")
      ),
      columns = c(
        "gene_name", "gene_id",
        "tx_id", "tx_name", "tx_biotype",
        "exon_id", "exon_idx",
        "seq_name", "seq_strand",
        "exon_seq_start", "exon_seq_end"
      ),
      return.type = "data.frame"
    )
    
    shiny::validate(
      shiny::need(
        nrow(ex) > 0,
        "No protein-coding exons found in EnsDb for this gene."
      )
    )
    
    common_tx <- base::intersect(
      base::unique(ex$tx_id),
      base::unique(appris_pc_gene$`Transcript ID`)
    )
    
    shiny::validate(
      shiny::need(
        length(common_tx) > 0,
        "No overlap between EnsDb transcript IDs and APPRIS transcript IDs."
      )
    )
    
    ex <- ex %>%
      dplyr::filter(tx_id %in% common_tx)
    
    ex2 <- ex %>%
      tibble::as_tibble() %>%
      dplyr::rename(
        transcript_id = tx_id,
        transcript_name = tx_name,
        transcript_biotype = tx_biotype,
        start = exon_seq_start,
        end = exon_seq_end,
        strand = seq_strand
      ) %>%
      dplyr::filter(transcript_biotype == "protein_coding")
    
    exon_lookup <- ex2 %>%
      dplyr::mutate(
        exon_5p = dplyr::if_else(strand == "+", start, end)
      ) %>%
      dplyr::distinct(start, end, strand, exon_5p) %>%
      dplyr::arrange(exon_5p) %>%
      dplyr::mutate(
        gene_exon_id = dplyr::row_number()
      ) %>%
      dplyr::select(start, end, strand, gene_exon_id)
    
    ex2 <- ex2 %>%
      dplyr::left_join(
        exon_lookup,
        by = c("start", "end", "strand")
      )
    
    ex_plot <- ex2 %>%
      dplyr::filter(!is.na(transcript_name)) %>%
      dplyr::group_by(transcript_name) %>%
      dplyr::mutate(n_exons = dplyr::n()) %>%
      dplyr::ungroup()
    
    exon_nodes_v2 <- ex2 %>%
      dplyr::mutate(
        exon_5p = dplyr::if_else(strand == "+", start, end),
        exon_3p = dplyr::if_else(strand == "+", end, start)
      ) %>%
      dplyr::distinct(
        start, end, strand,
        gene_exon_id,
        exon_5p, exon_3p
      ) %>%
      dplyr::arrange(exon_5p) %>%
      dplyr::mutate(
        exon_order = gene_exon_id,
        exon_label = paste0("E", gene_exon_id),
        exon_node = paste0(start, "-", end)
      )
    
    splice_edges_v2 <- ex2 %>%
      dplyr::mutate(
        exon_5p = dplyr::if_else(strand == "+", start, end),
        exon_3p = dplyr::if_else(strand == "+", end, start),
        exon_node = paste0(start, "-", end)
      ) %>%
      dplyr::arrange(transcript_id, exon_5p) %>%
      dplyr::group_by(transcript_id, transcript_name) %>%
      dplyr::mutate(
        next_exon_node = dplyr::lead(exon_node)
      ) %>%
      dplyr::filter(!is.na(next_exon_node)) %>%
      dplyr::ungroup() %>%
      dplyr::count(
        exon_node,
        next_exon_node,
        name = "n_transcripts"
      ) %>%
      dplyr::left_join(
        exon_nodes_v2 %>%
          dplyr::select(
            exon_node,
            x_start = exon_order
          ),
        by = "exon_node"
      ) %>%
      dplyr::left_join(
        exon_nodes_v2 %>%
          dplyr::select(
            next_exon_node = exon_node,
            x_end = exon_order
          ),
        by = "next_exon_node"
      )
    
    list(
      gene_symbol = gene_symbol,
      appris_pc_gene = appris_pc_gene,
      ex_plot = ex_plot,
      exon_nodes_v2 = exon_nodes_v2,
      splice_edges_v2 = splice_edges_v2
    )
  })
  
  output$combined_plot <- shiny::renderPlot({
    
    dat <- gene_data()
    
    g_structure <- ggplot2::ggplot(
      dat$ex_plot,
      ggplot2::aes(
        xstart = start,
        xend = end,
        y = transcript_name
      )
    ) +
      ggtranscript::geom_range(fill = "steelblue") +
      ggtranscript::geom_intron(
        data = ggtranscript::to_intron(
          dat$ex_plot,
          group_var = "transcript_name"
        ),
        ggplot2::aes(
          xstart = start,
          xend = end,
          y = transcript_name
        ),
        arrow.min.intron.length = 500
      ) +
      ggplot2::geom_text(
        ggplot2::aes(
          x = (start + end) / 2,
          label = gene_exon_id
        ),
        color = "red",
        size = 4,
        fontface = "bold",
        vjust = -1.3
      ) +
      ggplot2::theme_bw() +
      ggplot2::labs(
        title = paste0(dat$gene_symbol, " transcript exon structures"),
        x = "Genomic position",
        y = "Transcript"
      )
    
    g_arc <- ggplot2::ggplot() +
      ggplot2::geom_rect(
        data = dat$exon_nodes_v2,
        ggplot2::aes(
          xmin = exon_order - 0.35,
          xmax = exon_order + 0.35,
          ymin = 0.8,
          ymax = 1.0
        ),
        fill = "grey85",
        color = "black"
      ) +
      ggplot2::geom_text(
        data = dat$exon_nodes_v2,
        ggplot2::aes(
          x = exon_order,
          y = 0.75,
          label = exon_label
        ),
        size = 3
      ) +
      ggplot2::geom_curve(
        data = dat$splice_edges_v2,
        ggplot2::aes(
          x = x_start,
          y = 1.05,
          xend = x_end,
          yend = 1.05,
          linewidth = n_transcripts
        ),
        curvature = -0.35,
        arrow = grid::arrow(length = grid::unit(0, "mm")),
        alpha = 0.7
      ) +
      ggplot2::scale_linewidth(range = c(0.2, 2.5)) +
      ggplot2::coord_cartesian(
        ylim = c(0.7, 2.0),
        clip = "off"
      ) +
      ggplot2::theme_void() +
      ggplot2::labs(
        title = paste0(dat$gene_symbol, " gene-level splice graph")
      )
    
    g_structure / g_arc +
      patchwork::plot_layout(heights = c(4, 1))
  })
  
  output$tx_table <- shiny::renderTable({
    
    dat <- gene_data()
    dat$exon_nodes_v2
    #dat$appris_pc_gene %>%
    #  dplyr::select(
    #    `Gene name (HGNC)`,
    #    `Transcript ID`,
    #    `Transcript type`,
    #    `Transcript support level`
    #    #`APPRIS annotation`
    #  )
  })
}

shiny::shinyApp(ui = ui, server = server)