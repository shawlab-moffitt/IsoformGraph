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
library(tidyr)

edb <- EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86

appris_zip <- "appris.zip"

if (!file.exists(appris_zip)) {
  download.file(
    "https://raw.githubusercontent.com/shawlab-moffitt/appris_primary_transcripts/main/Gencode48_Ensembl114_appris_data.principal_score.txt.zip",
    destfile = appris_zip,
    mode = "wb"
  )
}

appris <- readr::read_tsv(
  unz(
    appris_zip,
    "Gencode48_Ensembl114_appris_data.principal_score.txt"
  ),
  col_names = TRUE,
  show_col_types = FALSE
)

ui <- shiny::fluidPage(
  shiny::titlePanel("IsoformGraph: Exon Splicing Graph"),
  
  shiny::sidebarLayout(
    shiny::sidebarPanel(
      shiny::textInput("gene_symbol", "Gene symbol", value = "PTPRC"),
      
      shiny::numericInput(
        "plot_width",
        "Plot width (pixels)",
        value = 1200,
        min = 400,
        max = 5000,
        step = 100
      ),
      
      shiny::numericInput(
        "plot_height",
        "Plot height (pixels)",
        value = 900,
        min = 400,
        max = 5000,
        step = 100
      ),
      
      shiny::actionButton("run", "Generate plot")
    ),
    
    shiny::mainPanel(
      shiny::uiOutput("plot_ui"),
      shiny::h4("Exon nodes and skipped-exon annotation"),
      shiny::tableOutput("tx_table")
    )
  )
)

server <- function(input, output, session) {
  
  output$plot_ui <- shiny::renderUI({
    shiny::plotOutput(
      "combined_plot",
      width = paste0(input$plot_width, "px"),
      height = paste0(input$plot_height, "px")
    )
  })
  
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
        `Transcript support level` >= 1
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
        exon_3p = dplyr::if_else(strand == "+", end, start),
        exon_node = paste0(start, "-", end)
      ) %>%
      dplyr::distinct(
        exon_node,
        start,
        end,
        strand,
        gene_exon_id,
        exon_5p,
        exon_3p
      ) %>%
      dplyr::arrange(exon_5p) %>%
      dplyr::mutate(
        exon_order = dplyr::row_number(),
        exon_label = paste0("E", gene_exon_id)
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
            exon_start = start,
            exon_end = end,
            x_start = exon_order,
            donor_label = exon_label
          ),
        by = "exon_node"
      ) %>%
      dplyr::left_join(
        exon_nodes_v2 %>%
          dplyr::select(
            next_exon_node = exon_node,
            next_exon_start = start,
            next_exon_end = end,
            x_end = exon_order,
            acceptor_label = exon_label
          ),
        by = "next_exon_node"
      ) %>%
      dplyr::mutate(
        edge_id = dplyr::row_number(),
        junction_start = pmin(exon_end, next_exon_start),
        junction_end = pmax(exon_end, next_exon_start),
        junction_label = paste0(donor_label, " -> ", acceptor_label),
        junction_coord = paste0(exon_node, " -> ", next_exon_node)
      )
    
    skipped_hits <- splice_edges_v2 %>%
      dplyr::left_join(
        exon_nodes_v2 %>%
          dplyr::select(
            skipped_exon_node = exon_node,
            skipped_exon_label = exon_label,
            skipped_start = start,
            skipped_end = end
          ),
        by = character()
      ) %>%
      dplyr::filter(
        skipped_exon_node != exon_node,
        skipped_exon_node != next_exon_node,
        skipped_start > junction_start,
        skipped_end > junction_start,
        skipped_start < junction_end,
        skipped_end < junction_end
      )
    
    splice_edges_v2 <- splice_edges_v2 %>%
      dplyr::left_join(
        skipped_hits %>%
          dplyr::select(
            edge_id,
            skipped_exon_node,
            skipped_exon_label,
            skipped_start,
            skipped_end
          ),
        by = "edge_id"
      ) %>%
      dplyr::mutate(
        skips_exon = !is.na(skipped_exon_node)
      )
    
    skip_annotation <- skipped_hits %>%
      dplyr::group_by(skipped_exon_node) %>%
      dplyr::summarise(
        skip_status = "skipped_by_junction",
        skipping_junction = paste(unique(junction_label), collapse = "; "),
        skipping_junction_coord = paste(unique(junction_coord), collapse = "; "),
        n_skipping_junctions = dplyr::n_distinct(edge_id),
        .groups = "drop"
      )
    
    exon_nodes_v2 <- exon_nodes_v2 %>%
      dplyr::left_join(
        skip_annotation,
        by = c("exon_node" = "skipped_exon_node")
      ) %>%
      dplyr::mutate(
        skip_status = tidyr::replace_na(skip_status, "not_skipped"),
        skipping_junction = tidyr::replace_na(skipping_junction, ""),
        skipping_junction_coord = tidyr::replace_na(skipping_junction_coord, ""),
        n_skipping_junctions = tidyr::replace_na(n_skipping_junctions, 0L)
      )
    
    list(
      gene_symbol = gene_symbol,
      appris_pc_gene = appris_pc_gene,
      ex_plot = ex_plot,
      exon_nodes_v2 = exon_nodes_v2,
      splice_edges_v2 = splice_edges_v2,
      skipped_hits = skipped_hits
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
          ymax = 1.0,
          fill = skip_status
        ),
        color = "black"
      ) +
      ggplot2::scale_fill_manual(
        values = c(
          "not_skipped" = "grey85",
          "skipped_by_junction" = "tomato"
        ),
        name = "Exon status"
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
    
    dat$exon_nodes_v2 %>%
      dplyr::select(
        exon_label,
        exon_node,
        start,
        end,
        strand,
        skip_status,
        n_skipping_junctions,
        skipping_junction,
        skipping_junction_coord
      )
  })
}

shiny::shinyApp(ui = ui, server = server)
