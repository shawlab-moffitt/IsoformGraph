# IsoformGraph

**IsoformGraph** is a Shiny application for visualizing protein-coding transcript exon structures and gene-level splice graphs.

🔗 **Live app:** https://shawlab-moffitt.shinyapps.io/isoformgraph/

## Overview

IsoformGraph allows users to enter a gene symbol and generate an exon-level visualization of transcript structure using Ensembl annotation and APPRIS transcript information. The app highlights exon organization across protein-coding transcripts and summarizes splice connections between ordered exons.

The current app uses `CD44` as the default example gene and includes a transcript support level filter. The Shiny interface contains a gene-symbol input, transcript-support cutoff, a generate button, a combined exon/splice-graph plot, and an exon-node table. :contentReference[oaicite:0]{index=0}

## Features

- Visualize exon structures across protein-coding transcripts
- Generate gene-level splice graphs
- Label exons by genomic order
- Filter transcripts using APPRIS transcript annotation
- Display exon node information used to construct the splice graph
- Run locally or access through the hosted Shiny app

## Live Web Application

The IsoformGraph app is publicly available here:

**https://shawlab-moffitt.shinyapps.io/isoformgraph/**

## Installation

Clone the repository:

```r
git clone https://github.com/shawlab-moffitt/IsoformGraph.git
