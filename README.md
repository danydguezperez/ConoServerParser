# ConoServerParser: Flexible Shiny Application for Conotoxin FASTA Analysis

## Overview

ConoServerParser is an R/Shiny application designed to complement **ConoServer** by addressing
practical limitations encountered when working with conotoxin sequence data.  Built for
researchers studying the diverse venom peptides produced by marine cone snails, ConoServerParser
streamlines the process of parsing ConoServer’s pipe‑delimited FASTA headers, performing
metadata‑driven filtering and visualisation, and exporting customised sequence sets.

The application was developed and described in the [PagBioMicS blog post](https://www.pagbiomics.com/blog/---conoserverparser--a-flexible-shiny-application-for-filtering--exploring--and-exporting-conotoxin-sequences-from-conoserver)
by Dany Domínguez‑Pérez (June 12 2025).  A video tutorial demonstrating how to use the app
is also available on [YouTube](https://www.youtube.com/watch?v=ZgRrB305xRg).

### Key features

- **FASTA parsing and metadata extraction** – reads ConoServer FASTA files and interprets the
  pipe‑delimited headers to extract fields such as ConoID, species, gene superfamily,
  pharmacological family, cysteine framework and evidence level【256574700114101†L149-L168】.
- **Interactive filtering and exploration** – provides a Shiny interface where you can filter
  sequences by any field (e.g. superfamily, pharmacological family, evidence) and preview the
  resulting selections【256574700114101†L149-L167】.
- **Customisable FASTA export** – generate new FASTA files with headers composed of selected
  metadata fields and custom delimiters; export sequences grouped by superfamily or
  pharmacological family in ZIP archives【256574700114101†L149-L174】.
- **Mature peptide generation** – convert ConoPrec CSV outputs into FASTA files of mature
  peptides, overcoming multi‑sequence limitations in downstream tools【256574700114101†L159-L176】.
- **Dynamic visual summaries** – create interactive bar plots summarising sequence counts across
  superfamilies, pharmacological families or evidence types【256574700114101†L162-L167】.
- **Multi‑format table export** – download filtered metadata tables in TSV, CSV or XLSX formats【256574700114101†L156-L179】.

## Repository contents

This repository contains the source code and example data used by the Shiny app:

| File/folder | Description |
|---|---|
| `ConoServerParser.R` | Main Shiny application script. Running this file launches the app locally. |
| `ConoServerParser_precursors_2025-06-02.fasta` | Example precursor FASTA file downloaded from ConoServer (dated 2025‑06‑02). |
| `data/` | Contains auxiliary data files used by the application. |
| `data/conoserver_protein.fa` | A larger protein FASTA file retrieved from ConoServer. |
| `data/Y_superfamily.csv` | Reference table mapping ConoServer ConoIDs to gene superfamily designations. |

Deployment‑specific files (`rsconnect/`) and temporary files (`.RData`, `.Rhistory`, `.DS_Store`) are
omitted from this repository to keep it clean and portable.

## Running the application

To run ConoServerParser locally you will need a recent version of **R** (≥ 4.x) and the
packages listed below.  Launching the app opens a web browser window with the Shiny user
interface.

1. Install the required R packages:

   ```r
   install.packages(c(
     "shiny", "tidyverse", "data.table", "Biostrings", "plotly", 
     "readxl", "openxlsx"
   ))
   ```

2. Clone this repository and set your working directory to its root:

   ```bash
   git clone https://github.com/<your‑username>/ConoServerParser.git
   cd ConoServerParser
   ```

3. Launch the app from within R:

   ```r
   shiny::runApp("ConoServerParser.R")
   ```

4. Use the sidebar controls in the web interface to upload ConoServer FASTA files, adjust
   filters, explore sequence distributions, and export customised FASTA or tabular files.

## Online application

If you prefer not to run the app locally, a hosted version is available at
[ConoServerParser on shinyapps.io](https://danysaurio.shinyapps.io/ConoServerParser)
(availability may vary depending on usage limits).

## Citation

If you use ConoServerParser in your research, please cite the PagBioMicS blog article and
acknowledge Dany Domínguez‑Pérez:

> Domínguez‑Pérez, D. (2025). *ConoServerParser: A Flexible Shiny Application for Filtering,
> Exploring, and Exporting Conotoxin Sequences from ConoServer*. PagBioMicS Blog, June 12 2025【256574700114101†L20-L33】.

## License

This software is provided under an open‑source licence.  Feel free to modify and adapt the
app for academic or non‑commercial use; please attribute the original author when
distributing derivative works.
