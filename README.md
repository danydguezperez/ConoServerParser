# 🧪 ConoServerParser: A Flexible Shiny Application for Filtering, Exploring, and Exporting Conotoxin Sequences from ConoServer

## Overview

ConoServerParser is an R/Shiny application designed to complement **ConoServer** by addressing
practical limitations encountered when working with conotoxin sequence data.  Built for
researchers studying the diverse venom peptides produced by marine cone snails, ConoServerParser
streamlines the process of parsing ConoServer’s pipe‑delimited FASTA headers, performing
metadata‑driven filtering and visualisation, and exporting customised sequence sets.

Conus snails (genus *Conus*) produce venom peptides known as conopeptides, which are highly
specific for ion channels and receptors【953987192211277†L20-L25】.  This specificity makes
conotoxins invaluable for physiological studies and potential drug development【953987192211277†L20-L25】.
With more than 700 *Conus* species and thousands of distinct peptides, conotoxin research spans
evolutionary biology and neurophysiology【953987192211277†L26-L37】.  Conopeptides are classified
into disulfide‑rich and disulfide‑poor peptides and categorised by pharmacological family,
gene superfamily and cysteine framework【953987192211277†L39-L49】.  These classifications and
the sheer diversity of conotoxins underscore the need for tools like ConoServerParser to
facilitate downstream analysis of ConoServer data.

The application was developed and described in the [PagBioMicS blog post](https://www.pagbiomics.com/blog/---conoserverparser--a-flexible-shiny-application-for-filtering--exploring--and-exporting-conotoxin-sequences-from-conoserver)
by Dany Domínguez‑Pérez (June 12 2025).  A video tutorial demonstrating how to use the app
is also available on [YouTube](https://www.youtube.com/watch?v=ZgRrB305xRg).

### Key features

ConoServerParser brings together a suite of tools designed to make ConoServer data easier to
work with.  The application implements the capabilities described in the PagBioMicS article
and groups overlapping items where appropriate:

 - **FASTA parsing and metadata extraction** – Reliably interprets ConoServer’s
   pipe‑delimited headers to extract fields such as ConoID, species, gene superfamily,
   pharmacological family, cysteine framework and evidence level into a clean tabular
   format【953987192211277†L149-L152】.
 - **Customisable FASTA export and header preview** – Allows you to construct new FASTA
   headers using selected metadata fields and custom delimiters, preview how the final
   headers will look and avoid oversized or unnecessary fields.  You can also group
   downloads by superfamily or pharmacological category and export grouped archives as
   ZIP files【953987192211277†L153-L158】【953987192211277†L168-L174】.
 - **Multi‑format table export** – Enables downloading of filtered sequence metadata tables
   in CSV, TSV or XLSX formats for downstream analysis【953987192211277†L156-L158】【953987192211277†L179-L179】.
 - **Mature peptide FASTA generation** – Automatically extracts and formats mature
   sequences from ConoPrec CSV outputs, overcoming the multi‑sequence limitations of
   downstream tools【953987192211277†L159-L161】【953987192211277†L175-L176】.
 - **Interactive filtering and exploration** – Provides an intuitive interface with
   dropdowns and search tools to select, filter and exclude sequences by any field
   (superfamily, pharmacological family, evidence level, etc.)【953987192211277†L162-L164】【953987192211277†L177-L177】.
 - **Dynamic visual summaries and data visualisation** – Generates interactive bar plots
   comparing sequence distributions across superfamilies, pharmacological families or
   evidence types and other visualisations to summarise your data【953987192211277†L165-L166】【953987192211277†L178-L178】.

## Repository contents

This repository contains the source code and example data used by the Shiny app:

| File/folder | Description |
|---|---|
| `ConoServerParser.R` | Main Shiny application script. Running this file launches the app locally. |
| `ConoServerParser_precursors_2025-06-02.fasta` | Example precursor FASTA file exported from ConoServer (dated 2025‑06‑02) containing sequence data used to demonstrate the app. |
| `data/` | Contains auxiliary data files used by the application. |
| `data/conoserver_protein.fa` | FASTA file containing a larger set of protein sequences from ConoServer. |
| `data/Y_superfamily.csv` | Reference table mapping ConoServer ConoIDs to Y superfamily designations, used for example analyses. |

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
> Exploring, and Exporting Conotoxin Sequences from ConoServer*. PagBioMicS Blog, June 12 2025【953987192211277†L11-L15】.

## License

This software is provided under an open‑source licence.  Feel free to modify and adapt the
app for academic or non‑commercial use; please attribute the original author when
distributing derivative works.
