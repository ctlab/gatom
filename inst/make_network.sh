#!/usr/bin/env bash

# This script creates carbon atom network from KEGG database

# 1. Downloading public rpair file
wget ftp://ftp.genome.jp/pub/db/rclass/rpair

# 2. Processing atom alignments (takes ~15 min)
mkdir network
python3 process_rpair.py rpair

# 3. Making network into an R object
# requires KEGGREST and data.table R packages
Rscript make_network.R 
