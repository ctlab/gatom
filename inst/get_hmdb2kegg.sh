#!/usr/bin/env bash

# Downloading Metabocards flat file
# 1GB in size!
wget "http://www.hmdb.ca/downloads/2_5/metabocards.zip" -O metabocards-2.5.zip

unzip -p metabocards-2.5.zip | \
    grep -A 1 -e "^# kegg_compound_id:" -e "BEGIN_METABOCARD"  | \
    grep -v -e -- -e '^$' -e "kegg_compound_id" |\
    paste - - | \
    sed "s/^#BEGIN_METABOCARD //" | \
    sed "s/Not available/NA/i" | \
    sed "s/^\(.*\)\t\(\w\+\); /\1\t\2\n\1\t/" | \
    sed "s/\(HMDB00656\t\w*\t\)NA/\1C05526/" | \
    grep -v 'NA' |\
    sed "1i HMDB\tKEGG" > hmdb2kegg.tsv
