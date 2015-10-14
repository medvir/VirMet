#!/bin/bash

sudo apt-get update -qq
sudo apt-get install -qq -y curl build-essential ncurses-dev byacc zlib1g-dev python-dev git cmake samtools ncbi-blast+

cd ~/
curl ftp://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/edirect.zip -O
unzip edirect.zip
export PATH=~/edirect:$PATH
