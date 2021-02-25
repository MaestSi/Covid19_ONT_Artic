#!/bin/bash

#
# Copyright 2021 Simone Maestri. All rights reserved.
# Simone Maestri <simone.maestri@univr.it>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

PIPELINE_DIR=$(realpath $( dirname "${BASH_SOURCE[0]}" ))
MINICONDA_DIR=$(which conda | sed 's/miniconda3.*$/miniconda3/')
conda config --add channels r
conda config --add channels bioconda

git clone https://github.com/artic-network/artic-ncov2019
git clone https://github.com/artic-network/fieldbioinformatics
cd fieldbioinformatics
sed -i 's/name: artic/name: Covid19_ONT_Artic_env/' environment.yml
conda env create -f environment.yml
source activate Covid19_ONT_Artic_env
python setup.py install
cd ..
conda install seqtk r bioconductor-biostrings r-ggplot2 r-scales NanoFilt
conda deactivate
conda create -n pycoQC_env pip
source activate pycoQC_env
python -m pip install pycoQC
cd "$MINICONDA_DIR/envs/Covid19_ONT_Artic_env/bin"
ln -s "$MINICONDA_DIR/envs/pycoQC_env/bin/pycoQC"
conda deactivate

cd $PIPELINE_DIR

echo -e "\n"
echo "Modify variables PIPELINE_DIR and MINICONDA_DIR in config_Covid19_ONT_Artic.R"
echo -e "PIPELINE_DIR <- \"$PIPELINE_DIR\""
echo -e "MINICONDA_DIR <- \"$MINICONDA_DIR\""
echo -e "\n"
