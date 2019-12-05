#!/bin/bash

set -e

CONDA_DIR="${1:-"${HOME}/work/miniconda"}"
CONDA_ENV="${2:-'ngs'}"
SETUP_DIR="${3:-'.'}"
SERVER="${4:-'BIH'}"

echo "$SETUP_DIR"

# get absolute path
CONDA_DIR="$(readlink -f "$CONDA_DIR")"
SETUP_DIR="$(readlink -f "$SETUP_DIR")"
echo "$SETUP_DIR"



# make NGS folders
cd "$SETUP_DIR"
mkdir NGS
mkdir NGS/links
mkdir NGS/tools
mkdir NGS/refgenome
mkdir NGS/known_sites


# get ngs repo
cd "$SETUP_DIR"
cd NGS
git clone "https://github.com/tjblaette/ngs"
if [[ "$SERVER" == 'BIH' ]]
then
	cd ngs
	git checkout bih
fi


# get BWA MEM
cd "$SETUP_DIR"
cd NGS/tools
git clone https://github.com/lh3/bwa bwa_0.7.10
cd bwa_0.7.10
git checkout 0.7.10
make


# get VarScan2
cd "$SETUP_DIR"
cd NGS/tools
mkdir varscan_v2.3.9
cd varscan_v2.3.9
wget https://sourceforge.net/projects/varscan/files/VarScan.v2.3.9.jar
cp "${SETUP_DIR}/NGS/ngs/run_jar.sh" .
mv run_jar.sh run_varscan_v2.3.9.sh


# get bpipe
cd "$SETUP_DIR"
cd NGS/tools
git clone https://github.com/ssadedin/bpipe bpipe_0.9.9.6
cd bpipe_0.9.9.6
git checkout 0.9.9.6
mkdir JAVA_TMP
./gradlew clean jar
rm -r JAVA_TMP


# get fastQC
cd "$SETUP_DIR"
cd NGS/tools
mkdir fastqc_0.11.8
cd fastqc_0.11.8
wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.8.zip
unzip fastqc*.zip
rm fastqc*.zip
mv FastQC/* .
rmdir FastQC/
chmod ug+x fastqc



# link tools
cd "$SETUP_DIR"
cd NGS/links
ln -s ../tools/bwa_0.7.10/bwa bwa
ln -s ../tools/varscan_v2.3.9/run_varscan_v2.3.9.sh varscan
ln -s ../tools/bpipe_0.9.9.6/bin/bpipe bpipe
ln -s ../tools/NGSQCToolkit_v2.3.3/QC/IlluQC_PRLL.pl ngsqc
ln -s ../tools/fastqc_0.11.8/fastqc fastqc
ln -s ../tools/GenomeAnalysisTK-3.4-46/run_jar.sh gatk
ln -s $(which STAR) star # conda executable is 'STAR', bpipe expects 'star'
# add links folder to PATH later, when adding miniconda


# download and install newest miniconda3
cd "$SETUP_DIR"
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b -f -p "${CONDA_DIR}"
rm Miniconda3-latest-Linux-x86_64.sh


# add miniconda3 and NGS/links to $PATH
if [[ "$SERVER" == 'BIH' ]]
then

echo -n '
case "${HOSTNAME}" in
    med-login*)
	;;
    *)
	export PATH=' >> ~/.bashrc
echo -n "${SETUP_DIR}/NGS/links:${CONDA_DIR}/bin" >> ~/.bashrc
echo ':$PATH
	;;
esac' >> ~/.bashrc

else

echo -n 'export PATH=' >> ~/.bashrc
echo -n "${SETUP_DIR}/NGS/links:${CONDA_DIR}/bin" >> ~/.bashrc
echo ':$PATH' >> ~/.bashrc

fi

echo $PATH
source ~/.bashrc
echo $PATH


# create a conda environment to install tools to
echo "Setting up conda environment, this might take a very long time..."
conda config --add channels bioconda
conda config --add channels conda-forge
conda create --name "${CONDA_ENV}" python=3.5 cutadapt="1.18" samtools="0.1.19" picard="2.18.26" gatk="3.8" bedtools="2.29.0" fgbio="0.8.0" star="2.4.2a" htseq="0.11.2"
source activate "${CONDA_ENV}"

echo -e "Installing conda R packages, this might take a very long time..."
conda install r-stringr r-ggthemes r-reshape2 r-ggplot2 r-pheatmap bioconductor-deseq2



# prompt GATK obtainment
echo -e "\nPlease obtain GATK 3.4.46 online or from files-for-bih-server"
echo -e "Save it to NGS/tools and check the link in NGS/links\n"


# prompt ANNOVAR obtainment
echo -e "\nPlease obtain ANNOVAR from http://annovar.openbioinformatics.org"
echo -e "Save it to NGS/tools and link it in NGS/links\n"


# get references and annotation files?
echo -e "\nPlease obtain reference genome sequence and annotation files"
echo -e "Save them to NGS/refgenome and NGS/known_sites respectively"
echo -e "Create the necessary index files for each genome using 'make_genomeFiles.sh'\n"


# report success!
echo "Successfully completed the automated part of your environment's setup!"
echo "Please follow the instructions above to complete your setup."


