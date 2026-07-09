# installation
conda create -n trash2 -c conda-forge r-base zlib gcc libcurl

cd /TRASH_2/src
chmod 700 TRASH.R
chmod 700 HORT.R

mkdir -p /TRASH_2/input
mkdir -p /TRASH_2/output

chmod 700 run_trash_2.sh

# manual run
conda activate trash2
Rscript TRASH.R