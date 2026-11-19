python HiCAT.py -i ./testdata/cen21.fa -t ./testdata/AlphaSat.fa


#!/bin/bash
TARGET_DIR="/home/marek/marek/programming/BFU/centromeres/persica"
MONOMER="/home/marek/marek/programming/BFU/centromeres/monomers/CANR4.fa"

for input_file in "$TARGET_DIR"/*.fa; do
    name="${input_file%.fa}"
    output_dir="${TARGET_DIR}/${name}"
    mkdir -p "$output_dir"

    echo "Processing $input_file -> $name"

    python HiCAT.py -th 12 -i "$input_file" -t "$MONOMER" -o "$name"
done
