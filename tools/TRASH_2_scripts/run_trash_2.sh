#! /bin/bash
#

CONDA_BASE=$(conda info --base 2>/dev/null)
if [ -n "$CONDA_BASE" ]; then
    source "$CONDA_BASE/etc/profile.d/conda.sh"
    conda activate trash2
fi

cd "$HOME/sochorova/TRASH_2/" || exit 1

INPUT_DIR="$HOME/sochorova/TRASH_2/input"
OUTPUT_DIR="$HOME/sochorova/TRASH_2/output"

for file in "$INPUT_DIR"/*.fa; do
    filename=$(basename "$file")
    output_name="${filename%.fasta}"
    mkdir -p "$OUTPUT_DIR/$output_name"

    Rscript src/TRASH.R -p 24 -f "$file" -o "$OUTPUT_DIR/$output_name"
done
