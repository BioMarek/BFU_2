#! /bin/bash
#

CONDA_BASE=$(conda info --base 2>/dev/null)
if [ -n "$CONDA_BASE" ]; then
    source "$CONDA_BASE/etc/profile.d/conda.sh"
    conda activate trash2
fi

cd "$HOME/marek/programming/BFU/TRASH_2/" || exit 1

INPUT_DIR="$HOME/marek/programming/BFU/TRASH_2/input"
OUTPUT_DIR="$HOME/marek/programming/BFU/TRASH_2/output"

for file in "$INPUT_DIR"/*.fasta; do
    filename=$(basename "$file")
    output_name="${filename%.fasta}"
    mkdir -p "$OUTPUT_DIR/$output_name"

    Rscript src/TRASH.R -f "$file" -o "$OUTPUT_DIR/$output_name"
done
