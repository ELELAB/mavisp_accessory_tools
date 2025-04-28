#!/bin/bash

source /usr/local/envs/py310/bin/activate

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

INPUT_FILE="$1"
OUTPUT_FILE="${2}.csv"

if [[ -z "$INPUT_FILE" || -z "$OUTPUT_FILE" ]]; then
    echo "Use: $0 <input_file> <output_file>"
    exit 1
fi

rm -f "$OUTPUT_FILE"
echo "protein,cofactors" > "$OUTPUT_FILE"

while IFS= read -r upid; do
    echo "Processing $upid"

    cofactors=$(python "$SCRIPT_DIR/find_cofactor.py" -u "$upid" 2>/dev/null)

    [[ -z "$cofactors" ]] && cofactors="Error"
    echo "$cofactors" >> "$OUTPUT_FILE"


done < "$INPUT_FILE"


