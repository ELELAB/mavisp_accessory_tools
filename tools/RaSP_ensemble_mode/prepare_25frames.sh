#!/bin/bash

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <input_pdb>"
    exit 1
fi

INPUT_PDB="$1"

if [ ! -f "$INPUT_PDB" ]; then
    echo "Error: input file $INPUT_PDB not found."
    exit 1
fi

MODEL_COUNT=$(grep -c '^MODEL' "$INPUT_PDB")

if [ "$MODEL_COUNT" -eq 0 ]; then
    echo "Error: no MODEL records found in $INPUT_PDB."
    echo "This usually means the PDB has a single structure or MODEL designations are missing."
    exit 1
fi

if [ "$MODEL_COUNT" -eq 1 ]; then
    echo "Error: only one MODEL record found in $INPUT_PDB."
    echo "The file appears to contain a single model, not multiple frames."
    exit 1
fi

echo "Found $MODEL_COUNT MODEL records. Splitting into separate frame directories..."

awk '
    /^MODEL/ {
        model++
        dir = sprintf("frame%d", model)
        system("mkdir -p " dir)
        outfile = sprintf("%s/model_%d.pdb", dir, model)
        next
    }

    /^ENDMDL/ {
        if (outfile != "") {
            close(outfile)
            outfile = ""
        }
        next
    }

    # Only write ATOM/HETATM/CRYST1/etc lines inside a model block
    /^ATOM/ || /^HETATM/ || /^CRYST1/ || /^TER/ {
        if (outfile != "")
            print >> outfile
    }
' "$INPUT_PDB"



