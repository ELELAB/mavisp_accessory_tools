#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <chain_id> <cores> <num_frames>"
    exit 1
fi

CHAIN_ID=$1
CORES=$2
NUM_FRAMES=$3

# Load RaSP conda environment
source /usr/local/miniconda3/etc/profile.d/conda.sh
conda activate /usr/local/envs/RaSP_workflow

# Generate cluster directories dynamically
directories=()
for ((i=1; i<=NUM_FRAMES; i++)); do
    directories+=("frame$i")
done

# Step 1: Add chain if missing
add_chain() {
    awk -v chain="$CHAIN_ID" '{
        if ($0 ~ /^ATOM/ || $0 ~ /^HETATM/) {
            if (substr($0, 22, 1) == " ") {
                print substr($0, 1, 21) chain substr($0, 23);
            } else {
                print $0;
            }
        } else {
            print $0;
        }
    }' "$1" > tmp_chain.pdb
}

# Step 2: Remove solvent and hydrogen
remove_solvent_and_hydrogen() {
    pdb_element "$1" > tmp_element.pdb
    pdb_delelem -H tmp_element.pdb > tmp_noH.pdb
    rm tmp_element.pdb
}

# Step 3: Fix amino acid names (specific to histidine and others if needed)
fix_amino_acid_names() {
    awk '{
        if ($0 ~ /^(ATOM|HETATM)/) {
            res = substr($0,18,4)
            if (res ~ /HIE|HID|HIP|HISD|HISE|HISH/) res="HIS "
            else if (res=="LYN ") res="LYS "
            else if (res=="ASH ") res="ASP "
            else if (res=="GLH ") res="GLU "
            else if (res=="CYX ") res="CYS "
            $0 = substr($0,1,17) res substr($0,22)
        }
        print
    }' "$1" > tmp_fixed_names.pdb
}

# Step 4: Pad the file twice and ensure lines are 80 characters
pad_file() {
    sed -E 's/.{1,79}$/&                                                                            /' "$1" | cut -c1-80 > tmp_padded1.pdb
    sed -E 's/.{1,79}$/&                                                                            /' tmp_padded1.pdb | cut -c1-80 > "$2"
    rm tmp_padded1.pdb
}

# Step 5: Remove ending lines such as END/TER/ENDMOL
remove_end_lines() {
    awk '/^(ATOM|HETATM)/' "$1" > tmp_no_end_lines.pdb
}

for dir in "${directories[@]}"; do
    # Check if the directory exists
    echo "=== Processing cluster: $dir ==="
    if [ -d "$dir" ]; then
        cd "$dir"

        # Find the PDB file in the directory
        pdb_file=$(ls *.pdb 2>/dev/null | head -n 1)
        if [ -z "$pdb_file" ]; then
            echo "No PDB file found in $dir. Skipping."
            cd - > /dev/null
            continue
        fi

        OUTPUT_PDB="${pdb_file%.pdb}_processed.pdb"

        # Execute processing steps
        fix_amino_acid_names "$pdb_file"
        add_chain tmp_fixed_names.pdb
        remove_solvent_and_hydrogen tmp_chain.pdb
        remove_end_lines tmp_noH.pdb
        pad_file tmp_no_end_lines.pdb "$OUTPUT_PDB"

        # Clean up intermediate files
        rm tmp_chain.pdb tmp_noH.pdb tmp_fixed_names.pdb tmp_no_end_lines.pdb

        echo "Processing complete. Output written to $OUTPUT_PDB"

        # Run RaSP_workflow command
        RaSP_workflow -i $OUTPUT_PDB -r cpu -p /usr/local/envs/RaSP_workflow/RaSP_workflow/src/ -o . -n $CORES -c $CHAIN_ID

        # Run RaSP_postprocess command
        RaSP_postprocess -i output/predictions/*

        cd - > /dev/null
    else
        echo "Directory $dir does not exist. Skipping commands."
    fi
done

# Deactivate environment
conda deactivate
