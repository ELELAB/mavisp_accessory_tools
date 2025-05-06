#!/bin/bash

source /usr/local/envs/py310/bin/activate

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# --- Argument parsing ---
INPUT_FILE=""
OUTPUT_PREFIX=""
FILTER_FILE=""


usage() {
  echo "Usage: $0 -i <input_file> [-o <output_prefix>] [-f <filter_file>]"
  exit 1
}

while getopts ":i:o:f:h" opt; do
  case $opt in
    i) INPUT_FILE="$OPTARG" ;;
    o) OUTPUT_PREFIX="$OPTARG" ;;
    f) FILTER_FILE="$OPTARG" ;;
    h) usage ;;
    \?) echo "Invalid option: -$OPTARG" >&2; usage ;;
    :) echo "Option -$OPTARG requires an argument." >&2; usage ;;
  esac
done

# Check required input
if [[ -z "$INPUT_FILE" ]]; then
  echo "Error: -i <input_file> is required."
  usage
fi

# Default output prefix if not set
if [[ -z "$OUTPUT_PREFIX" ]]; then
  OUTPUT_PREFIX="summary_output"
fi


OUTPUT_FILE="${OUTPUT_PREFIX}.csv"
UNIQUE_FILE="${OUTPUT_PREFIX}_unique_heteroatoms.txt"
FILTERED_OUTPUT_FILE="${OUTPUT_PREFIX}_filtered.csv"


# Clean previous outputs
rm -f "$OUTPUT_FILE" "$UNIQUE_FILE" "$FILTERED_OUTPUT_FILE"
echo "protein,heteroatoms" > "$OUTPUT_FILE"
if [[ -n "$FILTER_FILE" ]]; then
  echo "protein,cofactors" > "$FILTERED_OUTPUT_FILE"
fi

declare -A UNIQUE_COF
declare -A VALID_COF

# Load filter list (if provided)
if [[ -n "$FILTER_FILE" ]]; then
  if [[ ! -f "$FILTER_FILE" ]]; then
    echo "Error: Filter file '$FILTER_FILE' not found."
    exit 1
  fi
  while IFS= read -r lid; do
    [[ -n "$lid" ]] && VALID_COF["$(echo "$lid" | tr '[:lower:]' '[:upper:]')"]=1
  done < "$FILTER_FILE"
fi


# Process each UniProt ID
while IFS= read -r upid; do
  echo "Processing $upid"
  cofactors=$(python "$SCRIPT_DIR/find_cofactor.py" -u "$upid" 2>/dev/null)

  [[ -z "$cofactors" ]] && cofactors="${upid},Error"
  echo "$cofactors" >> "$OUTPUT_FILE"

  # Extract cofactor list
  line_cofactors=$(echo "$cofactors" | cut -d',' -f2)

  if [[ "$line_cofactors" != "None" && "$line_cofactors" != "Error" ]]; then
    for c in $line_cofactors; do
      UNIQUE_COF["$c"]=1
    done
  fi

# Make filtered output 
  if [[ -n "$FILTER_FILE" ]]; then
    if [[ "$line_cofactors" == "None" || "$line_cofactors" == "Error" ]]; then
      echo "$upid,None" >> "$FILTERED_OUTPUT_FILE"
      continue
    fi

    filtered_hits=()
    for c in $line_cofactors; do
      c_upper=$(echo "$c" | tr '[:lower:]' '[:upper:]')
      if [[ -n "${VALID_COF[$c_upper]}" ]]; then
        filtered_hits+=("$c_upper")
      fi
    done

    if [[ ${#filtered_hits[@]} -eq 0 ]]; then
      echo "$upid,None" >> "$FILTERED_OUTPUT_FILE"
    else
      echo "$upid,${filtered_hits[*]}" >> "$FILTERED_OUTPUT_FILE"
    fi
  fi


done < "$INPUT_FILE"

# Write unique heteroatoms
for cofactor in "${!UNIQUE_COF[@]}"; do
  echo "$cofactor"
done | sort > "$UNIQUE_FILE"

