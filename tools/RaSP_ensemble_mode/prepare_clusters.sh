#!/bin/bash

# Usage: ./prepare_clusters.sh <input_directory> <num_clusters>
# Example: ./prepare_clusters.sh /data/user/shared_projects/mavisp/TPMT/simulations_analysis/free/AF2_18-245/replicate1/CHARMM36/md/6.clustering/ 4

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_directory> <num_clusters>"
    exit 1
fi

INPUT_DIR=$1
NUM_CLUSTERS=$2

# Loop over clusters
for (( i=1; i<=$NUM_CLUSTERS; i++ ))
do
    CLUSTER_FILE="$INPUT_DIR/clusters_$i.pdb"
    DEST_DIR="frame$i"

    # Check if cluster file exists
    if [ ! -f "$CLUSTER_FILE" ]; then
        echo "Error: $CLUSTER_FILE not found. Stopping."
        exit 1
    fi

    # Create destination directory if it doesn't exist
    mkdir -p "$DEST_DIR"

    # Copy file to destination directory
    cp "$CLUSTER_FILE" "$DEST_DIR/"
    echo "Copied $CLUSTER_FILE to $DEST_DIR/"
done

