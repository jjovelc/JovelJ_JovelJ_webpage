#!/bin/bash

# This script will remove columns from a table
# if their average abundance is below a threshold ($3)
#
# USAGE: bash filter_columns.sh <infile> <outfile> <threshold>
#

inputFile=$1
outputFile=$2
threshold=$3

# Count the number of rows
numRows=$(wc -l < "$inputFile")

# Find columns to keep, those with average 10 or more
keepCols=$(awk -F'\t' -v numRows="$numRows" '
{
    for (i=2; i<=NF; i++) {
        sum[i] += $i
    }
}
END {
    printf "1,"; # Always keep the first column
    for (i=2; i<=NF; i++) {
        if (sum[i] / numRows >= 10) {
            printf i ","
        }
    }
}' "$inputFile" | sed 's/,$//')  # Remove trailing comma

# Keep only the identified columns
awk -F'\t' -v keep="$keepCols" '{
    split(keep, keepArray, ",");
    output = "";
    for (i=1; i<=length(keepArray); i++) {
        colIndex = keepArray[i] + 0;
        output = output $colIndex;
        if (i != length(keepArray)) {
            output = output "\t";
        }
    }
    print output;
}' "$inputFile" > "$outputFile"

echo "Process complete. Filtered data saved to $outputFile."
