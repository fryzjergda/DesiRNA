#!/bin/bash

# Function to find the latest created directory
function find_latest_directory() {
    pattern=$1
    latest_dir=$(ls -dt $pattern*/ | head -n 1)
    echo $latest_dir
}

# Function to run the program and check the output
function run_test() {
    # Get the absolute path of the project's root directory
    BASE_DIR=$(git rev-parse --show-toplevel)
    
    input_file=$1
    expected_phrase="Design solved succesfully!"
    input_path="${BASE_DIR}/example_files/inputs/$input_file"
    oligo_option=$2
    p_option=$3
    

    # Navigate to the project root directory
    cd ..

    # Run the Python program with specified options	
    
    echo DesiRNA.py -f "$input_path" -t 2 -o "$oligo_option" -p "$p_option"
    python3 ${BASE_DIR}/DesiRNA.py -f "$input_path" -t 2 -o "$oligo_option" -p "$p_option"

    # Find the latest directory matching the pattern
    output_dir=$(find_latest_directory "${input_file%.txt}*")
    echo $(find_latest_directory "${input_file%.txt}*")
    output_stats_file=$(find "$output_dir" -name "*stats")
    echo $(find "$output_dir" -name "*stats")

    # Check if the output stats file exists and contains the expected phrase
    if [ -f "$output_stats_file" ] && grep -q "$expected_phrase" "$output_stats_file"; then
        echo "Test with input $input_file and options -o $oligo_option -p $p_option PASSED"
        rm -rf "$output_dir" # Clean up the directory after the test
    else
        echo "Test with input $input_file and options -o $oligo_option -p $p_option FAILED"
        exit 1
    fi

    # Navigate back to the test directory
    cd - > /dev/null
}

# Running tests with different inputs and options
run_test "Standard_design_input.txt" "on" "2004"
run_test "Alternative_structures_design_input.txt" "off" "2004"
run_test "Homodimer_design_input.txt" "off" "1999"
run_test "RNA_RNA_complex_design_input.txt" "off" "2004"
run_test "Pseudoknot_design_input.txt" "off" "2004"
run_test "Seed_sequence_design_input.txt" "off" "1999"
# Add more test cases here with different inputs and options...

echo "All tests passed!"
