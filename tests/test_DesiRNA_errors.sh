#!/bin/bash

BASE_DIR=$(dirname "$0")/..
#echo "Base directory: $BASE_DIR"

# Function to find the latest created directory
function find_latest_directory() {
    pattern=$1
    latest_dir=$(ls -dt $pattern*/ | head -n 1)
    echo $latest_dir
}


function run_test {
    command=$1
    expected_message=$2
    generates_dir=$3 # Flag indicating if this test generates an output directory

    echo "Running test: $command"
    DesiRNA_output=$($command 2>&1)
    if [[ $DesiRNA_output == *"$expected_message"* ]] && [[ $DesiRNA_output != *"Exception occurred"* ]]; then
        echo "Test PASSED"
    else
        echo "Test FAILED"
        echo "Output: $DesiRNA_output"
        failed_tests=true
    fi
    echo ""
    
}

failed_tests=false

# Test for -h option
run_test "python3 ${BASE_DIR}/DesiRNA.py -h" "Standard Options:" "no"

# Test for -H option
run_test "python3 ${BASE_DIR}/DesiRNA.py -H" "Advanced Options:" "no"

# Test for input file with not allowed characters in structures
run_test "python3 ${BASE_DIR}/DesiRNA.py -f ${BASE_DIR}/tests/Error_inputs/Error_input_sec_struct_character.txt" "Not allowed characters in structures. Check input file." "yes"

# Test for input file with no opening bracket
run_test "python3 ${BASE_DIR}/DesiRNA.py -f ${BASE_DIR}/tests/Error_inputs/Error_input_no_opening_bracket.txt" "There is no opening bracket for nt position" "yes"

# Test for input file with no closing bracket
run_test "python3 ${BASE_DIR}/DesiRNA.py -f ${BASE_DIR}/tests/Error_inputs/Error_input_no_closing_bracket.txt" "There is no closing bracket for nt position" "yes"

# Test for input file with not allowed characters in sequence restraints
run_test "python3 ${BASE_DIR}/DesiRNA.py -f ${BASE_DIR}/tests/Error_inputs/Error_input_sequence_character.txt" "Not allowed characters in sequence restraints. Check input file." "yes"

# Test for input file with different lengths of secondary structure and sequence restraints
run_test "python3 ${BASE_DIR}/DesiRNA.py -f ${BASE_DIR}/tests/Error_inputs/Error_input_different_length.txt" "Secondary structure and sequence restraints are of different length. Check input file." "yes"

# Test for input file that cannot mutate sequence due to constraints
run_test "python3 ${BASE_DIR}/DesiRNA.py -f ${BASE_DIR}/tests/Error_inputs/Error_input_cannot_mutate.txt" "Cannot mutate the sequence due to sequence constraints." "yes"

# Test for input file with wrong sequence restraints
run_test "python3 ${BASE_DIR}/DesiRNA.py -f ${BASE_DIR}/tests/Error_inputs/Error_input_wrong_seq_restraints.txt" "cannot pair with nucleotide" "yes"

# Test for input file with too many structures
run_test "python3 ${BASE_DIR}/DesiRNA.py -f ${BASE_DIR}/tests/Error_inputs/Error_input_too_much_structures.txt" "Too much structures in the input. Can only design RNA complexes of max two sequences." "no"

# Test for invalid scoring function option
run_test "python3 ${BASE_DIR}/DesiRNA.py -f ${BASE_DIR}/tests/Error_inputs/Error_input_wrong_seq_restraints.txt -sf EK:0.1" "EK is not an available option for scoring function. Check your command." "no"

# Test for incorrect ACGU content sum
run_test "python3 ${BASE_DIR}/DesiRNA.py -f ${BASE_DIR}/tests/Error_inputs/Error_input_wrong_seq_restraints.txt -acgu on -acgu_content 20,40,40,10" "The ACGU content should sum up to 100, check your command." "no"

rm tmp_log.txt
rm -rf *_R10_e100*

if [ "$failed_tests" = false ]; then
    echo "All tests passed!"
else
    echo "Some tests failed!"
fi