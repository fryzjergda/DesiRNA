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
    t_option=$2
    acgu_option=$3
    p_option=$4
    tmin_option=$5
    tmax_option=$6
    sf_option=$7
    nd_option=$8
    acgu_content_option=$9
    oligo_option=${10}
    dimer_option=${11}
    tm_option=${12}
    tm_perc_max_option=${13}
    tm_perc_min_option=${14}
    re_seq_option=${15}

    # Navigate to the project root directory
    cd ..

    # Run the Python program with specified options	
    
    echo ${BASE_DIR}/DesiRNA.py -f "$input_path" -t 2 -acgu "$acgu_option" -p "$p_option" -tmin "$tmin_option" -tmax "$tmax_option" -sf "$sf_option" -nd "$nd_option" -acgu_content "$acgu_content_option" -o "$oligo_option" -d "$dimer_option" -tm "$tm_option" -tm_perc_max "$tm_perc_max_option" -tm_perc_min "$tm_perc_min_option" -re_seq "$re_seq_option"
    python3 ${BASE_DIR}/DesiRNA.py -f "$input_path" -t 2 -acgu "$acgu_option" -p "$p_option" -tmin "$tmin_option" -tmax "$tmax_option" -sf "$sf_option" -nd "$nd_option" -acgu_content "$acgu_content_option" -o "$oligo_option" -d "$dimer_option" -tm "$tm_option" -tm_perc_max "$tm_perc_max_option" -tm_perc_min "$tm_perc_min_option" -re_seq "$re_seq_option"

    # Find the latest directory matching the pattern
    output_dir=$(find_latest_directory "${input_file%.txt}*")
#    echo $(find_latest_directory "${input_file%.txt}*")
    output_stats_file=$(find "$output_dir" -name "*stats")
#    echo $(find "$output_dir" -name "*stats")

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
run_test "Standard_design_input.txt" "2" "off" "1999" "10" "150" "Ed-Epf:1.0" "off" "" "off" "off" "on" "0.7" "0.0" "same"
run_test "Alternative_structures_design_input.txt" "2" "off" "1999" "10" "150" "Ed-Epf:1.0" "off" "" "off" "off" "on" "0.7" "0.0" "same"
run_test "Homodimer_design_input.txt" "10" "off" "1999" "10" "150" "Ed-Epf:1.0" "off" "" "off" "on" "on" "0.7" "0.0" "same"
run_test "RNA_RNA_complex_design_input.txt" "2" "off" "1999" "10" "150" "Ed-Epf:1.0" "off" "" "off" "off" "on" "0.7" "0.0" "same"
run_test "Pseudoknot_design_input.txt" "2" "off" "1999" "10" "150" "Ed-Epf:1.0" "off" "" "off" "off" "on" "0.7" "0.0" "same"
run_test "Seed_sequence_design_input.txt" "2" "off" "1999" "10" "150" "Ed-Epf:1.0" "off" "" "off" "off" "on" "0.7" "0.0" "same"
run_test "Standard_design_input.txt" "2" "on" "1999" "10" "150" "Ed-Epf:1.0" "off" "" "off" "off" "on" "0.7" "0.0" "same"
run_test "Standard_design_input.txt" "2" "off" "2004" "10" "150" "Ed-Epf:1.0" "off" "" "off" "off" "on" "0.7" "0.0" "same"
run_test "Standard_design_input.txt" "2" "off" "1999" "20" "160" "Ed-Epf:1.0" "off" "" "off" "off" "on" "0.7" "0.0" "same"
run_test "Standard_design_input.txt" "2" "off" "1999" "10" "150" "Ed-Epf:0.9,1-MCC:0.05,sln_Epf:0.01,Ed-MFE:0.01,1-precision:0.01,1-recall:0.01,Edef:0.01" "off" "" "off" "off" "on" "0.7" "0.0" "same"
run_test "Standard_design_input.txt" "2" "on" "1999" "10" "150" "Ed-Epf:1.0" "off" "10,40,40,10" "off" "off" "on" "0.7" "0.0" "same"
run_test "Standard_design_input.txt" "10" "off" "1999" "10" "150" "Ed-Epf:1.0" "off" "" "on" "off" "on" "0.7" "0.0" "same"
run_test "Standard_design_input.txt" "2" "off" "1999" "10" "150" "Ed-Epf:1.0" "off" "" "off" "off" "on" "0.8" "0.1" "same"
run_test "Standard_design_input.txt" "2" "off" "1999" "10" "150" "Ed-Epf:1.0" "off" "" "off" "off" "on" "0.7" "0.0" "different"



echo "All tests passed!"
