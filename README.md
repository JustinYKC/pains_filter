# pains_filter
A small tool to screen given virtual chemical libraries by the PAINS filter based on SMARTS.

# Requirements
pains_filter requires the following packages:
- rdkit (>=2019.09.x.x)

# Usage
### Parameters
*   `fliter_pains` *(Required)* Apply PAINS filter to query compounds.
    *   `-q_mol, --query_file` *(Required)* The path of the input file that contains query compounds to be performed with PAINS filter.
    *   `-mod_file, --model_file` *(Required)* The path of the output file that contains result of Non-PAINS compounds.