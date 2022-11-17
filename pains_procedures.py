import argparse
from painsFilter import painsFilter


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter molecules with PAINS")
    subparsers = parser.add_subparsers(dest= 'subcmd', help='subcommands', metavar='SUBCOMMAND')
    subparsers.required = True

    parser_f1 = subparsers.add_parser("fliter_pains", help="Apply PAINS filter to query compounds")
    parser_f1.add_argument("-q_mol", "--query_file", help="The path of the input file that contains query compounds to be performed with PAINS filter", type=argparse.FileType('r'), dest= "query file", required=True)
    parser_f1.add_argument("-out", "--output_file", help="The path of the output file that contains result of Non-PAINS compounds", type=argparse.FileType('w'), dest= "output file", required=True)
    
    args = parser.parse_args()
    pains = PainsFilter()

    if args.subcmd == "fliter_pains":
        pains.filter_by_pains(input_smiles_path, output_smiles_path)