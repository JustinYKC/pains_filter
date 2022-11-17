import argparse
from pathlib import Path
from painsFilter import PainsFilter


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter molecules with PAINS")
    subparsers = parser.add_subparsers(dest= 'subcmd', help='subcommands', metavar='SUBCOMMAND')
    subparsers.required = True

    parser_f1 = subparsers.add_parser("fliter_pains", help="Apply PAINS filter to query compounds")
    parser_f1.add_argument("-q_mol", "--query_file", help="The path of the input file that contains query compounds to be performed with PAINS filter", dest= "infile", required=True)
    parser_f1.add_argument("-out", "--output_file", help="The path of the output file that contains result of Non-PAINS compounds", dest= "outfile", required=True)
    
    args = parser.parse_args()
    pains = PainsFilter()

    if args.subcmd == "fliter_pains":
        pains.filter_by_pains(Path(args.infile), Path(args.outfile))