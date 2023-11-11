from rdkit import Chem
from pathlib import Path
import argparse
from tqdm import tqdm


class PainsFilter(object):
    """
    Pain filter class, used to build a pain filter identify PAINS.

    """

    def __init__(self):
        self._define_pains_substructures()

    def _define_pains_substructures(self):
        """
        A private function used to read pre-defined substructures for PAINS from SMRATS: 'pains.txt'.

        """
        with open("./pains.txt", "r") as _pains_substructure_smart:
            pains_substructures_smart = _pains_substructure_smart.readlines()
            self._pains_substructure_list = [
                Chem.MolFromSmarts(k) for k in pains_substructures_smart
            ]
        print(
            f"The number of substructure fliters for PAINS: {len(self._pains_substructure_list)}"
        )

    def _check_input(self, query_path: Path) -> str:
        """
        A private function used to check the type of the file containing query compounds.

        Parameters
        ----------
        _query_path: (File path object, required)
            The path of the file containing query compounds


        """
        if query_path.suffix == ".smi":
            return "smi"
        elif query_path.suffix == ".sdf":
            return "sdf"
        else:
            raise ValueError(
                "The query file is invalid, please use the valid format in terms of '.sdf', or '.smi'"
            )

    def _pains_filt(self, mol: Chem, sub_list: list) -> bool:
        """
        A private function used to check whether a query compound contains PAINS substrucutes.

        Parameters
        ----------
        _mol: (RDkit molecule object, required)
            The RDkit molecule object of a query compound

        _sub_list: (list, required)
            The list has all PAINS substructures

        Returns
        -------
        It returns True if any PAINS substructure found in a query compound, otherwise it returns False.

        """
        for ss in sub_list:
            if mol.HasSubstructMatch(ss):
                return True
        return False

    def filter_by_pains(
        self,
        query_path: Path,
        result_path: Path,
        delimiter: str = "\t",
        smiles_column: int = 0,
        name_column: int = 1,
    ):
        """
        A function performs the check of PAINS for query compounds.

        Parameters
        ----------
        query_path: (File path object, required)
            The path of the file containing query compounds.

        result_path: (File path object, required)
            The path of the output file that contains the result of Non-PAINS compounds including which SMILES and name.

        delimiter: (str, optional)
            The delimiter used to split SMIILES and compound names in a query SMILES file. Default is '\t'.

        smiles_column: (int, optional)
            The number of the column which contains SMILES. Default is 0.

        name_column: (int, optional)
            The number of the column which contains compound name. Default is 1.


        """
        input_file_format = self._check_input(query_path)
        if input_file_format == "smi":
            suppl = Chem.SmilesMolSupplier(
                str(query_path), delimiter, smiles_column, name_column
            )
        else:
            suppl = Chem.SDMolSupplier(str(query_path))
        _cpd_without_pains_dict = {}
        n_pains_found = 0
        print(f"Total number of query cpds for PAINS filter: {len(suppl)}")
        with open(str(result_path), "w") as output:
            for mol in (
                pbar := tqdm(
                    suppl,
                    desc="Filtering pains",
                    bar_format="{l_bar}{bar}{r_bar} PAINS found: {postfix} | Elapsed: {elapsed} | {rate_fmt}",
                    postfix=n_pains_found,
                )
            ):
                if mol:
                    mol_addH = Chem.AddHs(mol)  # Add hydrogen
                    if self._pains_filt(mol_addH, self._pains_substructure_list):
                        n_pains_found += 1
                        pbar.postfix = n_pains_found
                        pbar.update()
                    else:
                        try:
                            mol_title = mol.GetProp("_Name")
                        except:
                            mol_title = ""
                        output.write(f"{Chem.MolToSmiles(mol).strip()} {mol_title}\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter molecules with PAINS")
    subparsers = parser.add_subparsers(
        dest="subcmd", help="subcommands", metavar="SUBCOMMAND"
    )
    subparsers.required = True

    parser_f1 = subparsers.add_parser(
        "filter_pains", help="Apply PAINS filter to query compounds"
    )
    parser_f1.add_argument(
        "-q_mol",
        "--query_file",
        help="The path of the input file that contains query compounds to be performed with PAINS filter",
        dest="infile",
        required=True,
    )
    parser_f1.add_argument(
        "-out",
        "--output_file",
        help="The path of the output file that contains result of Non-PAINS compounds",
        dest="outfile",
        required=True,
    )

    args = parser.parse_args()
    pains = PainsFilter()

    if args.subcmd == "filter_pains":
        pains.filter_by_pains(Path(args.infile), Path(args.outfile))
