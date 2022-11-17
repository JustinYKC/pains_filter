from rdkit import Chem
from rdkit.Chem import AllChem
from pathlib import Path, PurePath
import pandas as pd
import sys

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
        with open("./pains.txt", 'r') as _pains_substructure_smart:
            pains_substructures_smart = _pains_substructure_smart.readlines()
            self._pains_substructure_list = [Chem.MolFromSmarts(k) for k in pains_substructures_smart]
        print (f"The number of substructure fliters for PAINS: {len(self._pains_substructure_list)}")
    def _check_input(self, _query_path:Path):
        """
        A private function used to check the type of the file containing query compounds.

        Parameters
        ----------
        _query_path: (File path object, required)
            The path of the file containing query compounds

 
        """ 
        if _query_path.suffix == ".smi": return "smi"
        elif _query_path.suffix == ".sdf": return "sdf"
        else: 
            print ("The query file is invalid, please use the valid format in terms of 'SDF', or 'SMILES'")
            sys.exit(-1)
    def _pains_filt(self, _mol:Chem, _sub_list:list):
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
        for ss in _sub_list:
            if _mol.HasSubstructMatch(ss):
                print ("PAIN found")
                return True
        return False
    def filter_by_pains(self, query_path:Path, result_path:Path, delimiter:str="\t", smiles_column:int=0, name_column:int=1):
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
            suppl = Chem.SmilesMolSupplier(str(query_path), delimiter, smiles_column, name_column)
        else:
            suppl = Chem.SDMolSupplier(str(query_path))
        _cpd_without_pains_dict={}
        print (len(suppl))
        for index, mol in enumerate(suppl):
            if mol:
                mol_addH = Chem.AddHs(mol) #Add hydrogen
                if (self._pains_filt(mol_addH, self._pains_substructure_list) == False):
                    _cpd_without_pains_dict[mol.GetProp("_Name")] = Chem.MolToSmiles(mol).strip()

        print (f"Total number of cpds without PAINS: {len(_cpd_without_pains_dict)}")
        with open(str(result_path), 'w') as output:
            for _cpd_name, _cpd_smiles in _cpd_without_pains_dict.items():
                output.write(f"{_cpd_smiles.strip()}\t{_cpd_name.strip()}\n")

if __name__ == "__main__":
    work_dir = Path("/home/justin/Vichem_paper/kinase_group/data/for_cpd_diversity_new")
    lib_sheet_df = pd.read_csv(PurePath(work_dir, Path("library_sheet.csv")), sep="\t", header=0)
    print (lib_sheet_df)
    for index, row in lib_sheet_df.iterrows():
        dir_name = row["Name"] 
        abbra_name = row["Abbreviation"]
        input_smiles_path = PurePath(work_dir, Path(dir_name), Path(f"standardised_smiles_{abbra_name}_withoutReplicates.smi"))
        output_smiles_path = PurePath(work_dir, Path(dir_name), Path(f"standardised_smiles_{abbra_name}_withoutReplicates_withoutPAINS.smi"))
        pains = PainsFilter()
        pains.filter_by_pains(input_smiles_path, output_smiles_path)
        #pains.filter_by_pains(  Path("/home/justin/Vichem_paper/kinase_group/data/for_cpd_diversity_new/APExBIO_DiscoveryProbe_Kinase_Inhibitor_Library/standardised_smiles_APExBIO-Ki_withoutReplicates.smi"),
        #                        Path("/home/justin/Vichem_paper/kinase_group/data/for_cpd_diversity_new/APExBIO_DiscoveryProbe_Kinase_Inhibitor_Library/standardised_smiles_APExBIO-Ki_withoutReplicates_withoutPAINS.smi"))