import csv
import time 
import selfies
import rdkit
import random
import numpy as np
import random
from rdkit import Chem
from selfies import encoder, decoder
from rdkit.Chem import MolFromSmiles as smi2mol
from rdkit.Chem import AllChem
from rdkit.DataStructs.cDataStructs import TanimotoSimilarity
from rdkit.Chem import Mol
from rdkit.Chem.AtomPairs.Sheridan import GetBPFingerprint, GetBTFingerprint
from rdkit.Chem.Pharm2D import Generate, Gobbi_Pharm2D
from rdkit.Chem import Draw

from rdkit.Chem import MolToSmiles as mol2smi
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

from stonedselfies import StonedSelfie



if __name__ == '__main__':

    csv_filename = "stoned_selfies_output.csv"



    num_random_samples = 1000     
    num_mutation_ls    = [1, 2, 3, 4, 5]

    with open(csv_filename, mode='w', newline='') as file:
        writer = csv.writer(file)

        with open('smiles_atomwise.txt') as smiles_from_txt:
            for smi in smiles_from_txt:


                fp_type = 'ECFP4' # question what is this?

                total_time = time.time()

                print("\n\nWorking with smile:", smi)
                mol = Chem.MolFromSmiles(smi)
                if mol == None: 
                    raise Exception('Invalid starting structure encountered')
                
                start_time = time.time()
                randomized_smile_orderings  = [StonedSelfie.randomize_smiles(mol) for _ in range(num_random_samples)]

                # Convert all the molecules to SELFIES
                selfies_ls = [encoder(x) for x in randomized_smile_orderings]

                all_smiles_collect = []
                all_smiles_collect_broken = []

                start_time = time.time()
                for num_mutations in num_mutation_ls: 
                    # Mutate the SELFIES: 
                    selfies_mut = StonedSelfie.get_mutated_SELFIES(selfies_ls.copy(), num_mutations=num_mutations)

                    # Convert back to SMILES: 
                    smiles_back = [decoder(x) for x in selfies_mut]
                    all_smiles_collect = all_smiles_collect + smiles_back
                    all_smiles_collect_broken.append(smiles_back)


                print('Mutation obtainment time (back to smiles): ', time.time()-start_time)


                # Work on:  all_smiles_collect
                start_time = time.time()
                canon_smi_ls = []
                for item in all_smiles_collect: 
                    mol, smi_canon, did_convert = StonedSelfie.sanitize_smiles(item)
                    if mol == None or smi_canon == '' or did_convert == False: 
                        raise Exception('Invalid smile string found')
                    canon_smi_ls.append(smi_canon)
                canon_smi_ls        = list(set(canon_smi_ls))
                print('Unique mutated structure obtainment time: ', time.time()-start_time)

                start_time = time.time()
                canon_smi_ls_scores = StonedSelfie.get_fp_scores(canon_smi_ls, target_smi=smi, fp_type=fp_type)
                print('Fingerprint calculation time: ', time.time()-start_time)
                print('Total time: ', time.time()-total_time)

                # Molecules with fingerprint similarity > 0.8
                indices_thresh_8 = [i for i,x in enumerate(canon_smi_ls_scores) if x > 0.8]
                mols_8 = [Chem.MolFromSmiles(canon_smi_ls[idx]) for idx in indices_thresh_8]

                # Molecules with fingerprint similarity > 0.6
                indices_thresh_6 = [i for i,x in enumerate(canon_smi_ls_scores) if x > 0.6 and x < 0.8]
                mols_6 = [Chem.MolFromSmiles(canon_smi_ls[idx]) for idx in indices_thresh_6]

                # Molecules with fingerprint similarity > 0.4
                indices_thresh_4 = [i for i,x in enumerate(canon_smi_ls_scores) if x > 0.4 and x < 0.6]
                mols_4 = [Chem.MolFromSmiles(canon_smi_ls[idx]) for idx in indices_thresh_4]



                #for string in mols_8:
                #    Write each string and its length as a new row
                print("Writing to csv file") #[smi]+[canon_smi_ls[idx] for idx in indices_thresh_8])
                writer.writerow([smi.strip()]+[canon_smi_ls[idx] for idx in indices_thresh_8])
                #writer.writerow([1,2])
                file.flush()