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

class _FingerprintCalculator:
    ''' Calculate the fingerprint for a molecule, given the fingerprint type
    Parameters: 
        mol (rdkit.Chem.rdchem.Mol) : RdKit mol object (None if invalid smile string smi)
        fp_type (string)            :Fingerprint type  (choices: AP/PHCO/BPF,BTF,PAT,ECFP4,ECFP6,FCFP4,FCFP6)  
    Returns:
        RDKit fingerprint object
    '''

    def get_fingerprint(self, mol: Mol, fp_type: str):
        method_name = 'get_' + fp_type
        method = getattr(self, method_name)
        if method is None:
            raise Exception(f'{fp_type} is not a supported fingerprint type.')
        return method(mol)

    def get_AP(self, mol: Mol):
        return AllChem.GetAtomPairFingerprint(mol, maxLength=10)

    def get_PHCO(self, mol: Mol):
        return Generate.Gen2DFingerprint(mol, Gobbi_Pharm2D.factory)

    def get_BPF(self, mol: Mol):
        return GetBPFingerprint(mol)

    def get_BTF(self, mol: Mol):
        return GetBTFingerprint(mol)

    def get_PATH(self, mol: Mol):
        return AllChem.RDKFingerprint(mol)

    def get_ECFP4(self, mol: Mol):
        return AllChem.GetMorganFingerprint(mol, 2)

    def get_ECFP6(self, mol: Mol):
        return AllChem.GetMorganFingerprint(mol, 3)

    def get_FCFP4(self, mol: Mol):
        return AllChem.GetMorganFingerprint(mol, 2, useFeatures=True)

    def get_FCFP6(self, mol: Mol):
        return AllChem.GetMorganFingerprint(mol, 3, useFeatures=True)

class StonedSelfie(object):

    @staticmethod
    def get_fingerprint(mol: Mol, fp_type: str):
        ''' Fingerprint getter method. Fingerprint is returned after using object of 
            class '_FingerprintCalculator'
            
        Parameters: 
            mol (rdkit.Chem.rdchem.Mol) : RdKit mol object (None if invalid smile string smi)
            fp_type (string)            :Fingerprint type  (choices: AP/PHCO/BPF,BTF,PAT,ECFP4,ECFP6,FCFP4,FCFP6)  
        Returns:
            RDKit fingerprint object
            
        '''
        return _FingerprintCalculator().get_fingerprint(mol=mol, fp_type=fp_type)

    @staticmethod
    def randomize_smiles(mol):
        '''Returns a random (dearomatized) SMILES given an rdkit mol object of a molecule.

        Parameters:
        mol (rdkit.Chem.rdchem.Mol) :  RdKit mol object (None if invalid smile string smi)
        
        Returns:
        mol (rdkit.Chem.rdchem.Mol) : RdKit mol object  (None if invalid smile string smi)
        '''
        if not mol:
            return None

        Chem.Kekulize(mol)
        
        return rdkit.Chem.MolToSmiles(mol, canonical=False, doRandom=True, isomericSmiles=False,  kekuleSmiles=True)

    @staticmethod
    def get_fp_scores(smiles_back, target_smi, fp_type): 
        '''Calculate the Tanimoto fingerprint (using fp_type fingerint) similarity between a list 
        of SMILES and a known target structure (target_smi). 
        
        Parameters:
        smiles_back   (list) : A list of valid SMILES strings 
        target_smi (string)  : A valid SMILES string. Each smile in 'smiles_back' will be compared to this stucture
        fp_type (string)     : Type of fingerprint  (choices: AP/PHCO/BPF,BTF,PAT,ECFP4,ECFP6,FCFP4,FCFP6) 
        
        Returns: 
        smiles_back_scores (list of floats) : List of fingerprint similarities
        '''
        smiles_back_scores = []
        target    = Chem.MolFromSmiles(target_smi)

        fp_target = StonedSelfie.get_fingerprint(target, fp_type)

        for item in smiles_back: 
            mol    = Chem.MolFromSmiles(item)
            fp_mol = StonedSelfie.get_fingerprint(mol, fp_type)
            score  = TanimotoSimilarity(fp_mol, fp_target)
            smiles_back_scores.append(score)
        return smiles_back_scores

    @staticmethod
    def sanitize_smiles(smi):
        '''Return a canonical smile representation of smi
        
        Parameters:
        smi (string) : smile string to be canonicalized 
        
        Returns:
        mol (rdkit.Chem.rdchem.Mol) : RdKit mol object                          (None if invalid smile string smi)
        smi_canon (string)          : Canonicalized smile representation of smi (None if invalid smile string smi)
        conversion_successful (bool): True/False to indicate if conversion was  successful 
        '''
        try:
            mol = smi2mol(smi, sanitize=True)
            smi_canon = mol2smi(mol, isomericSmiles=False, canonical=True)
            return (mol, smi_canon, True)
        except:
            return (None, None, False)

    @staticmethod
    def get_mutated_SELFIES(selfies_ls, num_mutations): 
        ''' Mutate all the SELFIES in 'selfies_ls' 'num_mutations' number of times. 
        
        Parameters:
        selfies_ls   (list)  : A list of SELFIES 
        num_mutations (int)  : number of mutations to perform on each SELFIES within 'selfies_ls'
        
        Returns:
        selfies_ls   (list)  : A list of mutated SELFIES
        
        '''
        for _ in range(num_mutations): 
            selfie_ls_mut_ls = []
            for str_ in selfies_ls: 
                
                str_chars = StonedSelfie.get_selfie_chars(str_)
                max_molecules_len = len(str_chars) + num_mutations
                #print("max_molecules_len: ", max_molecules_len)
                selfie_mutated, _ = StonedSelfie.mutate_selfie(str_, max_molecules_len)
                selfie_ls_mut_ls.append(selfie_mutated)
            
            selfies_ls = selfie_ls_mut_ls.copy()
        return selfies_ls

    @staticmethod
    def get_selfie_chars(selfie):
        '''Obtain a list of all selfie characters in string selfie
        
        Parameters: 
        selfie (string) : A selfie string - representing a molecule 
        
        Example: 
        >>> get_selfie_chars('[C][=C][C][=C][C][=C][Ring1][Branch1_1]')
        ['[C]', '[=C]', '[C]', '[=C]', '[C]', '[=C]', '[Ring1]', '[Branch1_1]']
        
        Returns:
        chars_selfie: list of selfie characters present in molecule selfie
        '''
        chars_selfie = [] # A list of all SELFIE sybols from string selfie
        while selfie != '':
            chars_selfie.append(selfie[selfie.find('['): selfie.find(']')+1])
            selfie = selfie[selfie.find(']')+1:]
        return chars_selfie

    @staticmethod
    def mutate_selfie(selfie, max_molecules_len, write_fail_cases=False):
        '''Return a mutated selfie string (only one mutation on slefie is performed)
        
        Mutations are done until a valid molecule is obtained 
        Rules of mutation: With a 33.3% propbabily, either: 
            1. Add a random SELFIE character in the string
            2. Replace a random SELFIE character with another
            3. Delete a random character
        
        Parameters:
        selfie            (string)  : SELFIE string to be mutated 
        max_molecules_len (int)     : Mutations of SELFIE string are allowed up to this length
        write_fail_cases  (bool)    : If true, failed mutations are recorded in "selfie_failure_cases.txt"
        
        Returns:
        selfie_mutated    (string)  : Mutated SELFIE string
        smiles_canon      (string)  : canonical smile of mutated SELFIE string
        '''
        valid=False
        fail_counter = 0
        chars_selfie = StonedSelfie.get_selfie_chars(selfie)
        
        while not valid:
            fail_counter += 1
                    
            alphabet = list(selfies.get_semantic_robust_alphabet()) # 34 SELFIE characters 

            choice_ls = [1, 2, 3] # 1=Insert; 2=Replace; 3=Delete
            random_choice = np.random.choice(choice_ls, 1)[0]
            
            # Insert a character in a Random Location
            if random_choice == 1: 
                random_index = np.random.randint(len(chars_selfie)+1)
                random_character = np.random.choice(alphabet, size=1)[0]
                
                selfie_mutated_chars = chars_selfie[:random_index] + [random_character] + chars_selfie[random_index:]

            # Replace a random character 
            elif random_choice == 2:                         
                random_index = np.random.randint(len(chars_selfie))
                random_character = np.random.choice(alphabet, size=1)[0]
                if random_index == 0:
                    selfie_mutated_chars = [random_character] + chars_selfie[random_index+1:]
                else:
                    selfie_mutated_chars = chars_selfie[:random_index] + [random_character] + chars_selfie[random_index+1:]
                    
            # Delete a random character
            elif random_choice == 3: 
                random_index = np.random.randint(len(chars_selfie))
                if random_index == 0:
                    selfie_mutated_chars = chars_selfie[random_index+1:]
                else:
                    selfie_mutated_chars = chars_selfie[:random_index] + chars_selfie[random_index+1:]
                    
            else: 
                raise Exception('Invalid Operation trying to be performed')

            selfie_mutated = "".join(x for x in selfie_mutated_chars)
            sf = "".join(x for x in chars_selfie)
            
            try:
                smiles = decoder(selfie_mutated)
                mol, smiles_canon, done = StonedSelfie.sanitize_smiles(smiles)
                if len(selfie_mutated_chars) > max_molecules_len or smiles_canon=="":
                    done = False
                if done:
                    valid = True
                else:
                    valid = False
            except:
                valid=False
                if fail_counter > 1 and write_fail_cases == True:
                    f = open("selfie_failure_cases.txt", "a+")
                    f.write('Tried to mutate SELFIE: '+str(sf)+' To Obtain: '+str(selfie_mutated) + '\n')
                    f.close()
        
        return (selfie_mutated, smiles_canon)
    


