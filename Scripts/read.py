# -*- coding: utf-8 -*-
"""
Created on Thu Sep 26 09:46:45 2019

@author: hcji
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Jun 25 08:24:41 2019

@author: hcji
"""

import json
import numpy as np
from scipy.sparse import csr_matrix, save_npz
from tqdm import tqdm
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdMolDescriptors import CalcExactMolWt
from pycdk.pycdk import MolFromSmiles, parser_formula, MolToFormula, getMolecularDescriptor
from DeepEI.utils import ms2vec, fp2vec, get_cdk_fingerprints, get_cdk_descriptors

from matchms.importing import load_from_msp
import sys

all_mol = list(load_from_msp(sys.argv[1]))

def read_mol(i):
    name = all_mol[i].get('name')
    smiles = all_mol[i].get('smiles')
    retention = ''
    peakindex = all_mol[i].peaks.mz
    peakintensity = all_mol[i].peaks.intensities
    RI = {}
    RI['SemiStdNP'] = np.nan
    RI['StdNP'] = np.nan
    RI['StdPolar'] = np.nan
    if retention != '':
        retention = retention.split(' ')
        for r in retention:
            if 'SemiStdNP' in r:
                RI['SemiStdNP'] = float(r.split('=')[1].split('/')[0])
            if ('StdNP' in r) and ('Semi' not in r):
                RI['StdNP'] = float(r.split('=')[1].split('/')[0])
            if 'StdPolar' in r:
                RI['StdPolar'] = float(r.split('=')[1].split('/')[0])
    output = {'name': name, 'smiles': smiles, 'RI': RI, 'peakindex': peakindex, 'peakintensity': peakintensity}
    return output


def collect():
    all_smiles = []
    Peak_data = []
    RI_data = []
    Morgan_fp = []
    CDK_fp = []
    CDK_des = []
    MolWt = []
    # for i in tqdm(range(20)):
#    for i in tqdm(range(len(all_mol))):
    for i in range(len(all_mol)):
        try:
            m = read_mol(i)
        except:
            continue
        print(i,'read_mol ok')
        '''
        if  'TMS derivative' in m['name']:
            derive = 1
        else:
            derive = 0
        '''
        try:
            mol = Chem.MolFromSmiles(m['smiles'])
            molwt = CalcExactMolWt(mol)
            if molwt > 2000:
                continue
            smiles = Chem.MolToSmiles(mol)
            # check element
            elements = parser_formula(MolToFormula(MolFromSmiles(smiles)))
            for e in elements:
                if e not in ['C', 'H', 'O', 'N', 'S', 'P', 'Si', 'F', 'Cl', 'Br', 'I']:
                    raise ValueError ('contain uncommon element')
            morgan_fp = np.array(AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=4096))
            print('try fingerprints')
            cdk_fp = get_cdk_fingerprints(smiles)
            # cdk_fp = fp2vec(cdk_fp)
            print('try descriptors')
            cdk_des = np.array(get_cdk_descriptors(smiles))
            # cdk_des = getMolecularDescriptor(MolFromSmiles(smiles)).values()
            # cdk_des  = np.array(list(itertools.chain(*cdk_des)))
            ri = list(m['RI'].values())
            print('call ms2vec')
            peak_vec = ms2vec(m['peakindex'], m['peakintensity'])
            print('OK')
        except:
            continue
        
        all_smiles.append(smiles)
        Peak_data.append(peak_vec)
        RI_data.append(ri)
        Morgan_fp.append(morgan_fp)
        CDK_fp.append(cdk_fp)
        CDK_des.append(cdk_des)
        MolWt.append(molwt)
 
    # save
    np.save('DeepEI/data/retention.npy', np.array(RI_data))
    np.save('DeepEI/data/descriptor.npy', np.array(CDK_des))
    np.save('DeepEI/data/molwt.npy', np.array(MolWt))
    
    Peak_data = csr_matrix(np.array(Peak_data))
    Morgan_fp = csr_matrix(np.array(Morgan_fp))
    CDK_fp = csr_matrix(np.array(CDK_fp))
    save_npz('DeepEI/data/peakvec.npz', Peak_data)
    save_npz('DeepEI/data/morgan.npz', Morgan_fp)
    save_npz('DeepEI/data/fingerprints.npz', CDK_fp)
    
    with open('DeepEI/data/all_smiles.json', 'w') as t:
        json.dump(all_smiles, t)


if __name__ == '__main__':
    collect()
