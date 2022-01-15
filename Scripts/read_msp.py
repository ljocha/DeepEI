from DeepEI.utils import ms2vec, get_cdk_fingerprints, get_cdk_descriptors
from matchms.importing import load_from_msp
import json
import numpy as np
from scipy.sparse import csr_matrix, save_npz
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdMolDescriptors import CalcExactMolWt

# incompatible with DeepEI/utils.py
#from pycdk.pycdk import MolFromSmiles, parser_formula, MolToFormula

from concurrent.futures import ProcessPoolExecutor
import os
from argparse import ArgumentParser

p = ArgumentParser()
p.add_argument('--ncores','-n',type=int,help='number of cores',default=1)
p.add_argument('--dest','-d',type=str,help='destination directory',default='.')
p.add_argument('infile',type=str,help='input file')

args = p.parse_args()
file_msp = args.infile
ncores = args.ncores
dest = args.dest

if not os.path.isdir(dest):
	print(f"{dest} does not exist")
	exit(1)

def process_mol(nm):
	n,m = nm
	try:
		osmiles = m.get('smiles')
		mol = Chem.MolFromSmiles(osmiles)
		name = m.get('name')
		peakindex = m.peaks.mz
		peakintensity = m.peaks.intensities

		molwt = CalcExactMolWt(mol)
		if molwt > 2000:
			return {}
		smiles = Chem.MolToSmiles(mol)
# XXX: pycdk
#		elements = parser_formula(MolToFormula(MolFromSmiles(smiles)))
#		for e in elements:
#			if e not in ['C', 'H', 'O', 'N', 'S', 'P', 'Si', 'F', 'Cl', 'Br', 'I']:
#				print(f"{osmiles}: uncommon element {e}, skipping")
#				return {}
		morgan_fp = np.array(AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=4096))
		cdk_fp = get_cdk_fingerprints(smiles)
		cdk_des = np.array(get_cdk_descriptors(smiles))
# XXX
#		ri = list(m['RI'].values())
		peak_vec = ms2vec(peakindex,peakintensity)

		print(f"{n}:{osmiles}: done")
		return { 
			'smiles': smiles,
			'name': name,
			'peak_vec': peak_vec,
#			'ri': ri,
			'morgan_fp': morgan_fp,
			'cdk_fp': cdk_fp,
			'cdk_des': cdk_des,
			'molwt': molwt,
		}
	except BaseException as e:
		print(f"{osmiles}: {e}")
		return {}

print(f"Loading {file_msp}...")
all_mol = load_from_msp(file_msp)
print("done")

with ProcessPoolExecutor(max_workers=ncores) as pool:
	all_output = pool.map(process_mol, enumerate(all_mol))

# filter out empty entries
all_output = list(filter(lambda x: x,all_output))

all_smiles = list(map(lambda x: x['smiles'], all_output))
Peak_data = np.array(list(map(lambda x: x['peak_vec'], all_output)))
# RI_data = map(lambda x: x['smiles'], all_output)
Morgan_fp = np.array(list(map(lambda x: x['morgan_fp'], all_output)))
CDK_fp = np.array(list(map(lambda x: x['cdk_fp'], all_output)))
CDK_des = np.array(list(map(lambda x: x['cdk_des'], all_output)))
MolWt = np.array(list(map(lambda x: x['molwt'], all_output)))

print("writing output ...")
os.chdir(dest)

# np.save('retention.npy', np.array(RI_data))
np.save('descriptor.npy', CDK_des)
np.save('molwt.npy', MolWt)

Peak_data = csr_matrix(Peak_data)
Morgan_fp = csr_matrix(Morgan_fp)
CDK_fp = csr_matrix(CDK_fp)

save_npz('peakvec.npz', Peak_data)
save_npz('morgan.npz', Morgan_fp)
save_npz('fingerprints.npz', CDK_fp)

with open('all_smiles.json', 'w') as t:
	json.dump(all_smiles, t)

print("done")
