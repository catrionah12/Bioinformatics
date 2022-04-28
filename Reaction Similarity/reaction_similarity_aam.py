from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
from rdkit import DataStructs
import pandas as pd
import time
from collections import Counter
from rdkit.Chem import rdFMCS
from sys import argv
import csv
import argparse

class Pair():
    def __init__(self, pair, mol1_fp, mol2_fp):
        self.pair = pair
        self.mol1_fp = mol1_fp
        self.mol2_fp = mol2_fp
    def get_sim(self):
        self.sim = DataStructs.DiceSimilarity(self.mol1_fp, self.mol2_fp)
    def get_ctr_fp(self):
        self.ctr_fp = rc_gen.GetCountFingerprint(self.pair[0], fromAtoms = s_centres) + rc_gen.GetCountFingerprint(self.pair[1], fromAtoms = p_centres) 
    def map_atoms(self):
        mcs = rdFMCS.FindMCS(self.pair, bondCompare=Chem.rdFMCS.BondCompare.CompareAny, ringMatchesRingOnly = False, timeout = 1)
        patt = Chem.MolFromSmarts(mcs.smartsString)
        for mol in self.pair:
            hit_ats = list(mol.GetSubstructMatch(patt))
            for index, atom in enumerate(mol.GetAtoms()):
                if index in hit_ats:
                    if atom.GetAtomMapNum() == 0:
                        atom.SetAtomMapNum(hit_ats.index(index)+1)
        
class Molecule():
    def __init__(self, mol):
        self.mol = mol
    def get_fp(self):
        self.fp = gen.GetCountFingerprint(self.mol)

class Reaction():
    def __init__(self, substrates, products, matches):
        self.substrates = substrates
        self.products = products
        self.matches = matches

class AtomNeighbors():
    def __init__(self, atom):
        self.atom = atom
        self.neighbors = set()
        self.bonds = Counter()
    def add_bonds(self, atom):
        for bond in atom.GetBonds():
            self.bonds.update([bond.GetBondType()])
   
def fp_match(a, b):
    
    num_pairs = min([len(a), len(b)])      # Max number of pairs that can be made
    pairs = []
    for mol in a:
        for mol2 in b:
            pair = Pair([mol.mol,mol2.mol], mol.fp, mol2.fp)     # Make a pair object for every combination of molecules
            pair.get_sim()
            pairs.append(pair)                   
    sorted_pairs = sorted(pairs, key = lambda a:a.sim, reverse = True)   #Sort by similarity
    matches = []
    match_pairs(matches, sorted_pairs, num_pairs)
    return matches

def ctr_fp_match(reac_pairs1, reac_pairs2):
    num_pairs = min([len(reac_pairs1), len(reac_pairs2)])  
    pairs = []
    for match1 in reac_pairs1:                        
        for match2 in reac_pairs2:
            pair = Pair([match1,match2], match1.ctr_fp, match2.ctr_fp)     #Make a pair of matched pairs, with their reaction centre fingerprints
            pair.get_sim()                                                  #Find similarity of the centre fingerprints
            pairs.append(pair)
    sorted_pairs = sorted(pairs, key = lambda a:a.sim, reverse = True) 
    matches = []
    match_pairs(matches, sorted_pairs, num_pairs)
    return matches

def match_pairs(matches, sorted_pairs, num_pairs):
    while len(matches) < num_pairs:
        top = sorted_pairs.pop(0)   #Add most similar pair to matches
        matches.append(top)
        for mol in top.pair:
            sorted_pairs[:] = [i for i in sorted_pairs if mol not in i.pair] #Remove other pairs containing matched molecules

def find_centres(pair, distance):
    
    atoms = []
    s_centres = []
    p_centres = []
    
    for pair_idx, mol in enumerate(pair.pair):                              #Go through each molecule in the pair 
        for atom in mol.GetAtoms():
            if atom.GetAtomMapNum() == 0:  #### include?                    #Add unmapped atoms to reaction centre list
                if pair_idx == 1:
                    p_centres.append(atom.GetIdx())
                else:
                    s_centres.append(atom.GetIdx())
            else:
                visited_atoms = []
                current_atom = AtomNeighbors(atom)
                current_atom = find_neighbours(distance, atom, current_atom, visited_atoms)[1]      #Get the neighbours and bonds of each atom
                if pair_idx == 1:
                    for substrate_atom in atoms:
                        if substrate_atom.atom.GetAtomMapNum() == current_atom.atom.GetAtomMapNum():            #See if neighbours or bonds have changed + add to reaction centre list if they have
                            if substrate_atom.neighbors != current_atom.neighbors or substrate_atom.bonds != current_atom.bonds or 0 in substrate_atom.neighbors or 0 in current_atom.neighbors:
                                s_centres.append(substrate_atom.atom.GetIdx())                          
                                p_centres.append(current_atom.atom.GetIdx())                                     
                else:
                    atoms.append(current_atom) 
    return s_centres, p_centres

def find_neighbours(distance, atom, original_atom, visited_atoms):          
    visited_atoms.append(atom.GetAtomMapNum())
    original_atom.add_bonds(atom)                                      #Add bonds to list
    for neighbour in atom.GetNeighbors():
        if neighbour.GetAtomMapNum() not in visited_atoms:   
            original_atom.neighbors.add(neighbour.GetAtomMapNum())     #Add neighbours to list 
            if distance != 0:
                distance -= 1
                distance = find_neighbours(distance, neighbour, original_atom, visited_atoms)[0]     #Find neighbours of neighbours until distance is reached
    distance +=1
    return distance, original_atom

parser = argparse.ArgumentParser(prog='Reaction Similarity')
parser.add_argument('fingerprint', type = str, help = 'Choose a type of fingerprint from \'m\' (Morgan), \'ap\' (Atom Pair) or \'t\' (topological)')
parser.add_argument('-file', type = str, help = 'Choose a reaction SMARTS file', default = 'SMARTS.csv')
args = parser.parse_args()
fp_type = args.fingerprint
reactions = []
sims = {}
distance = 1

# Make fingerprint generators for whole molecule and reaction centre fingerprints
if fp_type.lower() == 'm':
    gen = rdFingerprintGenerator.GetMorganGenerator(includeChirality = True)
    rc_gen = rdFingerprintGenerator.GetMorganGenerator(includeChirality = True, radius = 1)
elif fp_type.lower() == 'ap':
    gen = rdFingerprintGenerator.GetAtomPairGenerator(includeChirality = True)
    rc_gen = rdFingerprintGenerator.GetAtomPairGenerator(includeChirality = True, maxDistance = 2)
elif fp_type.lower() == 't':
    gen = rdFingerprintGenerator.GetRDKitFPGenerator()
    rc_gen = rdFingerprintGenerator.GetRDKitFPGenerator(maxPath = 2)
else:
    print('Invalid fingerprint chosen')
    exit()

f = open(args.file, 'r')
r = csv.reader(f)
next(r)

# Read in all molecules and generate all fingerprints
for row in r:
    try:
        reac = row[1].strip()
        lhs = reac.split('>>')[0]
        rhs = reac.split('>>')[1]
        substrates = [Molecule(Chem.MolFromSmiles(i)) for (i) in lhs.split('.')]
        products = [Molecule(Chem.MolFromSmiles(i)) for (i) in rhs.split('.')]
        for molecule in substrates:
            molecule.get_fp()
        for molecule in products:
            molecule.get_fp()    
        matches = fp_match(substrates, products)
        for pair in matches:
            pair.map_atoms()
            s_centres, p_centres = find_centres(pair, distance)
            pair.get_ctr_fp()   
        if len(lhs) == 0 or len(rhs) == 0 or len(reac.split('>>'))>2:
            print(f'{row[1]} not a SMILES/SMARTS reaction')
            reactions.append(None)
        else:
            reactions.append(Reaction(substrates, products, matches))
    except IndexError:
        print(f'{row[1]} not a SMILES/SMARTS reaction')
        reactions.append(None)
    
f.close()

# Compare fingerprints for each pair of reactions
for idx, reac1 in enumerate(reactions):
    for idx2, reac2 in enumerate(reactions[idx+1:]):

        if reac1 != None and reac2!= None:
            
            # Reaction centre similarity
            ctr_matches = ctr_fp_match(reac1.matches, reac2.matches)                    # Match up substrate/product pairs from each reaction
            transformation_level_sim = sum([(i).sim for (i) in ctr_matches])/len(ctr_matches)       # Average the similarity of each pair

            #Reaction similarity
            # Find molecules that match between the two reactions - in both directions
            matched_reactants = fp_match(reac1.substrates, reac2.substrates)
            matched_products = fp_match(reac1.products, reac2.products)
            matched_reactants_rev = fp_match(reac1.products, reac2.substrates)
            matched_products_rev = fp_match(reac1.substrates, reac2.products)
            
            #Calculate overall similarity score for each direction
            fwd_sim = (sum([DataStructs.DiceSimilarity(i.mol1_fp, i.mol2_fp) for i in matched_reactants])/len(matched_reactants) + sum([DataStructs.DiceSimilarity(i.mol1_fp, i.mol2_fp) for i in matched_products])/len(matched_products))/2 
            rev_sim = (sum([DataStructs.DiceSimilarity(i.mol1_fp, i.mol2_fp) for i in matched_reactants_rev])/len(matched_reactants_rev) + sum([DataStructs.DiceSimilarity(i.mol1_fp, i.mol2_fp) for i in matched_products_rev])/len(matched_products_rev))/2

            # Choose most similar direction + calculate overall similarity
            reaction_level_sim = max([fwd_sim, rev_sim])
            sims[f'R{idx+1}/R{idx2+idx+2}'] = (0.75*reaction_level_sim) + (0.25*transformation_level_sim)

        else:
            sims[f'R{idx+1}/R{idx2+idx+2}']  = 0

sims = sorted(sims.items(), key = lambda a:a[1], reverse = True)
with open('similarity.csv','w') as f:
    w = csv.writer(f)
    w.writerows(sims)