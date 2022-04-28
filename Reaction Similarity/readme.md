# Description

Calculates reaction similarity scores for all reaction pairs from a SMILES/SMARTS reaction file, combining reaction-level and transformation-level similarity scores.

## Requirements

This program requires Python >=3.8, and the RDKit package. RDKit can be installed using conda:

$	conda create -c conda-forge -n my-rdkit-env

$	conda activate my-rdkit-env 

## Usage

Run ‘reaction_similarity_aam.py’ to get overall similarity results using atom-to-atom mapping to calculate the transformation-level similarity.

Run ‘reaction_similarity_subtraction.py’ to get overall similarity results using fingerprint subtraction to calculate the transformation-level similarity.

#### Positional arguments

**fingerprint**

Choose a fingerprint type. Type ‘m’ to use the Morgan fingerprint, ‘ap’ to use the Atom Pair fingerprint, or ‘t’ to use the topological fingerprint. 

#### Options

**-file**

Choose an alternative reaction SMARTS file in csv format. By default, ‘SMARTS.csv’ will be used.

## Examples

‘reaction_similarity_subtraction.py m’ 

Calculate the reaction similarity of the reactions in ‘SMARTS.csv’ using the Morgan fingerprint and the fingerprint subtraction method.
