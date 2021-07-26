# find-ligand-for-CB7
MCTS implementation branched off ChemTS to find ligands for CB7's portal

## Set up the environment 
### Create a conda environment and activate it
```conda create --name test-chemts python=2.7```

```conda activate test-chemts```
### Install the rdkit
```conda install -c rdkit rdkit```
### Install openbabel 
```conda install -c openbabel openbabel```
### Install requirements
```python2 -m pip install -r requirements.txt```

## Run the MCTS, i is an integer
python2 mcts_ligand.py i

## Train the MCTS
python2 train_RNN.py
