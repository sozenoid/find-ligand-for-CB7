# find-ligand-for-CB7
MCTS implementation branched off ChemTS (https://github.com/tsudalab/ChemTS) to find ligands for CB7's portal

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
```cd cap_design_for_CB7_XTB```

```python2 mcts_ligand.py i```

## Train the MCTS
```cd train_RNN```

```python2 train_RNN.py```
