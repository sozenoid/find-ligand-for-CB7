# find-ligand-for-CB7
MCTS implementation branched off ChemTS (https://github.com/tsudalab/ChemTS) to find ligands for CB7's portal. The smiles dataset is located in ```data``` while the pretrained model and weights for smiles with at most 12 heavy atoms is located in RNN-model.

## Set up the environment 
### Create a conda environment and activate it
```conda create --name ligand-for-CB7 python=2.7```

```conda activate ligand-for-CB7```
### Install the rdkit
```conda install -c rdkit rdkit```
### Install openbabel 
```conda install -c openbabel openbabel```
### Install requirements
```python2 -m pip install -r requirements.txt```

### Edit the paths in xtb_binding.py
```wdir = "PATH_TO/cap_design_for_CB7_XTB"```

```obabel_path = "PATH_TO/obabel"```

```ledock_path = "" # Ledock should be in the same folder as wdir, leave blank```

```xtb_path = "PATH_TO/xtb-6.3.3/bin/xtb"```


## Run the MCTS, i is an integer
```cd cap_design_for_CB7_XTB```

```python2 mcts_ligand.py i```

The results will be displayed to STDOUT then stored and compressed in the 'outputs' folder. 

## Train the MCTS
Edit the train_RNN.py file to set up the maximum smiles length (maxlen=31 by default). Then start training. See the ChemTS for training issues and how to train using GPUs. 

```cd train_RNN```

```python2 train_RNN.py```
