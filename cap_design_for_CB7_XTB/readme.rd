conda create --name test-chemts python=2.7
conda activate test-chemts
conda install -c rdkit rdkit
conda install -c openbabel openbabel
python2 -m pip install -r requirements.txt
python2 mcts_ligand.py 1
python2 train_RNN.py
