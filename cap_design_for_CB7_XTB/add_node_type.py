from subprocess import Popen, PIPE
from math import *
import random
import numpy as np
from copy import deepcopy
from types import IntType, ListType, TupleType, StringTypes
import itertools
import time
import math
import argparse
import subprocess
from load_model import loaded_model
from keras.preprocessing import sequence
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import Descriptors
import sys
from rdkit.Chem import AllChem
from rdkit.Chem import MolFromSmiles, MolToSmiles
from make_smile import zinc_data_with_bracket_original, zinc_processed_with_bracket
from rdock_test import rdock_score
import sascorer
import pickle
import gzip
import networkx as nx
from rdkit.Chem import rdmolops


maxlen=31

def expanded_node(model,state,val):

    all_nodes=[]

    end="\n"

    #position=[]
    position=[]
    position.extend(state)
    total_generated=[]
    new_compound=[]
    get_int_old=[]
    for j in range(len(position)):
        get_int_old.append(val.index(position[j]))

    get_int=get_int_old

    x=np.reshape(get_int,(1,len(get_int)))
    x_pad= sequence.pad_sequences(x, maxlen=maxlen, dtype='int32',
        padding='post', truncating='pre', value=0.)

    for i in range(30):
        predictions=model.predict(x_pad)
        #print "shape of RNN",predictions.shape
        preds=np.asarray(predictions[0][len(get_int)-1]).astype('float64')
        preds = np.log(preds) / 1.0
        preds = np.exp(preds) / np.sum(np.exp(preds))
        next_probas = np.random.multinomial(1, preds, 1)
        next_int=np.argmax(next_probas)
        #get_int.append(next_int)
        all_nodes.append(next_int)

    all_nodes=list(set(all_nodes))

    print all_nodes





#total_generated.append(get_int)
#all_posible.extend(total_generated)






    return all_nodes


def node_to_add(all_nodes,val):
    added_nodes=[]
    for i in range(len(all_nodes)):
        added_nodes.append(val[all_nodes[i]])

    print added_nodes

    return added_nodes



def chem_kn_simulation(model,state,val,added_nodes):
    all_posible=[]

    end="\n"
    #val2=['C', '(', ')', 'c', '1', '2', 'o', '=', 'O', 'N', '3', 'F', '[C@@H]', 'n', '-', '#', 'S', 'Cl', '[O-]', '[C@H]', '[NH+]', '[C@]', 's', 'Br', '/', '[nH]', '[NH3+]', '4', '[NH2+]', '[C@@]', '[N+]', '[nH+]', '\\', '[S@]', '5', '[N-]', '[n+]', '[S@@]', '[S-]', '6', '7', 'I', '[n-]', 'P', '[OH+]', '[NH-]', '[P@@H]', '[P@@]', '[PH2]', '[P@]', '[P+]', '[S+]', '[o+]', '[CH2-]', '[CH-]', '[SH+]', '[O+]', '[s+]', '[PH+]', '[PH]', '8', '[S@@+]']
    for i in range(len(added_nodes)):
        #position=[]
        position=[]
        position.extend(state)
        position.append(added_nodes[i])
        #print state
        #print position
        #print len(val2)
        total_generated=[]
        new_compound=[]
        get_int_old=[]
        for j in range(len(position)):
            get_int_old.append(val.index(position[j]))

        get_int=get_int_old

        x=np.reshape(get_int,(1,len(get_int)))
        x_pad= sequence.pad_sequences(x, maxlen=maxlen, dtype='int32',
            padding='post', truncating='pre', value=0.)
        while not get_int[-1] == val.index(end):
            predictions=model.predict(x_pad)
            #print "shape of RNN",predictions.shape
            preds=np.asarray(predictions[0][len(get_int)-1]).astype('float64')
            preds = np.log(preds) / 1.0
            preds = np.exp(preds) / np.sum(np.exp(preds))
            next_probas = np.random.multinomial(1, preds, 1)
            #print predictions[0][len(get_int)-1]
            #print "next probas",next_probas
            #next_int=np.argmax(predictions[0][len(get_int)-1])
            next_int=np.argmax(next_probas)
            a=predictions[0][len(get_int)-1]
            next_int_test=sorted(range(len(a)), key=lambda i: a[i])[-10:]
            get_int.append(next_int)
            x=np.reshape(get_int,(1,len(get_int)))
            x_pad = sequence.pad_sequences(x, maxlen=maxlen, dtype='int32',
                padding='post', truncating='pre', value=0.)
            if len(get_int)>maxlen:
                break
        total_generated.append(get_int)
        all_posible.extend(total_generated)


    return all_posible



def predict_smile(all_posible,val):


    new_compound=[]
    for i in range(len(all_posible)):
        total_generated=all_posible[i]

        generate_smile=[]

        for j in range(len(total_generated)-1):
            generate_smile.append(val[total_generated[j]])
        generate_smile.remove("&")
        new_compound.append(generate_smile)

    return new_compound


def make_input_smile(generate_smile):
    new_compound=[]
    for i in range(len(generate_smile)):
        middle=[]
        for j in range(len(generate_smile[i])):
            middle.append(generate_smile[i][j])
        com=''.join(middle)
        new_compound.append(com)
    #print new_compound
    #print len(new_compound)

    return new_compound



def check_node_type(new_compound, back_name):
    node_index=[]
    valid_compound=[]
    all_smile=[]
    distance=[]

    score=[]
    for i in range(len(new_compound)):
        try:
            ko = Chem.MolFromSmiles(new_compound[i])
            
        except:
            ko=None
        
        if ko!=None:
            try:
                molscore=MolFromSmiles(new_compound[i])
            except:
                molscore=None
            if molscore!=None:
                SA_score = -sascorer.calculateScore(molscore)
            else:
                SA_score=1000
            cycle_list = nx.cycle_basis(nx.Graph(rdmolops.GetAdjacencyMatrix(MolFromSmiles(new_compound[i]))))
            if len(cycle_list) == 0:
                 cycle_length = 0
            else:
                cycle_length = max([ len(j) for j in cycle_list ])
            if cycle_length <= 6:
                cycle_length = 0
            if cycle_length==0:
                m, pc_covered ,rdid=rdock_score(new_compound[i])
                if m<10**10:
                    node_index.append(i)
                    valid_compound.append(new_compound[i])
                    #score.append(m+SA_score/5.0)
                    score.append(m*pc_covered**2)
                    """___________________SAVE THE MOL RESULTS_______________________________"""
                    with open("outputs/BACKUP-{}".format(back_name), "a") as a:
                        a.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(rdid, str(new_compound[i]), m, SA_score, pc_covered, m*pc_covered))


    return node_index,score,valid_compound
if __name__=="__main__":
    import os
    check_node_type(["c1ccccc1","c1ccccc1"] , "hey")
