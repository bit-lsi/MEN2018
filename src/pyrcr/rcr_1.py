# -*- coding: utf-8 -*-

"""Group 1 Implementation. Mohammed and Mahmudda"""
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 19 13:17:18 2019

@author: INSTHassaM3
"""
# -*- coding: utf-8 -*-

import os
import tempfile
import unittest
import numpy as np
from pybel import BELGraph
from pybel.dsl import gene, protein, rna
from pybel.manager import Manager
from pybel.testing.utils import n
import pandas as pd
import networkx as nx

HGNC = 'HGNC'

protein_a = protein(namespace=HGNC, name='a')
protein_b = protein(namespace=HGNC, name='b')
gene_c = gene(namespace=HGNC, name='c')
rna_d = rna(namespace=HGNC, name='d')
protein_e = protein(namespace=HGNC, name='e')
gene_f = gene(namespace=HGNC, name='f')
protein_g = protein(namespace=HGNC, name='g')
protein_h = protein(namespace=HGNC, name='h')
protein_i = protein(namespace=HGNC, name='i')
protein_j = protein(namespace=HGNC, name='j')


def __make_graph_1() -> BELGraph:
    graph = BELGraph(
        name='Lab course example',
        version='1.1.0',
        description='',
        authors='LSI',
        contact='lsi@uni-bonn.de',
    )

    graph.add_node_from_data(protein_a)
    graph.add_node_from_data(protein_b)
    graph.add_node_from_data(gene_c)
    graph.add_node_from_data(rna_d)

    graph.add_increases(
        protein_a,
        protein_b,
        citation='1',
        evidence='Evidence 1',
        annotations={'Annotation': 'foo'}
    )

    graph.add_increases(
        rna_d,
        protein_a,
        citation='2',
        evidence='Evidence 2',
        annotations={'Annotation': 'foo'}
    )

    graph.add_decreases(
        gene_c,
        protein_b,
        citation='3',
        evidence='Evidence 3',
        annotations={'Annotation': 'foo'}
    )

    return graph

def make_graph(data_file):
    graph_1 = __make_graph_1()
    vals = pd.read_csv(data_file,',')
    vals.columns = ['val','gene']
    vals = {k:v for k,v in zip(vals['gene'].values,vals['val'].values)}
    
    nodes = {}
    
    for i in graph_1.nodes():
        nodes[i] = vals[i.name]
    
    nx.set_node_attributes(graph_1,nodes,name='data')
    return graph_1
#chnge to dict
relations = {
    'INCREASES': +1,
    'DIRECTLY_INCREASES': +1,
    'DECREASES': -1,
    'DIRECTLY_DECREASES': -1,
    'RATE_LIMITING_STEP_OF': 0, 
    'CAUSES_NO_CHANGE': 0,
    'REGULATES': 0,
    'NEGATIVE_CORRELATION': -1,
    'POSITIVE_CORRELATION': +1,
    'ASSOCIATION': 0,
    'HAS_MEMBER': 0,
    'HAS_PRODUCT': 0,
    'HAS_COMPONENT': 0,
    'HAS_VARIANT': 0,
    'HAS_REACTANT': 0,
    'TRANSLATED_TO': 0,
    'TRANSCRIBED_TO': 0,
    'IS_A': 0,
    'SUBPROCESS_OF': 0,
    'ANALOGOUS_TO': 0,
    'BIOMARKER_FOR': 0,
    'PROGONSTIC_BIOMARKER_FOR': 0,
    'EQUIVALENT_TO': 0
}

def get_all(graph):
    for node_f in graph.nodes(data=True):
        nodef_val = node_f[1]['data']
        for node_t in graph.nodes(data=True):
            nodet_val = node_t[1]['data']
            if not nx.has_path(graph,node_f[0],node_t[0]):
                continue
            elif node_f == node_t:
                continue
            
            edg_val=1
            path = nx.shortest_path(graph,node_f[0],node_t[0])
            path_len=nx.shortest_path_length(graph,node_f[0],node_t[0])
            #print(path)
            for i in range(path_len):
                edge_data =list(graph.get_edge_data(path[i],path[i+1]).values())[0]['relation'].upper()
                edg_val *= relations[edge_data] if edge_data in relations else 0
            print(path)
            print(np.sign(edg_val*nodef_val)==np.sign(nodet_val))

def is_correct(graph,node_from,node_to):
    assert graph.has_node(node_from)
    assert graph.has_node(node_to)
    assert node_from != node_to
    assert nx.has_path(graph,node_from,node_to)
    
    # Get dara
    data = nx.get_node_attributes(graph,'data')
    nodef_val = data[node_from]
    nodet_val = data[node_to]
    
    edg_val=1
    path = nx.shortest_path(graph,node_from,node_to)
    path_len=nx.shortest_path_length(graph,node_from,node_to)
    for i in range(path_len):
        edge_data = list(graph.get_edge_data(path[i],path[i+1]).values())[0]['relation'].upper()
        edg_val *= relations[edge_data] if edge_data in relations else 0
    return np.sign(edg_val*nodef_val)==np.sign(nodet_val)


graph_1 = make_graph('data/test.csv')
print(is_correct(graph_1,protein_a,protein_b))
print(get_all(graph_1))