# -*- coding: utf-8 -*-

"""Group 1 Implementation. Mohammed and Mahmudda"""
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 19 13:17:18 2019

@author: INSTHassaM3
"""
# -*- coding: utf-8 -*-

import numpy as np
import networkx as nx
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