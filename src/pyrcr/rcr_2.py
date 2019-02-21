# -*- coding: utf-8 -*-

"""Group 2 Implementation. Alfonso and Rana"""

import os
import tempfile
import unittest

from pybel import BELGraph
from pybel.dsl import gene, protein, rna
from pybel.manager import Manager
from pybel.testing.utils import n
from pybel_jupyter import to_jupyter
import networkx as nx
import pandas as pd
import csv

from pybel.constants import *
from pybel.dsl.node_classes import CentralDogma

relations = {
    INCREASES: 1,
    DIRECTLY_INCREASES: 1,
    DECREASES: -1,
    DIRECTLY_DECREASES: -1,
    RATE_LIMITING_STEP_OF: 0,
    CAUSES_NO_CHANGE: 0,
    REGULATES: 0,
    NEGATIVE_CORRELATION: -1,
    POSITIVE_CORRELATION: 1,
    ASSOCIATION: 0,
    HAS_MEMBER: 0,
    HAS_PRODUCT: 0,
    HAS_COMPONENT: 0,
    HAS_VARIANT: 0,
    HAS_REACTANT: 0,
    TRANSLATED_TO: 0,
    TRANSCRIBED_TO: 0,
    IS_A: 0,
    SUBPROCESS_OF: 0,
    ANALOGOUS_TO: 0,
    BIOMARKER_FOR: 0,
    PROGONSTIC_BIOMARKER_FOR: 0,
    EQUIVALENT_TO: 0
}


def path_validation(graph, nodes_dict, source_node, target_node):
    """
    this method takes a dict for experimental data for nodes and adds them as attributes in the graph
    then it gets the shortest path between two nodes and compares the experimental values (node attribute)
    with predicted values (edge relations). It returns True if they match and False if they don't
    """
    if source_node == target_node:
        return "No path between same nodes"
    if not nx.has_path(graph, source_node, target_node):
        return "path doesn't exist"
    nx.set_node_attributes(graph, nodes_dict)
    if source_node not in nodes_dict:
        return "Experimental data missing in path"
    if target_node not in nodes_dict:
        return "Experimental data missing in path"
    s_node_val=graph.node[source_node]['value']
    path= nx.shortest_path(graph, source_node, target_node)
    i=0
    nodes_val=1
    edges_val=1
    for n in path:
        while i<len(path)-1:
            if path[i] not in nodes_dict:
                return "Experimental data missing in path"
            nodes_val*=graph.node[path[i]]['value']
            edges_val*=edge_relation(graph, path[i], path[i+1])
            i+=1
    nodes_val*=graph.node[path[i]]['value']
    if edges_val==0:
        return "path cannot be evaluated"
    if nodes_val >0 and edges_val >0 or nodes_val <0 and edges_val <0:
        return True
    else :
        return False


def values_from_excel(filepath):
    """
    This method gets a csv file with entities and their experiemental values
    and returns the values into a dictionary
    """
    dict01={}
    df=pd.read_csv(filepath)
    for index, row in df.iterrows():
        dict01[row["Gene.symbol"]] = float(row["logFC"])
    return dict01

def values_to_nodes(exp_val, graph):
    """
    This method maps the values from a dictionary into attributes of nodes with the same key
    """
    dict02={}
    for key, val in exp_val.items():
        for i, data in graph.nodes(data=True):
            if not isinstance(CentralDogma, i):
                continue
            if i.name == key:
                dict02[i]={'value':val}
    return dict02

def edge_relation(graph, source_node, target_node):
    """
    This function gets the edge between two nodes and returns the relation as a integer
    according to the relations dictionary
    """
    for iden, edge_d in graph[source_node][target_node].items():
        return relations[edge_d['relation']]
