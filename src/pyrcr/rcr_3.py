# -*- coding: utf-8 -*-

"""Group 3 Implementation. Farid and Linling"""

import os
import tempfile
import unittest
from pybel_jupyter import to_jupyter
from pybel import BELGraph
from pybel.dsl import gene, protein, rna
from pybel.manager import Manager
from pybel.testing.utils import n
import pandas as pd
import itertools
import networkx as nx

HGNC = 'HGNC'

#dir_path = os.path.dirname(os.path.realpath(__file__))


class ManagerMixin(unittest.TestCase):
    def setUp(self):
        super(ManagerMixin, self).setUp()

        self.db_fd, self.db_file = tempfile.mkstemp()

        self.connection = 'sqlite:///' + self.db_file
        self.manager = Manager(connection=self.connection)

    def tearDown(self):
        os.close(self.db_fd)
        os.unlink(self.db_file)


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


def make_graph_1() -> BELGraph:
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


def make_graph_2() -> BELGraph:
    """Make an example graph."""
    graph = BELGraph(
                     name='Lab course example',
                     version='1.1.0',
                     description='',
                     authors='LSI',
                     contact='lsi@uni-bonn.de',
                     )

    graph.add_node_from_data(gene_f)
    graph.add_node_from_data(protein_e)
    graph.add_node_from_data(protein_b)

    graph.add_increases(
                     protein_e,
                     protein_b,
                     citation='1',
                     evidence='Evidence 1',
                     annotations={'Annotation': 'foo'},
                     )

    graph.add_increases(
                     gene_f,
                     protein_e,
                     citation='2',
                     evidence='Evidence 2',
                     annotations={'Annotation': 'foo2'}
                     )

    return graph


def make_graph_3() -> BELGraph:
    """Make an example graph.
        A -> B -| C
        D -| F -> C
        C -| F
        C -- G
        """
    graph = BELGraph(
                     name='Lab course example',
                     version='1.1.0',
                     description='',
                     authors='LSI',
                     contact='lsi@uni-bonn.de',
                     )

    graph.add_increases(protein_a, protein_b, n(), n())
    graph.add_decreases(protein_b, gene_c, n(), n())
    graph.add_decreases(rna_d, gene_f, n(), n())
    graph.add_increases(protein_e, gene_f, n(), n())
    graph.add_increases(gene_f, gene_c, n(), n())
    graph.add_association(gene_c, protein_g, n(), n())

    return graph


def make_graph_4() -> BELGraph:
    """Make an example graph.
        A -> B
        B -| C
        B -| D
        B -| E
        B -| F
        B -> G
        B -> H
        B -| H
        B -> I
        B -- J
        """
    graph = BELGraph(
                     name='Lab course example',
                     version='1.1.0',
                     description='',
                     authors='LSI',
                     contact='lsi@uni-bonn.de',
                     )

    graph.add_increases(protein_a, protein_b, n(), n())
    graph.add_decreases(protein_b, gene_c, n(), n())
    graph.add_decreases(protein_b, rna_d, n(), n())
    graph.add_decreases(protein_b, protein_e, n(), n())
    graph.add_decreases(protein_b, gene_f, n(), n())
    graph.add_increases(protein_b, protein_g, n(), n())
    graph.add_decreases(protein_b, protein_h, n(), n())
    graph.add_increases(protein_b, protein_h, n(), n())
    graph.add_increases(protein_b, protein_i, n(), n())
    graph.add_association(protein_b, protein_j, n(), n())

    return graph


example_1 = make_graph_1()
example_2 = make_graph_2()
example_3 = make_graph_3()
example_4 = make_graph_4()



def mapvalue(file, graph):
    df = pd.read_csv(file)
    gene_expression_dict = dict()
    key = list(df['Gene.symbol'])
    value = list(df['logFC'])
    for i in range(len(key)):
        gene_expression_dict[key[i]] = value[i]

    attrs2 = dict()
    for i in graph.nodes():
        if i.name in gene_expression_dict:
            attrs2[i] = {'value':gene_expression_dict[i.name]}

nx.set_node_attributes(graph, attrs2)
    graph.nodes(data=True)
    return graph


def graph_test(graph):
    graph = mapvalue('test.csv', graph)
    listnodes = list(graph.nodes())
    allpairs = itertools.permutations(listnodes,2)
    edgevaluedict =dict()
    edgevalues = nx.get_edge_attributes(graph,'relation')
    mapdict={'increases':1,'decreases':-1,'association':0}

    for n,e  in edgevalues.items():
        nodes=n[:2]
        edgevalue=mapdict[e]

        edgevaluedict[nodes]= edgevalue

    resultdict = dict()
    for n in listnodes:
        resultdict[n] = []

    for i in allpairs:
        if not nx.has_path(graph,i[0],i[1]):
            continue

        path = nx.shortest_path(graph,i[0],i[1])
        edgelist =[]
        firstvalue = graph.node[path[0]]['value']
        lastvalue = graph.node[path[-1]]['value']
        edgelist.append(firstvalue)
        edgelist.append(lastvalue)
        for m in range(1,len(path)):
            edgedata = edgevaluedict[(path[m-1],path[m])]
            edgelist.append(edgedata)
        checkvalue=1
        for r in edgelist:
            checkvalue = checkvalue*r
        if checkvalue < 0:

            resultdict[i[0]].append((path,'false'))
        elif checkvalue > 0:
            resultdict[i[0]].append((path,'true'))

    return resultdict


def search_node (graph, interesting_node, file_path):
    if interesting_node not in graph :
        print ('The node is not in the graph, please change anoher node')
    else:
        wholedict = graph_test(graph, file_path)
        return wholedict[interesting_node]


