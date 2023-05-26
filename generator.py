import networkx as nx
from networkx.classes.graph import Graph as NXGraph 

from dimod.binary_quadratic_model import BinaryQuadraticModel as BQM

from dwave.samplers.sa import SimulatedAnnealingSampler
SASampler = SimulatedAnnealingSampler()

from dimod.sampleset import SampleSet

from dimod.vartypes import SPIN

import random
from random import randint
            
SEED = 23451234
random.seed(SEED)

from networkx.generators.random_graphs import connected_watts_strogatz_graph as WSGraph
from networkx.generators.random_graphs import random_powerlaw_tree as PLTree
from networkx.generators.random_graphs import barabasi_albert_graph as BAGraph

def solve(graph: NXGraph) -> int :
    model = BQM(vartype=SPIN)
    print(graph)
    for n in graph.nodes.data():
        model.add_variable(n[0], n[1]['w'])
    for e in graph.edges.data():
        model.add_quadratic(e[0], e[1], e[2]['w'])
    result : SampleSet = SASampler.sample(bqm=model, num_reads=10000)
    minn = 0
    for r in result.record:
        if (minn > r.energy):
            minn = r.energy
    return minn

def output(graph: NXGraph, name: str):
    nx.drawing.nx_agraph.write_dot(graph, f'dotfiles_2/{name}.dot')
    sol = solve(graph)
    with open(f'instances_2/{name}.txt', 'w', encoding='utf-8') as f:
        f.write(f'{name}\n')
        f.write(f'{len(graph)}\n')
        for i in range(len(graph)):
            u = graph.nodes[i]['w']
            f.write(f'{u}\n')
        for e in graph.edges.data():
            f.write(f'{e[0]} {e[1]} {e[2]["w"]}\n')
        f.write('-1\n')
        f.write(f'{int(sol)}\n')

def add_weight(graph: NXGraph, lower: int, upper: int):
    for n in graph.nodes():
        u = 0
        while (True):
            u = randint(lower, upper)
            if not((u >= -1) and (u <= 1)):
                break
        graph.nodes[n]['w'] = u  

    for e in graph.edges():
        u = 0
        while (True):
            u = randint(lower, upper)
            if not((u >= -1) and (u <= 1)):
                break
        graph.edges[e]['w'] = u

nodes = [20,50,100]
weights = [10,30,60,100]


import math
for n in nodes:
    for w in weights:
        graph = WSGraph(n, math.floor(math.log2(n) * 2), 0.3, 1000, SEED)
        add_weight(graph, -w, w)
        output(graph,f'ws_{n}_{w}')
        
        graph = BAGraph(n, math.floor(math.log2(n) * 2), SEED)
        add_weight(graph, -w, w)
        output(graph,f'ba_{n}_{w}')
        
        graph = PLTree(n, seed=SEED, tries=1000)
        add_weight(graph, -w, w)
        output(graph,f'plt_{n}_{w}')

# graph = WSGraph(100, 12, 0.3, 100, SEED)

# add_weight(graph, -20, 20)

# output(graph, 'ws100')



# graph = BAGraph(100, 18, SEED)

# add_weight(graph, -20, 20)

# output(graph, 'ba100')

# nx.draw(graph)