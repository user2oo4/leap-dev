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
    nx.drawing.nx_agraph.write_dot(graph, f'dotfiles/{name}.dot')
    sol = solve(graph)
    with open(f'instances/{name}.txt', 'w', encoding='utf-8') as f:
        f.write(f'{len(graph)}\n')
        for i in range(len(graph)):
            u = graph.nodes[i]['w']
            f.write(f'{u}\n')
        for e in graph.edges.data():
            f.write(f'{e[0]} {e[1]} {e[2]["w"]}\n')
        f.write('-1\n')
        f.write(f'{sol}\n')

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


graph = WSGraph(50, 8, 0.2, 100, SEED)

add_weight(graph, -10, 10)

output(graph, 'ws50')


graph = WSGraph(20, 4, 0.2, 100, SEED)

add_weight(graph, -10, 10)

output(graph, 'ws20')

# nx.draw(graph)
    

