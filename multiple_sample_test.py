import random
from dwave.cloud import Client

from dwave.system import DWaveSampler, EmbeddingComposite

from dimod.binary.binary_quadratic_model import BinaryQuadraticModel

import dwave.inspector

from dimod.vartypes import SPIN

import dwave_networkx as dnx

# D-cstrWave 2X
C = dnx.chimera_graph(12, 12, 4)

CSampler = DWaveSampler(solver={'topology__type': 'chimera'})

# print(CSampler.properties)

# print('input chain strength:')
# cstr = int(input())
# print('input anneal time:')
# anneal_time = float(input())

cstr_list = [10000,1000000,100000000,10000000000,1000000000000]
at_list = [5,10,20,50,100,200,500]
edge_cache = []
graph_nodes = {}

for i in range(0,512):
    if (i in CSampler.nodelist):
        graph_nodes[i]=0

for i in range(509):
    a = int(input())
    b = int(input())
    edge_cache.append((a,b))

for cstr in cstr_list:
    graph_edges = {}

    for i in range(509):
        if (i == 0):
            graph_edges[edge_cache[i]] = -1
        else:
            graph_edges[edge_cache[i]] = -cstr

    model = BinaryQuadraticModel(
        graph_nodes,
        graph_edges,
        SPIN
    )

    for at in at_list:
    # print(model)
        print('Current params:')
        print(f'chain_strength = {cstr}')
        print(f'anneal_time = {at}')
        print('results:')

        sample_result = CSampler.sample(model, num_reads=100, annealing_time = at)
        # print(sample_result.record.dtype.names)
        for r in sample_result.record:
            if (r.energy == -508*cstr-1):
                print(f'Correct samples: {r.num_occurrences}')
            if (r.energy == -508*cstr+1):
                print(f'Near samples: {r.num_occurrences}')
        # print(sample_result.info)
        # dwave.inspector.show(sample_result)