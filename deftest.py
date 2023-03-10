import random
from dwave.cloud import Client
from dwave.system import DWaveSampler, EmbeddingComposite, LazyFixedEmbeddingComposite
from dimod.binary.binary_quadratic_model import BinaryQuadraticModel
import dwave.inspector
from dimod.vartypes import SPIN
import dwave_networkx as dnx
import networkx as nx
from minorminer import find_embedding
from dwave.embedding.chain_strength import uniform_torque_compensation as UTC
from dwave.embedding import embed_bqm, unembed_sampleset, EmbeddedStructure
from functools import partial
from copy import deepcopy
CSampler = DWaveSampler(solver={'topology__type': 'chimera'})
Sampler = CSampler
hardware_graph : list[ tuple[int, int] ] = Sampler.edgelist
nx_graph = nx.Graph(hardware_graph)

EPS = 1e-9
cached : dict = {}

def YanStrength( source : BinaryQuadraticModel, embedding : EmbeddedStructure, multiplier: float):
    global cached
    print('Default called')
    if (cached == {}):
        cs_array = UTC(source, embedding)
        cached = deepcopy(cs_array)
        for i in source.variables:
            cs_array[i] *= multiplier 
        return cs_array
    else:
        cs_array = deepcopy(cached)
        for i in source.variables:
            cs_array[i] *= multiplier
        return cs_array



cs_list = [1.4, 1.5, 1.7, 2.0, 2.5, 3.0, 4.0]
at_list = [10, 20, 50]

# cs_list = [10,25,50]
# at_list = [5,10,20]

s1_table = [[],[],[]]
s2_table = [[],[],[]]
avg_table = [[],[],[]]
ncb_table = [[],[],[]]

print('reading input')

n = int(input())

model = BinaryQuadraticModel(vartype=SPIN)

for i in range(n):
    model.set_linear(i,int(input()))

while(True):
    line = (input())
    listt = line.split(' ')
    # print(line)
    # print(listt)
    temp = list(map(int,listt))
    if (temp == [-1]):
        break
    # print(temp)
    u,v,a = temp
    model.set_quadratic(u,v,a)

optimal_solution = int(input())

# TEST

composite = LazyFixedEmbeddingComposite(Sampler)

for i in range(len(cs_list)):
    mult = cs_list[i]
    finalStrength = partial(YanStrength, multiplier = mult)


    for j in range(len(at_list)):
        at = at_list[j]
        print('Current run:')
        print(f'Chain strength mult: {cs_list[i]}')
        print(f'annealing time = {at}')
        print(f'optimal solution = {optimal_solution}')

        sample_set = composite.sample(model, num_reads=100, chain_strength = finalStrength, annealing_time = at)

        s1_cnt = 0
        s2_cnt = 0
        ncb_cnt = 0
        avg = 0
        # print(sample_set)
        for r in sample_set.record:
            # print(dir(r))
            # print(optimal_solution)                 
            # print(r.energy)                 
            if ( r.energy < optimal_solution * 0.995 + EPS):
                s1_cnt += r.num_occurrences
                                    
            if ( r.energy < optimal_solution * 0.97 + EPS):
                s2_cnt += r.num_occurrences
            avg += r.energy * r.num_occurrences
            ncb_cnt += r.chain_break_fraction
        print(f'Number of samples within .5% = {s1_cnt}')
        print(f'Number of optimal within 3% = {s2_cnt}')
        print(f'Average solution value = {avg/100}')
        print(f'Percentage of chains broken = {ncb_cnt}')


        avg_table[j].append(avg/100)
        s1_table[j].append(s1_cnt)
        s2_table[j].append(s2_cnt)
        ncb_table[j].append(ncb_cnt)

    # dwave.inspector.show(sample_set)

for r in s1_table :
    print(r)
print()
print()

for r in s2_table :
    print(r)
print()
print()

for r in avg_table :
    print(r)
print()
print()

for r in ncb_table :
    print(r)