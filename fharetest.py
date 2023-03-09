from collections import OrderedDict
import random
from dwave.cloud import Client
from dwave.system import DWaveSampler, EmbeddingComposite, LazyFixedEmbeddingComposite
from dimod.binary.binary_quadratic_model import BinaryQuadraticModel
import dwave.inspector
from dimod.vartypes import SPIN
import dwave_networkx as dnx
import networkx as nx
from minorminer import find_embedding
from dwave.embedding import embed_bqm, unembed_sampleset, EmbeddedStructure
from functools import partial
from copy import deepcopy
CSampler = DWaveSampler(solver={'topology__type': 'chimera'})
Sampler = CSampler
hardware_graph : list[ tuple[int, int] ] = Sampler.edgelist
nx_graph = nx.Graph(hardware_graph)


from queue import PriorityQueue
# from typing import PriorityQueue


EPS = 1e-9
cached : dict = {}

def FastHareStrength( source : BinaryQuadraticModel, embedding : EmbeddedStructure, multiplier: float):
    global cached
    print('fhare called')
    if (cached == {}):
        
        embedded = embedding.embed_bqm(source, 1000000)
        to_logical: dict[int,int] = {}
        cs_array : dict[int,int] = {}

        j_field : dict[int,int] = {}
        h_field : dict[int,int] = {}

        gt : dict[int,set[int]] = {}

        leaves : PriorityQueue[tuple[int,int]] = PriorityQueue()

        for i in source.variables:
            cs_array[i] = 0
            for j in embedding[i]:  
                to_logical[j] = i
                j_field[j] = 0
                h_field[j] = abs(embedded.linear[j])
                gt[j] = set()

        for u in embedded.quadratic.keys():
            if (embedded.quadratic[u] != -1000000):
                # print(u)
                j_field[u[0]] += abs(embedded.quadratic[u])
                j_field[u[1]] += abs(embedded.quadratic[u])
            else:
                gt[u[0]].add(u[1])
                gt[u[1]].add(u[0])
        
        for i in source.variables:
            print(f'processing logical qubit {i}')
            leaves = PriorityQueue()
            for j in embedding[i]:
                if (len(gt[j]) == 1):
                    leaves.put(item=(j_field[j],j))
            for cnt in range(len(embedding[i])-1):
                if (leaves.empty()):
                    print(f'error: queue empty')
                w, u = leaves.get()
                v = list(gt[u])[0]
                print(f'edge {u} {v} {w}')
                gt[u].remove(v)
                gt[v].remove(u)

                if (cs_array[i] < w):
                    cs_array[i] = w
                
                j_field[v] += abs(w)
                if (len(gt[v]) == 1):
                    leaves.put(item=(j_field[v],v))

        cached = deepcopy(cs_array)
        for i in source.variables:
            cs_array[i] *= multiplier 
        return cs_array
    else:
        cs_array = deepcopy(cached)
        for i in source.variables:
            cs_array[i] *= multiplier
        return cs_array



cs_list = [0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8]
at_list = [10,20,50]

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
    finalStrength = partial(FastHareStrength, multiplier = mult)


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