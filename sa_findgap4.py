import random
import dimod
import dwave.inspector
import csv
import math
import numpy
import dwave_networkx as dnx
import networkx as nx
import time
from dwave.cloud import Client
from dwave.system import DWaveSampler, EmbeddingComposite, LazyFixedEmbeddingComposite
from dimod.binary.binary_quadratic_model import BinaryQuadraticModel
from dimod.vartypes import SPIN
from minorminer import find_embedding
from dwave.embedding.chain_strength import uniform_torque_compensation as UTC
from dwave.embedding import embed_bqm, unembed_sampleset, EmbeddedStructure
from functools import partial
from copy import deepcopy

from dwave.samplers import SimulatedAnnealingSampler
SASampler = SimulatedAnnealingSampler()
CSampler = DWaveSampler(solver={'topology__type': 'pegasus'})

hw = dnx.chimera_graph(16,16)




def fe(S, T, **kwargs):
    return find_embedding(S, T, random_seed=123123)

C_Composite = LazyFixedEmbeddingComposite(CSampler, find_embedding=fe)


print('reading input')

iname: str = input()

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


print(f'Optimal solution: {optimal_solution}')

EPS = 1e-9


cache : dict[float,dict[any,float]] = {}

embedding : dict = find_embedding(model.quadratic.keys(),hw.edges.keys())
nx_graph = nx.Graph(hw.edges.keys())

def compute(j: float, runs: int) -> dict[str,float]:
    # print('Final run:')
    # print(f'strength: {j}')

    embedded_model = embed_bqm(model, embedding, nx_graph, j)
    sample_set = SASampler.sample(embedded_model, num_reads = runs)

    s1_cnt: float = 0
    s2_cnt: float = 0
    ncb_cnt: float = 0
    avg: float = 0
    best: float = 0
    best_sample = 0

    break_frac : dict = {}
    
    for variable in embedding.keys():
        break_frac[variable] = 0

    for r in sample_set.record:           
        if ( r.energy < optimal_solution * 0.995 + EPS):
            s1_cnt += r.num_occurrences
                                
        if ( r.energy < optimal_solution * 0.97 + EPS):
            s2_cnt += r.num_occurrences
        avg += r.energy * r.num_occurrences
        if (r.energy < best):
            best = r.energy
            best_sample = r
        for variable in embedding.keys():
            sus = 0
            mapping = embedding[variable]
            for i in range(len(embedding[variable])-1):
                u = sample_set.variables.index(mapping[i])
                v = sample_set.variables.index(mapping[i+1])
                if not (r.sample[u] == r.sample[v]):
                    sus = 1
            if (sus):
                break_frac[variable] += 1


    res: dict[str,float] = {}
    res['p5no'] = s1_cnt / runs
    res['3no'] = s2_cnt / runs
    res['avg'] = avg / runs
    res['best'] = best
    for v in break_frac.keys():
        res[v] = break_frac[v] / runs
    # print(res)
    return res

def compute_dwave(j: float, composite: dimod.Sampler, runs: int) -> dict[str,float]:
    # print('Final run:')
    # print(f'strength: {j}')

    sample_set = composite.sample(model, num_reads = runs, chain_strength = j)

    s1_cnt: float = 0
    s2_cnt: float = 0
    ncb_cnt: float = 0
    avg: float = 0
    best: float = 0
    best_sample = 0

    for r in sample_set.record:           
        if ( r.energy < optimal_solution * 0.995 + EPS):
            s1_cnt += r.num_occurrences
                                
        if ( r.energy < optimal_solution * 0.97 + EPS):
            s2_cnt += r.num_occurrences
        avg += r.energy * r.num_occurrences
        if (r.energy < best):
            best = r.energy
            best_sample = r
        ncb_cnt += r.chain_break_fraction

    res: dict[str,float] = {}
    res['p5no'] = s1_cnt / runs
    res['3no'] = s2_cnt / runs
    res['avg'] = avg / runs
    res['break'] = ncb_cnt / runs
    res['best'] = best
    # print(res)
    return res


DEF_BS_RUNS: int = 15


DEF_BBOUND : float = 0.1

def get_cs(bound: float) -> float:
    start_low = 0.01
    start_high: float = 2 * UTC(model)
    l: float = 0
    r: float = 0
    
    l = start_low
    r = start_high
    for i in range(DEF_BS_RUNS):
        mid = (l+r)/2
        sus = 0
        compute_res = compute(mid, 1000)
        print(mid)
        for v in model.variables:
            print(f'{v} {compute_res[v]}')
            if (compute_res[v] > bound + 1e-9):
                sus = 1
        if (sus):
            l = mid
        else:
            r = mid
    res = (l + r) / 2
    return res

start_time = time.time()

result : list[float] = []

for i in range(10):
    result.append(get_cs((i+1)*(0.02)))


f = open(f'findgap4_runtime.csv', mode='a')
print(f'Run time: {(time.time() - start_time)} seconds')
f.write(f'{iname},{(time.time() - start_time)}\n')
f.close()

f = open(f'findgap4_results/{iname}.csv', 'w', encoding='utf-8')
f.write('abs,rel,p5no,3no,best,avg,break\n')
sus: float = UTC(model)
def output(j):
    print(j)
    print(j/sus)
    u = compute_dwave(j, C_Composite, 1000)
    print(u)
    f.write(f'{j},{j/sus},{u["p5no"]},{u["3no"]},{u["best"]},{u["avg"]},{u["break"]}\n')

for cs in result:
    output(cs)
output(sus)

f.close()
print(UTC(model))   