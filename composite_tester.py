import random
import dimod
import dwave.inspector
import csv
import math
import dwave_networkx as dnx
import networkx as nx
import time
from dwave.cloud import Client
from dwave.system import DWaveSampler, EmbeddingComposite, LazyFixedEmbeddingComposite
from dimod.binary.binary_quadratic_model import BinaryQuadraticModel
from dimod.vartypes import SPIN
from dimod import SampleSet
from minorminer import find_embedding
from dwave.embedding.chain_strength import uniform_torque_compensation as UTC
from dwave.embedding import embed_bqm, unembed_sampleset, EmbeddedStructure
from functools import partial
from copy import deepcopy
from dwave.samplers import SimulatedAnnealingSampler


SASampler = SimulatedAnnealingSampler()
# CSampler = DWaveSampler(solver={'topology__type': 'pegasus'})

from helpers import HillClimbChimeraSampler


def fe(S, T, **kwargs):
    print('find_embedding called')
    return find_embedding(S, T, random_seed=123123)

FCSampler = HillClimbChimeraSampler()

# C_Composite = LazyFixedEmbeddingComposite(CSampler, find_embedding=fe)
FC_Composite = LazyFixedEmbeddingComposite(FCSampler, find_embedding=fe)





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
cached : dict = {}


def compute(j: float, composite: dimod.Sampler, runs: int) -> dict[str,float]:
    print('Final run:')
    print(f'strength: {j}')
    
    print(time.time())
    sample_set = composite.sample(model, num_reads = runs, chain_strength = j, seed = 123)
    print(composite.embedding)
    print(time.time())

    s1_cnt: float = 0
    s2_cnt: float = 0
    ncb_cnt: float = 0
    avg: float = 0
    best: float = 0
    for r in sample_set.record:           
        if ( r.energy < optimal_solution * 0.995 + EPS):
            s1_cnt += r.num_occurrences
                                
        if ( r.energy < optimal_solution * 0.97 + EPS):
            s2_cnt += r.num_occurrences
        avg += r.energy * r.num_occurrences
        best = min(best, r.energy)
        ncb_cnt += r.chain_break_fraction
    
    res: dict[str,float] = {}
    res['p5no'] = s1_cnt / runs
    res['3no'] = s2_cnt / runs
    res['avg'] = avg / runs
    res['break'] = ncb_cnt / runs
    res['best'] = best
    print(res)
    return res


start_time = time.time()

compute(100000, FC_Composite, 1)

print(f'Run time: {(time.time() - start_time)} seconds')
