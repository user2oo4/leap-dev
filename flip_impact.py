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

from dwave.samplers import SimulatedAnnealingSampler
SASampler = SimulatedAnnealingSampler()

CSampler = DWaveSampler(solver={'topology__type': 'chimera'})
Sampler = CSampler
hardware_graph : list[ tuple[int, int] ] = Sampler.edgelist
nx_graph = nx.Graph(hardware_graph)

EPS = 1e-9


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

# PREPROCESS

flip_mask : dict[int,int] = {}
for i in model.variables:
    flip_mask[i] = 1
    
def flip(model, mask):
    for (u,v) in model.quadratic.keys():
        w = model.get_quadratic(u,v)
        w *= mask[u]
        w *= mask[v]
        model.set_quadratic(u,v,w)

    for u in model.variables:
        w = model.get_linear(u)
        w *= mask[u]
        model.set_linear(u,w)

# TEST
print('Testing')
## this is the stupidest fucking workaround ever but it has to be done
## since passing random_seed from LFEC class doesn't work for some BS reason

def fe(S, T, **kwargs):
    return find_embedding(S, T, random_seed=123123)

f = open(f'results/flip_impact.csv', 'w', encoding='utf-8')
f.write('node,sum,delta\n')
composite = LazyFixedEmbeddingComposite(Sampler, find_embedding=fe)

neutral_avg = 0

for i in range(0,1):
    print(f'Current run:{i-1}')
    print(f'optimal solution = {optimal_solution}')

    sample_set = composite.sample(model, num_reads=100, annealing_time = 50)

    s1_cnt = 0
    s2_cnt = 0
    ncb_cnt = 0
    avg = 0
    # print(sample_set)
    for r in sample_set.record:
        avg += r.energy * r.num_occurrences
    avg /= 100
    print(f'Average solution value = {avg}')
    neutral_avg = avg

sum_w = {}
for u in model.variables:    
    sum_w[u] = 0

for (u,v) in model.quadratic.keys():
    w = model.get_quadratic(u,v)
    sum_w[u] += w
    sum_w[v] += w

for u in model.variables:
    w = model.get_linear(u)    
    sum_w[u] += w



for i in model.variables:
    print(f'Current run:{i}')
    print(f'optimal solution = {optimal_solution}')

    sample_set = composite.sample(model, num_reads=100, annealing_time = 50)

    s1_cnt = 0
    s2_cnt = 0
    ncb_cnt = 0
    avg = 0
    # print(sample_set)
    for r in sample_set.record:
        avg += r.energy * r.num_occurrences
    avg /= 100
    print(f'Average solution value = {avg}')
    f.write(f'{i},{sum_w[i]},{avg-neutral_avg}\n')
    # break
    # dwave.inspector.show(sample_set)
f.close()