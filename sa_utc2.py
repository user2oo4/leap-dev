import random
import dimod
import dwave.inspector
import csv
import math
import dwave_networkx as dnx
import networkx as nx
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
CSampler = DWaveSampler(solver={'topology__type': 'chimera'})

hw = dnx.chimera_graph(16,16)

class FakeChimeraSampler(dimod.Sampler, dimod.Structured):
    @property
    def properties(self) -> dict[str, any]:
        return SASampler.properties
    
    
    @property
    def parameters(self) -> dict[str, any]:
        return SASampler.parameters
    
    @property
    def nodelist(self):
        return hw.nodes.keys()
    
    
    @property
    def edgelist(self):
        return hw.edges.keys()
    
    def sample(self, bqm, **parameters):
        return SASampler.sample(bqm, **parameters)

FCSampler = FakeChimeraSampler()
Sampler = FCSampler

EPS = 1e-9
cached : dict = {}

def UTCI(bqm: BinaryQuadraticModel, embedding : EmbeddedStructure = None):
    deg: dict[any, int] = {}
    sum2: dict[any, float] = {}
    strength: dict[any, float] = {}
    for u in bqm.variables:
        deg[u] = 0
        sum2[u] = 0
    for (u,v) in bqm.quadratic.keys():
        w = bqm.quadratic[(u,v)]
        deg[u] += 1
        deg[v] += 1
        sum2[u] += w*w
        sum2[v] += w*w
    for u in bqm.variables:
        strength[u] = math.sqrt(sum2[u])
    return strength
    


def UTCIStrength( source : BinaryQuadraticModel, embedding : EmbeddedStructure, multiplier: dict[any, float]):
    global cached
    print('Default called')
    if (cached == {}):
        cs_array = UTCI(source, embedding)

        cached = deepcopy(cs_array)
        print(cached)
        for i in source.variables:
            cs_array[i] *= multiplier[i] 
        return cs_array
    else:
        cs_array = deepcopy(cached)
        for i in source.variables:
            cs_array[i] *= multiplier[i]
        return cs_array

s1_table = [[]]
s2_table = [[]]
avg_table = [[]]
ncb_table = [[]]

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

print(f'Optimal solution: {optimal_solution}')
# TEST

## this is the stupidest fucking workaround ever but it has to be done
## since passing random_seed from LFEC class doesn't work for some BS reason

def fe(S, T, **kwargs):
    return find_embedding(S, T, random_seed=123123)


composite = LazyFixedEmbeddingComposite(Sampler, find_embedding=fe)

current_mult : dict[any, float] = {}

for i in model.variables:
    current_mult[i] = 1

f = open(f'results/utci_scale_ba20.csv', 'w', encoding='utf-8')
f.write('node,scale,point5no,3no,avgsol,chainbreak\n')

i = 1.5
while(True):
    
    print('Current run:')

    file = open('results/individual_scale_ba20_dwave2.csv', newline='')
    reader = csv.reader(file)
    for row in reader:
        if (row[0] != 'node'):
            current_mult[int(row[0])] = 1 * i
    print(current_mult)
    finalStrength = partial(UTCIStrength, multiplier = current_mult)
    sample_set = composite.sample(model, num_reads=1000, chain_strength = finalStrength, annealing_time = 50)

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
    print(f'Number of samples within 3% = {s2_cnt}')
    print(f'Average solution value = {avg/1000}')
    print(f'Percentage of chains broken = {ncb_cnt/10}')
    # dwave.inspector.show(sample_set)
    f.write(f'{0},{i},{s1_cnt},{s2_cnt},{avg/1000},{ncb_cnt/10}\n')
    i = i - 0.02
    if (i<0.099):
        break
f.close()