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

EPS = 1e-9
cached : dict[any,float] = {}
cached_avg : float = {}

def YanStrength( source : BinaryQuadraticModel, embedding : EmbeddedStructure, multiplier: float):
    global cached
    print('YanStrength called')
    if (cached == {}):
        
        embedded = embedding.embed_bqm(source, 1000000)
        to_logical: dict = {}
        cs_array : dict = {}

        j_field : dict = {}
        h_field : dict = {}

        for i in source.variables:
            cs_array[i] = 0
            for j in embedding[i]:  
                to_logical[j] = i
                j_field[j] = 0
                h_field[j] = 0

        for u in embedded.quadratic.keys():
            if (embedded.quadratic[u] != -1000000):
                # print(u)
                j_field[u[0]] += abs(embedded.quadratic[u])
                j_field[u[1]] += abs(embedded.quadratic[u])
        
        for u in embedded.variables:
            h_field[u] = embedded.linear[u]
        
        for i in source.variables:
            print(f'calculating {i}')
            m = len(embedding[i])
            p2 = 2 ** (m)
            print(f'{m} {p2}')

            for mask in range(1, p2):
                # print(f'{i} : {mask} / {p2}')
                sum1 : int = 0
                sum2 : int = - model.linear[i]
                setSize : int = 0

                for j in range(m):
                    v = embedding[i][j]
                    setSize += (mask&1)
                    if (mask&1):
                        sum1 += h_field[v] - j_field[v]
                        sum2 += h_field[v]
                    else:
                        sum2 -= j_field[v]
                    mask = mask >> 1

                sum1 = abs(sum1)
                sum2 = abs(sum2)

                cs_array[i] = max(cs_array[i], min(sum1,sum2)/setSize) 
            cs_array[i] = cs_array[i]
        cached = deepcopy(cs_array)
        print(cached)
        for i in source.variables:
            cs_array[i] *= multiplier 
        return cs_array
    else:
        cs_array = deepcopy(cached)
        for i in source.variables:
            cs_array[i] *= multiplier
        return cs_array



cs_list = [0.7, 0.8, 0.9, 0.95, 1, 1.05, 1.1, 1.2, 1.3]
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

## this is the stupidest fucking workaround ever but it has to be done
## since passing random_seed from LFEC class doesn't work for some BS reason

def fe(S, T, **kwargs):
    return find_embedding(S, T, random_seed=123123)


composite = LazyFixedEmbeddingComposite(Sampler, find_embedding=fe)

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