import random
from dwave.cloud import Client
from dwave.system import DWaveSampler, EmbeddingComposite
from dimod.binary.binary_quadratic_model import BinaryQuadraticModel
import dwave.inspector
from dimod.vartypes import SPIN
import dwave_networkx as dnx
from networkx import Graph
import networkx
from queue import PriorityQueue

seed : int = 15465
random.seed(seed)

from dwave_networkx.generators import chimera_graph

from dwave.samplers import SimulatedAnnealingSampler
SASampler = SimulatedAnnealingSampler()

Sampler = SASampler


def solve(graph: Graph) -> int :
    model = BinaryQuadraticModel(vartype=SPIN)
    print(graph)
    for n in graph.nodes.data():
        model.add_variable(n[0], n[1]['w'])
    for e in graph.edges.data():
        model.add_quadratic(e[0], e[1], e[2]['w'])
    result = SASampler.sample(bqm=model, num_reads=10000, seed=seed)
    minn = 0
    for r in result.record:
        if (minn > r.energy):
            minn = r.energy
    return minn

def output(graph: Graph, name: str):
    networkx.drawing.nx_agraph.write_dot(graph, f'dotfiles/{name}.dot')
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

n = 4
m = 2
cg : Graph = chimera_graph(n,m)

model = BinaryQuadraticModel(vartype=SPIN)

logical : dict[any, any] = {} 

pedge : PriorityQueue[tuple[int, int, int]] = PriorityQueue()

print("Generating initial graph")

for u in cg:
    w = random.randint(1,3)
    if (random.randint(0,1) == 1):
        w = -w
    model.add_variable(u,w)
    logical[u] = u

for (u,v) in cg.edges():
    
    w = random.randint(1,3)
    if (random.randint(0,1) == 1):
        w = -w
    
    model.add_quadratic(u, v, w)

    pedge.put((abs(w), u, v))

print("Sampling initial graph")

result = Sampler.sample(model, seed=seed, num_reads=10000)

print("Begin merging")

min_energy = 1000000

for r in result.record:
    min_energy = min(min_energy, r.energy)

print(min_energy)

nodecnt = n*m*8

while ((not(pedge.empty()))and(nodecnt > 25)) :
    tup = pedge.get()

    u = tup[1]
    v = tup[2]

    u_ind = result.variables.index(u)
    v_ind = result.variables.index(v)

    fail = 0
    for r in result.record:
        if ((r.energy == min_energy) and (r.sample[u_ind] != r.sample[v_ind])):
            fail = 1

    # Yep, this is unoptimal as fuck. Too bad!
    if ((fail == 0) and (logical[u] != logical[v])):
        fail = 1
        ua = logical[u] 
        va = logical[v]
        a = min(ua,va)
        nodecnt -= 1
        for i in logical.keys():
            if ((logical[i] == ua)or(logical[i] == va)):
                logical[i] = a
    
print("Constructing new graph")

toInd : dict[any,any] = {}

for i in logical.values():
    toInd[i] = -1

cnt = 0

for i in logical.values():
    if (toInd[i] == -1):
        toInd[i] = cnt
        cnt += 1

for i in logical.keys():
    logical[i] = toInd[logical[i]]

j: dict[tuple[any,any], any] = {}
h: dict[any, any] = {}

for i in model.linear.keys():
    h[logical[i]] = 0
for i in model.linear.keys():
    h[logical[i]] += model.linear[i]

    
for u in model.quadratic.keys():
    if (logical[u[0]] != logical[u[1]]):
        j[(min(logical[u[0]],logical[u[1]]), max(logical[u[0]], logical[u[1]]))] = 0
for u in model.quadratic.keys():
    if (logical[u[0]] != logical[u[1]]):
        j[(min(logical[u[0]],logical[u[1]]), max(logical[u[0]], logical[u[1]]))] += model.quadratic[u]

print("Finalizing new graph")

finalGraph = Graph()

for i in h.keys():
    finalGraph.add_node(i, w = h[i])

for u in j.keys():
    finalGraph.add_edge(u[0], u[1], w = j[u])

print("Outputting new graph")

output(finalGraph, "rm25")
