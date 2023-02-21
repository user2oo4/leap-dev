import random
from dwave.cloud import Client
from dwave.system import DWaveSampler, EmbeddingComposite
from dimod.binary.binary_quadratic_model import BinaryQuadraticModel
import dwave.inspector
from dimod.vartypes import SPIN
import dwave_networkx as dnx
import networkx as nx
from minorminer import find_embedding
from dwave.embedding import embed_bqm

CSampler = DWaveSampler(solver={'topology__type': 'chimera'})
Sampler = CSampler
hardware_graph = Sampler.edgelist
nx_graph = nx.Graph(hardware_graph)

EPS = 1e-9

with Client.from_config() as client:
    CSolver = client.get_solver(name='DW_2000Q_6')
    Solver = CSolver

    cs_list = [25]
    at_list = [100]

    

    # cs_list = [10,25,50]
    # at_list = [5,10,20]

    os_table = [[],[],[],[],[]]
    ncb_table = [[],[],[],[],[]]

    print('reading input')

    n = int(input())

    model = BinaryQuadraticModel(vartype=SPIN)

    for i in range(n):
        model.set_linear(i,0)

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
    logical_graph = model.quadratic.keys()
    print('finding embedding')
    embedding : dict = find_embedding(logical_graph,hardware_graph)
    print(embedding)
    chain_edge_cnt = 0
    for r in embedding.values():
        chain_edge_cnt = chain_edge_cnt + len(r) - 1

    for cs in cs_list:
        embedded_model = embed_bqm(model, embedding, nx_graph, cs)
        optimal_energy = -cs * chain_edge_cnt + optimal_solution

        for j in range(len(at_list)):
            at = at_list[j]
            print('Current run:')
            print(f'chain strength = {cs}')
            print(f'annealing time = {at}')
            print(f'optimal solution = {optimal_solution}')
            print(f'optimal energy = {optimal_energy}')

            sample_result = CSolver.sample_bqm(embedded_model, num_reads=100, annealing_time = at)
            sample_result.wait_sampleset()
            sample_set = sample_result.sampleset

            os_cnt = 0
            ncb_cnt = 0
            print(sample_set)
            for r in sample_set.record:
                # print(dir(r))
                sus = 1
                for logical in embedding.values():
                    for i in range(len(logical)-1):
                        u = sample_set.variables.index(logical[i])
                        v = sample_set.variables.index(logical[i+1])
                        if not (r.sample[u] == r.sample[v]):
                            sus = 0
                if sus:
                    ncb_cnt += r.num_occurrences
                    print(r.energy, end=' ')
                    if (abs(r.energy - optimal_solution) < EPS):
                        os_cnt += r.num_occurrences

            print(f'Number of optimal samples = {os_cnt}')
            print(f'Number of no chain breaks = {ncb_cnt}')

            # os_table[j].append(os_cnt)
            # ncb_table[j].append(ncb_cnt)

            dwave.inspector.show_bqm_response(model, dict(embedding=embedding, chain_strength=cs), sample_result)
    # for r in os_table :
    #     print(r)
    # print()
    # print()
    # for r in ncb_table :
    #     print(r)