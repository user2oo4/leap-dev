import random
from dwave.cloud import Client
from dwave.system import DWaveSampler, EmbeddingComposite
from dimod.binary.binary_quadratic_model import BinaryQuadraticModel
import dwave.inspector
from dimod.vartypes import SPIN
import dwave_networkx as dnx
import networkx as nx

from dwave.samplers import SimulatedAnnealingSampler
SASampler = SimulatedAnnealingSampler()


Sampler = SASampler

EPS = 1e-9

with Client.from_config() as client:
    CSolver = client.get_solver(name='DW_2000Q_6')
    Solver = CSolver

    cs_list = [10,25,50,100,250,500,1000,3000,10000]
    at_list = [5,10,20,50,100]

    

    # cs_list = [10,25,50]
    # at_list = [5,10,20]

    os_table = [[],[],[],[],[]]
    ncb_table = [[],[],[],[],[]]

    print('reading input')

    n = int(input())

    model = BinaryQuadraticModel(vartype=SPIN)

    for i in range(n):
        model.set_linear(i,int(input()))
    
    minn = 1000000
    minCnt = 0

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
    
    sampleset = Sampler.sample(model, num_reads = 1000)

    for r in sampleset.record:
        if (r.energy < minn):
            minn = r.energy
            minCnt = r.num_occurrences
        elif (r.energy == minn):
            minCnt += r.num_occurrences
    # print(sampleset)
    print(minn)
    print(minCnt)
