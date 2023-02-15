import random
from dwave.cloud import Client

from dwave.system import DWaveSampler, EmbeddingComposite

from dimod.binary.binary_quadratic_model import BinaryQuadraticModel

import dwave.inspector

from dimod.vartypes import SPIN

import dwave_networkx as dnx

# D-Wave 2X
C = dnx.chimera_graph(12, 12, 4)

sus=[]

CSampler = DWaveSampler(solver={'topology__type': 'zephyr'})

for e in CSampler.edgelist:
    if ((e[0]>=0 and e[0]<512) and (e[1]>=0 and e[1]<512)):
        # print(e)
        sus.append(e)


random.shuffle(sus)

for e in sus:
    print(f'{e[0]} {e[1]}')