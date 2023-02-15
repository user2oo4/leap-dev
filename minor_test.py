import minorminer
import dwave_networkx as dnx
from dimod.generators import and_gate

from dwave.system import DWaveSampler
sampler_manual = DWaveSampler(solver={'topology__type': 'chimera'})

bqm_and1 = and_gate('in1', 'in2', 'out')
C1 = dnx.chimera_graph(1)
embedding_and = minorminer.find_embedding(list(bqm_and1.quadratic.keys()), C1)
print(embedding_and)

sampler_manual.sample_ising