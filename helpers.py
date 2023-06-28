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
from climb_sample import hc_sample
from dwave.samplers import SimulatedAnnealingSampler


from climb_sample import hc_sample

from dwave.embedding import (target_to_source, unembed_sampleset, embed_bqm,
                             chain_to_quadratic, EmbeddedStructure)

SASampler = SimulatedAnnealingSampler()
# CSampler = DWaveSampler(solver={'topology__type': 'pegasus'})

hw = dnx.chimera_graph(16,16)

class HillClimbChimeraSampler(dimod.Sampler, dimod.Structured):
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
    
    def sample(self, bqm: BinaryQuadraticModel, num_reads: int, seed = -1, **parameters) -> SampleSet:
        
        print(f'start: {time.time()}')
        new_bqm, mapping = bqm.relabel_variables_as_integers(inplace=False)
        print(f'relabel: {time.time()}')
        state = random.getstate()
        if (seed == -1):
            random.seed()
        else:
            random.seed(seed)
        raw_samples : list[dict[int,int]] = []
        for i in range(num_reads):
            cur_seed = random.randint(0,998244353)
            # print(cur_seed)
            arr = hc_sample(new_bqm, cur_seed)
            # print(arr)
            dic : dict[int,int] = {}
            for (k,v) in mapping.items():
                # print(f'{k},{v}\n')
                dic[v] = arr[k]
            raw_samples.append(dic)
        random.setstate(state)
        
        print(f'end of runs: {time.time()}')
        return SampleSet.from_samples_bqm(raw_samples, bqm)
