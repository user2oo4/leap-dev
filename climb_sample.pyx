# distutils: language=c++
# cython: language_level=3

from dimod.binary.binary_quadratic_model import BinaryQuadraticModel
from libcpp.set cimport set as cset
from libcpp.vector cimport vector as cvector 
from cython.operator cimport dereference as deref, preincrement as inc
import random

cdef extern from "raw_climb_sample.cpp":
    void raw_hc_sample(
            int n,
            int seed,
            int* h,
            (cvector[int])* j,
            (cvector[int])* weight,
            int* var) 


def hc_sample(bqm: BinaryQuadraticModel, seed: int):
    state = random.getstate()
    random.seed(seed)
    cdef int[2048] var, h, delta
    cdef (cvector[int])[2048] j, weight
    # print('Hello from cython')
    for (k,v) in bqm.linear.items():
        h[k] = v
    for ((u,v),w) in bqm.quadratic.items():
        # print(f'{u} {v} {w}')
        j[u].push_back(v)
        weight[u].push_back(w)
        
        j[v].push_back(u)
        weight[v].push_back(w)

    raw_hc_sample(bqm.num_variables, seed, h, j, weight, var)

    # for i in range(bqm.num_variables):
    #     var[i] = random.randint(0,1)*2-1

    # for i in range(bqm.num_variables):
    #     delta[i] = h[i] * var[i]
    #     for k in range(j[i].size()):
    #         delta[i] += var[i] * var[j[i][k]] * weight[i][k]
    #     if (delta[i] > 0):
    #         pos.insert(i)
    
    # cdef cset[int].iterator it
    # while (pos.size() != 0):
    #     it = pos.begin()
    #     tar = deref(it)
    #     pos.erase(it)
    #     delta[tar] = -delta[tar]
    #     var[tar] = -var[tar]
    #     for i in range(j[tar].size()):          
    #         delta[j[tar][i]] += 2 * var[tar] * var[j[tar][i]] * weight[tar][i]
    #         if (delta[j[tar][i]] <= 0):
    #             pos.erase(j[tar][i])
    #         else:
    #             pos.insert(j[tar][i])

    result : list[int] = []
    for i in range(bqm.num_variables):
        if (i<=10):
            print(f'{i}, {var[i]}')
        result.append(var[i])
    random.setstate(state)
    return result