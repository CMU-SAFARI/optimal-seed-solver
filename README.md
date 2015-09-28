# optimal-seed-solver

# Background
Optimal Seed Solver (OSS) is a dynamic-programming algorithm that finds the optimal seeds of a read,
which renders the minimum total seed frequency.

Seed selection is an important step for pigeonhole based seed-and-extend read mappers. In seed
selection, a read is broken into multiple non-overlapping substrings called seeds. Seeds are used as
anchors to index into the reference genome. To tolerate _e_ errors, a read is typically broken into
_e+1_ seeds such that by the pigeonhole principle, at least one seed is error free.

The overall frequency of the selected seeds has a direct impact on the performance of the mapper. If
seeds are frequent, the mapper has to verify many potential mappings through edit-distance
calculation, which is a time-consuming process.

OSS aims to find the least set of _e+1_ non-overlapping seeds from a read using dynamic programing
method. For details about the algorithm of OSS, please consult our manuscript online at: http://arxiv.org/abs/1506.08235

# 
