# optimal-seed-solver

## Background:
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
method. For details about the algorithm of OSS, please consult our manuscript online at:
http://arxiv.org/abs/1506.08235

## Code files:

This repository contains many code files. Besides hosting the OSS code, it also contains other seed
selection implementations, code files for main classes and misc code files. Below we summarize each code file.

### Seed selection mechanisms:

+ **optimalSolverLN.h, optimalSolverLN.cc**: contains the code of OSS.
+ **optimalSolver.h, optimalSolver.cc**: contains the code of an optimal solver that has
quadruple-complexity. We used this code to verify the integrity of OSS.
+ **basicSolver.h, basicSolver.cc**: contains the code of a basic solver, which selects seeds at fix
positions with fixed lengths.
+ **fastHASHSolver.h, fastHASHSolver.cc**: contains the code simulating the fastHASH seed selection
mechanism.
+ **hobbesSolver.h, hobbesSolver.cc**: contains the code simulating the hobbes seed selection
mechanism.
+ **spacedSeedSolver.h, spacedSeedSolver.cc**: contains the code simulating the spaced seed
selection mechanism.
+ **thresholdSolver.h, thresholdSolver.cc**: contains the code simulating the Gem mapper's seed
selection mechanism.

### MISC code files:

+ **KmerHash.h, KmerHash.cc**: a hash function that transforms a ASCII string into a 64-bit int number.
+ **HashTree.h, HashTree.cc**: a suffix tree to index the reference genome. Notice that this
implementation is **not optimized for performance**. It was initially built for gathering statistics of
the reference genome. For human reference genome it needs 350 GB main memory for indexing. Please
consider replacing HashTree with BWT implementations if you plan to integrate OSS into a production
mapper.
+ **RefDB.h, RefDB.cc**: a database that caches the reference genome.

### Main classes:

Each seed selection mechanism has a complementary main-code file. They are named in a fashion of
"test(mechanism).cc".

## How to run:

There are two steps to run OSS (and similarly other seed selection mechanisms):

1. Index the reference genome:
`$ ./testHashTree referenceName.fasta treeFile.tree`

2. Run OSS:
`$ ./testOptimalSeedLN treeFile.tree readFile.fastq`
