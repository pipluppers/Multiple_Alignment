# Multiple_Alignment
Thanks to Dr. Tao Jiang for overseeing this project at UCR for CS 238: Algorithmic Techniques in Computational Biology. 

Input: 5x5 Scoring Matrix, 3 sequences with alphabet {A,C,G,T}

Output: 3 sequences that maximize the score

Generate a multiple alignment such that the score between the sequences is maximized. The scoring matrix is provided by the user and the scoring criteria is SP (Sum-Of-Pairs). The SP score is based on the maximum pairwise alignment of each combination of the nucleotides at each location in the sequences.

Dynamic Programming will be used to determine the score between the sequences. This program uses two different versions of dynamic programming:
1) one that runs normally in cubic time and space
2) one that runs in cubic time but linear space

The second algorithm uses divide and conquer to split up the table. When one of the sequences is small enough (i.e. has a size of 0 or 1 nucleotide), the first algorithm will begin to run and compute the alignment. The space-saving dynamic programming is necessary because modern genome sequences can be extremely long and memory will have a large impact on the runtime.
