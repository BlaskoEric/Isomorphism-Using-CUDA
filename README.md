# Isomorphism-Using-CUDA

## Acknowledgement
* This program was a project for parallel algorithms. It is based off of the
research from Min-Young-Son, Young-Hak Kim, Byoung-Woo Oh at the Dept. of 
Computer Engineering, Kumoh National Institute of Technology. The original 
paper can be found at https://www.researchgate.net/publication/287222374_An_Efficient_Parallel_Algorithm_for_Graph_Isomorphism_on_GPU_using_CUDA

## Description
* This program detects if two graphs are isomorphic. The program was implemented
using two methods: CPU/GPU and GPU/GPU. The original research paper used a hybrid
combination of CPU and GPU to test for isomorphism. I implemented it in the same
fashin as to get identical results when using large graphs. I then implemented a 
new version which used the GPU for both parts of the algorithm which resulted in
roughly a quarter speed up with graphs larger than 5000 nodes.

### Algorithm Part 1
* The first part of this algorithm reduces the total number of candidates by determining
if nodes from both graphs contain the same number of edges and if they have an edge to themselves. 
If two nodes do not have the same number of edges or are different in having an edge to itself,
it is removed from a candidate matrix. The final candidate matrix is used in the GPU to perform
the final calculation.

### Algorithm Part 2
* The second part is performed in the GPU using the candidate matrix from part 1. In this 
process, the nodes are compared in finer detail. If a difference is detected between two nodes
it is removed from the candidate matrix. Once complete, if each candidate has a remaining value
of one, then the two graphs are isomorphic. 
