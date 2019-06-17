#ifndef MATRIX
#define MATRIX

#include <iostream>

class matrix
{

private:
    std::size_t height,width,sizeArray;
    bool *array;

public:
           __host__ matrix();
           __host__ matrix(std::size_t);
           __host__ matrix(std::size_t,std::size_t,int);
           __host__ matrix(const matrix &);
           __host__ matrix &operator=(const matrix &mat);
           __host__ matrix &operator=(const matrix *mat);
           __host__ ~matrix();
__device__ __host__ int nodeCount();
__device__ __host__ int edgeCount();
__device__ __host__ bool getValueAtIndex(int,int);
           __host__ void add_edge(int,int);
__device__ __host__ int getNodeG(int);
__device__ __host__ int getNodeH(int);
__device__ __host__ int getJPos(int);
__device__ __host__ int nodeDegree(int);
__device__ __host__ bool recursiveNode(int);
__device__ __host__ void assignValue(int,int,bool);
__device__ __host__ void assignValue(int,bool);
__device__ __host__ void displayArray();

};

#endif
