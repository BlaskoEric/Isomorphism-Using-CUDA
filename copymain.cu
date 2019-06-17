#include<iostream>
#include"matrix.h"
#include<chrono>
#include<string>
#include<stdlib.h>
#include<vector>

using namespace std;
using namespace std::chrono;

const int N = 6;
bool print;
bool runAll;

void displayResults(matrix * g, matrix *h, matrix * can, int * d)
{
    if(print)
    {
        printf("G1 Matrix:\n");
        g->displayArray();

        printf("G2 Matrix: \n");
        h->displayArray();
    
        printf("Candidate Matrix:\n");
        can->displayArray();

        printf("D Matrix: \n");
        for(int i = 0; i < N; i++)
            printf("%d ",d[i]);     
        printf("\n");
    } 
}
 
void MakeMatrix(matrix* mat1,matrix* mat2)
{
    mat1->add_edge(0,1);
    mat1->add_edge(0,2);
    mat1->add_edge(0,4);
    mat1->add_edge(2,3);
    mat1->add_edge(2,4);
    mat1->add_edge(3,3);
    mat1->add_edge(3,4);
    mat1->add_edge(4,5);
    mat1->add_edge(5,5);

    mat2->add_edge(0,1);
    mat2->add_edge(0,3);
    mat2->add_edge(0,4);
    mat2->add_edge(1,2);
    mat2->add_edge(1,3);
    mat2->add_edge(3,4);
    mat2->add_edge(3,5);
    mat2->add_edge(4,4);
    mat2->add_edge(5,5);
}

__global__ void launchAlgTwo(matrix * dG, matrix * dH, matrix * dcandidate, int * dD, int * dfinalMatrix,int COUNT)
{
    int original = blockIdx.x*blockDim.x + threadIdx.x;
    int idx = original;

    int G1 = dcandidate->getNodeG(dfinalMatrix[idx]);

    while(idx < COUNT)
    {
        int update = 1;
        int G2 = dcandidate->getNodeH(dfinalMatrix[idx]);
        //printf("original = %d, idx = %d\n",original,idx);
        
        while(update)
        {
            update = 0;
            for(int k = 0; k < N; k++)
            {
                if(dD[k] != -1)
                {
                    printf("Original:%d, dG[%d,%d], dH[%d,%d]\n",original,G1,k,G2,dD[k]);
                    if(dG->getValueAtIndex(G1,k) != dH->getValueAtIndex(G2,dD[k]))
                    {
                        if(dcandidate->getValueAtIndex(G1,G2) == 1)
                            update = 1;
                        dcandidate->assignValue(G1,G2,0);
                        //printf("Original:%d, dG[%d,%d], dH[%d,%d]\n",original,G1,k,G2,dD[k]);
                    }
                }
            }
            if(dcandidate->nodeDegree(idx) == 1)
                dD[idx] = dcandidate->getJPos(idx);
            else if (dcandidate->nodeDegree(idx) == 0)
                dD[idx] = 0;
            else
                dD[idx] = -2;
            
        }
        idx += blockDim.x * gridDim.x;
    }
}


__host__ void startAlgTwo(matrix * dG, matrix * dH, matrix * dcandidate, int * dD)
{
    int count = 0;
    vector<int> finalMatrix;
    for(int i = 0; i < N;i++)
    {
        if(dcandidate->nodeDegree(i) > 1)
        {
            for(int j = 0; j < N; j++)
            {
                if(dcandidate->getValueAtIndex(i,j) == 1)
                {
                    finalMatrix.push_back(i*N+j);
                    count++;
                }
            }
        }
    }

    int *dfinalMatrix = 0;

    cudaMalloc((void**) &dfinalMatrix,finalMatrix.size()*sizeof(int));
    cudaMemcpy(dfinalMatrix,&finalMatrix[0],finalMatrix.size() *sizeof(int), cudaMemcpyHostToDevice);
    dim3 block = (ceil(N/256),1,1);
    //dim3 block = (1,1,1);
    launchAlgTwo<<<block,256>>>(dG,dH,dcandidate,dD,dfinalMatrix,count);
    cudaDeviceSynchronize();

    
     displayResults(dG,dH,dcandidate,dD);
}

__global__ void AlgOneGPU(matrix * dG, matrix * dH, matrix * dcandidate, int * dD)
{
    int idx = blockIdx.x*blockDim.x + threadIdx.x;
    if(idx < N)
    {
        int G1Degree = dG->nodeDegree(idx);
        bool recursiveNode = dG->recursiveNode(idx);
        for(int j = 0; j < N; j++)
        {
            if(G1Degree == dH->nodeDegree(j) && recursiveNode == dH->recursiveNode(j))
            {
                dcandidate->assignValue(idx,j,1);
            }
        }
        if(dcandidate->nodeDegree(idx) == 0)
            dD[idx] = -2;
        else if (dcandidate->nodeDegree(idx) == 1)
        {
            dD[idx] = dcandidate->getJPos(idx);
        }
        else
            dD[idx] = -1;
    }
}

__host__ void AlgOneStartCPU()
{
    printf("Starting Algorithm 1:\n");
    auto start = high_resolution_clock::now();
   
    matrix hG(N,N);
    matrix hH(N,N);
    matrix hcandidate(N,N);

    int hD[N];
    for(int i = 0; i < N; i++) hD[i] = -1;

    MakeMatrix(&hG,&hH);

    int gdegree[N];
    int hdegree[N];
    bool gself[N];
    bool hself[N];

    for(int i = 0; i < N; i++)
    {
        gdegree[i] = hG.nodeDegree(i);
        hdegree[i] = hH.nodeDegree(i);
        gself[i] = hG.recursiveNode(i);
        hself[i] = hH.recursiveNode(i);
    }        

    if (hG.nodeCount() != hH.nodeCount() || hG.edgeCount() != hH.edgeCount())
    {
        printf("Graphs are not isomophic\n");
        return;
    }
    for (int i = 0; i < N; i++)
    {
        //int G1Degree = hG.nodeDegree(i);
        //bool recursiveNode = hG.recursiveNode(i);    
        for(int j = 0; j < N; j++)
        {
            if (gdegree[i] == hdegree[j] && gself[i] == hself[j])
                hcandidate.assignValue(i,j,1);
        }
        if(hcandidate.nodeDegree(i) == 0)
            hD[i] = -2;
        else if(hcandidate.nodeDegree(i) == 1)
            hD[i] = hcandidate.getJPos(i);
        else
            hD[i] = -1;
    }

    displayResults(&hG,&hH,&hcandidate,hD);

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);

    printf("Time for CPU: %d\n",duration.count());

    if(runAll)
    {
        matrix *dG;
        matrix *dH;
        matrix *dcandidate;
        int *dD;

        cudaMallocManaged(&dG,sizeof(hG));
        cudaMallocManaged(&dH,sizeof(hH));
        cudaMallocManaged(&dcandidate,sizeof(hcandidate));
        cudaMallocManaged(&dD,N*sizeof(int));    

        *dG = hG;
        *dH = hH;
        *dcandidate = hcandidate;
        *dD = hD[0];
    
        for(int i = 0; i < N; i++)
        {
            dD[i] = hD[i];
        }
        startAlgTwo(dG,dH,dcandidate,dD);
        cudaFree(dD);
     }
}

__host__ void AlgOneStartGPU()
{
    printf("Starting Algorithm 1:\n");
    auto start = high_resolution_clock::now();
    
    matrix hG(N,N);
    matrix hH(N,N);
    matrix hcandidate(N,N);

    int hD[N];
    for(int i = 0; i < N; i++) hD[i] = -1;

    MakeMatrix(&hG,&hH);

    int gdegree[N];
    int hdegree[N];
    bool gself[N];
    bool hself[N];

    for(int i = 0; i < N; i++)
    {
        gdegree[i] = hG.nodeDegree(i);
        hdegree[i] = hH.nodeDegree(i);
        gself[i] = hG.recursiveNode(i);
        hself[i] = hH.recursiveNode(i);
    }           

    matrix *dG;
    matrix *dH;
    matrix *dcandidate;
    int *dD;

    cudaMallocManaged(&dG,sizeof(dG));
    cudaMallocManaged(&dH,sizeof(dH));
    cudaMallocManaged(&dcandidate,sizeof(dcandidate));
    cudaMallocManaged(&dD,N*sizeof(int));    

    *dG = hG;
    *dH = hH;
    *dcandidate = hcandidate;
    *dD = hD[0];
 
    if(dG->nodeCount() == dH->nodeCount() && dG->edgeCount() == dH->edgeCount())
    {
        dim3 block = (ceil(N/1024),1,1);
        AlgOneGPU<<<block,1024>>>(dG,dH,dcandidate,dD);
        cudaDeviceSynchronize();
    }
    else
    {
        printf("Graphs are not Isomophic\n");
    }
    displayResults(dG,dH,dcandidate,dD);
 
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    printf("Time for GPU: %d\n",duration.count());

    if(runAll)
        startAlgTwo(dG,dH,dcandidate,dD);
    cudaFree(dD);
}



int main(int argc,char** argv)
{
    bool CPU = false;
    bool GPU = false;

    if(argc == 1)
    {
        CPU = true;
        GPU = false;
        printf("Invalid arg. Using CPU method\n");
    }
    if(argc > 1)
    {
        std::string method = argv[1];
        if(method == "-gpu")
        { 
            GPU = true;
            CPU = false;
            printf("Using GPU method\n");
        }
        else
        {
            GPU = false;
            CPU = true;
            printf("Using CPU method\n");
        }
    }
    if(argc >= 3)
    {
        std::string Print = argv[2];
        if(Print == "-off")
            print = false;
        else 
            print = true;
    }
    else
    {
        print = true;
    }
    if(argc >= 4)
    {
        string run = argv[3];
        if(run == "-alg1")
            runAll = false;
        else
            runAll = true;
    }
    else
    {
        runAll = true;
    }
    
    if(CPU)
    {
        AlgOneStartCPU();
    } 
    if(GPU)
    {
        AlgOneStartGPU();   
    }

  return 0;
}


