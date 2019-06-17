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

__global__ void AlgOneGPU(int * d_gDegree,int * d_hDegree,bool * d_gSelfEdge,bool * d_hSelfEdge,
                          int * dD, int * d_candidate)
{
    int idx = blockIdx.x*blockDim.x + threadIdx.x;
    int candidateCount = N;
    if(idx < N)
    {   
        int G1Degree = d_gDegree[idx];
        bool recursiveNode = d_gSelfEdge[idx];
        for(int j = 0; j < N; j++)
        {
            if(G1Degree != d_hDegree[j] || recursiveNode != d_hSelfEdge[j])
            {
                candidateCount -= 1;
                d_candidate[idx*N+j] = 0;
            }
        }
        if(candidateCount == 0)
            dD[idx] = -2;
        else if (candidateCount == 1)
        {
            for(int i = 0; i < N; i++)    
            {
                if(d_candidate[idx*N+i] == 1)
                    dD[idx] = i;
            }
        }
        else
            dD[idx] = -1;
    }
}

__host__ void AlgOneStartCPU()
{
    printf("Starting Algorithm 1:\n");
    auto start = high_resolution_clock::now();
   
    matrix hG(N,N,0);
    matrix hH(N,N,0);
    matrix hcandidate(N,N,1);
    int hD[N];

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
        hD[i] = -1;
    }        

    if (hG.nodeCount() != hH.nodeCount() || hG.edgeCount() != hH.edgeCount())
    {
        printf("Graphs are not isomophic\n");
        return;
    }
    for (int i = 0; i < N; i++)
    {
        int G1Degree = gdegree[i];
        bool recursiveNode = gself[i];    
        for(int j = 0; j < N; j++)
        {
            if (G1Degree != hdegree[j] || recursiveNode != hself[j])
                hcandidate.assignValue(i,j,0);
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
    
    matrix hG(N,N,0);
    matrix hH(N,N,0);
    matrix hcandidate(N,N,1);
    MakeMatrix(&hG,&hH);
   
    int gCount = 0;
    int hCount = 0;
    int *d_gDegree;
    int *d_hDegree;
    bool *d_gSelfEdge;
    bool *d_hSelfEdge;
    int *dD;
    int *d_candidate;
       
    cudaMallocManaged( &d_gDegree,N*sizeof(int));
    cudaMallocManaged( &d_hDegree,N*sizeof(int));
    cudaMallocManaged( &d_gSelfEdge,N*sizeof(bool));
    cudaMallocManaged( &d_hSelfEdge,N*sizeof(bool));
    cudaMallocManaged( &d_candidate,N*N*sizeof(int));
    cudaMallocManaged( &dD,N*sizeof(int));
   
    for(int i = 0; i < N; i++)
    {
        for(int j = 0; j < N; j++)
            d_candidate[i*N+j] = 1;
    }
        
    for(int i = 0; i < N; i++)
    {
        d_gDegree[i] = hG.nodeDegree(i);
        d_hDegree[i] = hH.nodeDegree(i);
        d_gSelfEdge[i] = hG.recursiveNode(i);
        d_hSelfEdge[i] = hH.recursiveNode(i);
        dD[i] = -1;
        gCount += hG.nodeDegree(i);
        hCount += hH.nodeDegree(i);
    }           
    if(gCount == hCount)
    {
        dim3 block = (ceil(N/256),1,1);
        AlgOneGPU<<<block,256>>>(d_gDegree,d_hDegree,d_gSelfEdge,d_hSelfEdge,dD,d_candidate);
        cudaDeviceSynchronize();
        
        cudaFree(d_gDegree);
        cudaFree(d_hDegree);
        cudaFree(d_gSelfEdge);
        cudaFree(d_hSelfEdge);
     }
    else
    {
        printf("Graphs are not Isomophic\n");
    }
    //displayResults(hG,hH,dcandidate,dD);

    if(print)
    {
        printf("G1 Matrix:\n");
        hG.displayArray();

        printf("G2 Matrix: \n");
        hH.displayArray();


        printf("Candidate Matrix:\n");
        for(int i = 0; i < N;i++)
        {
            for(int j = 0; j < N;j++)
            {
                printf("%d ",d_candidate[i*N+j]);
            }
            printf("\n");
        }

        printf("D Matrix: \n");
        for(int i = 0; i < N; i++)
            printf("%d ",dD[i]);     
        printf("\n"); 
    }
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    printf("Time for GPU: %d\n",duration.count());
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
        std::string run = argv[3];
        if(run == "-alg1")
            runAll = false;
        else
            runAll = true;
    }
    else
    {
        runAll = true;
    }
   
    printf("here"); 
    if(CPU)
    {
        AlgOneStartCPU();
    } 
    if(GPU)
    {
        printf("Here");
        AlgOneStartGPU();   
    }

  return 0;
}


