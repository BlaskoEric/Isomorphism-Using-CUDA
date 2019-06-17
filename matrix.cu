#include <iostream>
#include <cstdio>
#include "matrix.h"

        matrix::matrix(){
            height = 1;
            width = 1;
            sizeArray = height*width;
            cudaError_t err = cudaMallocManaged(&array,sizeArray*sizeof(array[0]));
            if (err != cudaSuccess)
            {
                //cout << "Memory allocation failed"<<endl;
                printf("Memory allocation failed");
            }
        }

        matrix::matrix(size_t h){
            height = h;
            width = 1;
            sizeArray = height*width;
            cudaError_t err = cudaMallocManaged(&array,sizeArray*sizeof(array[0]));
            if (err != cudaSuccess)
            {
                //cout << "Memory allocation failed"<<endl;
                printf("Memory allocation failed");
            }
        }

        matrix::matrix(size_t h,size_t w,int value){
            height = h;
            width = w;
            sizeArray = height*width;
            cudaError_t err = cudaMallocManaged(&array,sizeArray*sizeof(array[0]));

            if (err != cudaSuccess)
            {
                //cout << "Memory allocation failed"<<endl;
                printf("Memory allocation failed");
            }
            for(int i = 0; i < height*width; i++)
                array[i] = value;

        }

        matrix::matrix(const matrix &mat){
            height = mat.height;
            width = mat.width;
            sizeArray = mat.sizeArray;
            cudaError_t err = cudaMallocManaged(&array,sizeArray*sizeof(array[0]));
            if (err != cudaSuccess)
            {
                //cout << "Memory allocation failed"<<endl;
                printf("Memory allocation failed");
            }

            for(size_t i = 0;i<sizeArray;++i){
                array[i] = mat.array[i];
            }

        //copy(mat.array,mat.array+mat.sizeArray,array);
        }

        matrix &matrix::operator=(const matrix &mat){
            height = mat.height;
            width = mat.width;
            sizeArray = mat.sizeArray;
            cudaError_t err = cudaMallocManaged(&array,sizeArray*sizeof(array[0]));
            if (err != cudaSuccess)
            {
                //cout << "Memory allocation failed"<<endl;
                printf("Memory allocation failed");
            }

            for(size_t i = 0;i<sizeArray;++i){
                array[i] = mat.array[i];
            }

            //copy(mat.array,mat.array+mat.sizeArray,array);
            return *this;
        }

        matrix &matrix::operator=(const matrix *mat){
            height = mat->height;
            width = mat->width;
            sizeArray = mat->sizeArray;
            cudaError_t err = cudaMallocManaged(&array,sizeArray*sizeof(array[0]));
            if (err != cudaSuccess)
            {
                //cout << "Memory allocation failed"<<endl;
                printf("Memory allocation failed");
            }

            for(size_t i = 0;i<sizeArray;++i){
                array[i] = mat->array[i];
            }

            //copy(mat.array,mat.array+mat.sizeArray,array);
            return *this;
        }   

        matrix::~matrix(){
            cudaFree(array);
        }
        void matrix::assignValue(int i,int j, bool value){
            int l = i*width + j;
            array[l] = value;
        }

        void matrix::assignValue(int l, bool value){
            array[l] = value;
        }

        void matrix::displayArray(){
            size_t i,j,l;
            for(i=0;i<height;++i){
                for(j=0;j<width;++j){
                    l =i*width + j;
                    //cout<<array[l]<<"\t";
                    printf("%d ",array[l]);
                }
                //cout<<endl;
                printf("\n");
            }
            printf("\n");
        }

        void matrix::add_edge(int x, int y)
        {
            if(x > sizeArray || y > sizeArray)
                printf("Invalid edge!\n");
            else
            {
                array[x*width + y] = 1;
                array[y*width + x] = 1;
            }
        }

        int matrix::nodeCount()
        {
            return sizeArray;
        }

        int matrix::edgeCount()
        {
            int count = 0;
            int i,j;
            for(i = 0; i < width; i++)
            {
                for(j = 0; j < height; j++)
                {
                    if(array[i*width + j] == 1)
                        count++;
                }
            }
            return count;
        }
    
        int matrix::nodeDegree(int d)
        {
            int count = 0;
            int i;
            for(i = 0; i < width; i++)
            {
                if(array[d*width + i] == 1)
                    count++;
            }
            return count;
        }

        bool matrix::recursiveNode(int d)
        {
            if(array[d*width + d] == 1)
                return true;
            else
                return false;
        }

        bool matrix::getValueAtIndex(int i, int j)
        {
            return array[i*width+j];
        }

        int matrix::getJPos(int idx)
        {
            for(int i = 0; i < width; i++)
            {
                if(array[idx*width+i] == 1)
                    return i;
            }
            return 0;
        }

        int matrix::getNodeG(int i)
        {
            return i/width;
        }
    
        int matrix::getNodeH(int i)
        {
            return i%width;
        }

