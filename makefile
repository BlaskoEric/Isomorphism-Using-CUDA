Iso: main.o matrix.o
	nvcc -o Iso main.o matrix.o

main.o:	main.cu
	nvcc -dc main.cu

matrix.o: matrix.cu
	nvcc -dc matrix.cu
