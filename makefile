all: kcliques

kcliques:
	/usr/local/cuda/bin/nvcc --optimize 3 --gpu-architecture=compute_70 -o kcliques kcliques.cu
