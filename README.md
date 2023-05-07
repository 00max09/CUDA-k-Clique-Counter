# HPC-k-Clique-Counter
k Clique Counting CUDA implementation


# Compilation for Nvida Titan V

nvcc kcliques.cu --optimize 3 --gpu-architecture=compute_70 -ccbin "D:\Program Files\Microsoft Visual Studio\2022\Community\VC\Tools\MSVC\14.30.30705\bin\Hostx64\x64";