#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <stdio.h>
#include "common/errors.h"
#define MAX_DEGREE 1024
#define MAX_DEPTH 12
#define OFFSET (MAX_DEGREE + 2) * MAX_DEGREE / 64

__constant__ int offsets[MAX_DEGREE];
__global__ void sort_edges(int *input, int *input_vertex_places)
{
    __shared__ int to_sort[MAX_DEGREE];
    __shared__ int to_sort2[MAX_DEGREE];

    int main_vertex = blockIdx.x;
    int start_main_edge = input_vertex_places[main_vertex];
    int end_main_edge = input_vertex_places[main_vertex + 1];
    int len = end_main_edge - start_main_edge;

    if (threadIdx.x < end_main_edge - start_main_edge)
    {
        to_sort[threadIdx.x] = input[start_main_edge + threadIdx.x];
        printf("%d %d %d %d\n", main_vertex, start_main_edge, end_main_edge, threadIdx.x);
    }
    __syncthreads();
    for (int i = 2; i < 2 * len; i *= 2)
    {
        if (threadIdx.x * i < len)
        {
            int start = threadIdx.x * i;
            int end = min(start + i, len);
            int mid = min(start + i / 2, len);
            int mid2 = start + i / 2;
            int count = threadIdx.x * i;
            while (mid2 < end || start < mid)
            {
                if (start < mid && (mid2 == end || to_sort[start] < to_sort[mid2]))
                {
                    to_sort2[count] = to_sort[start];
                    start++;
                }
                else
                {
                    to_sort2[count] = to_sort[mid2];
                    mid2++;
                }
                count++;
                printf("%d %d %d %d \n", count, mid2, start, len);
            }
        }
        __syncthreads();
        if (threadIdx.x < len)
        {
            to_sort[threadIdx.x] = to_sort2[threadIdx.x];
        }
        __syncthreads();
    }
    if (threadIdx.x < end_main_edge - start_main_edge)
    {
        input[start_main_edge + threadIdx.x] = to_sort[threadIdx.x];
    }
}

__global__ void kcliques(int k, int *input, int *input_vertex_places, int *output)
{
    extern __shared__ int sh_mem[];
    int main_vertex = blockIdx.x;
    // printf("%d \n",4 * MAX_DEGREE * MAX_DEGREE /32);
    // return;

    //__shared__ int incidence_matrix[MAX_DEGREE][MAX_DEGREE / 32];
    //__shared__ int main_input[MAX_DEGREE]; // copy there

    int start_main_edge = input_vertex_places[main_vertex];
    int end_main_edge = input_vertex_places[main_vertex + 1];
    int len = end_main_edge - start_main_edge;
    if (threadIdx.x < len)
    {
        sh_mem[OFFSET + threadIdx.x] = input[start_main_edge + threadIdx.x];
        for (int i = 0; i <= len / 32; i++)
        {
            sh_mem[offsets[threadIdx.x] + i] = 0;
            printf("%d %d \n", threadIdx.x, i);
        }
    }
    __syncthreads();
    int actual_main_edge = 0;
    // printf("FOG");
    if (threadIdx.x < end_main_edge - start_main_edge)
    {
        printf("%d %d %d %d\n", main_vertex, start_main_edge, end_main_edge, threadIdx.x);

        int start_my_edge = input_vertex_places[input[start_main_edge + threadIdx.x]];   // if x to big don't do
        int end_my_edge = input_vertex_places[input[start_main_edge + threadIdx.x] + 1]; // if x to big don't do

        printf("xd %d %d %d %d\n", main_vertex, start_my_edge, end_my_edge, threadIdx.x);
        // return;
        for (int i = start_my_edge; i < end_my_edge; i++)
        {
           // printf("gowno");
            while (actual_main_edge < len && sh_mem[OFFSET + actual_main_edge] < input[i])
            {
                actual_main_edge++;
            }
            if (actual_main_edge < len && sh_mem[OFFSET + actual_main_edge] == input[i])
            {
                // printf("%d %d %d\n", main_vertex, threadIdx.x, actual_main_edge);
           //    printf("%d %d %d \n", threadIdx.x, actual_main_edge / 32, (1 << (actual_main_edge % 32)));

                // return;
                sh_mem[offsets[threadIdx.x] + actual_main_edge / 32] += (1 << (actual_main_edge % 32)); // to check move
            }
        }
        int stack[MAX_DEPTH];
        int res[MAX_DEPTH];
        for (int i = 0; i < k; i++)
        {
            res[i] = 0;
        }
        int stack_pointer = 2;
        int actual_stack_values[MAX_DEPTH][MAX_DEGREE / 32], pos = MAX_DEGREE / 32 - 1;
        for (int i = offsets[threadIdx.x + 1]; i >= offsets[threadIdx.x]; i--)
        {
            actual_stack_values[1][pos] = sh_mem[i];
            pos--;
        }
                    
        stack[0] = main_vertex;
        stack[1] = threadIdx.x;
        stack[2] = threadIdx.x + 1;
        while (stack_pointer >= 2)
        {
            printf("%d %d %d %d\n", stack[0], stack[1], stack[2], stack[3]);
            stack[stack_pointer]++;
            if (stack[stack_pointer] > MAX_DEGREE)
            {
                stack_pointer--;
                // for (int i = 0; i < stack_pointer; i++)
                // {
                //     for (int z = stack[stack_pointer] / 32; i < MAX_DEGREE / 32; i++)
                //     {
                //         if (i == 0)
                //         {
                //             actual_stack_values[threadIdx.x][z] = incidence_matrix[stack[i]][z];
                //         }
                //         else
                //         {
                //             actual_stack_values[threadIdx.x][z] &= incidence_matrix[stack[i]][z];
                //         }
                //     }
                // }
            }
            if (actual_stack_values[threadIdx.x][stack[stack_pointer] / 32] & 1 << stack[stack_pointer] % 32) // if 1 is on stack_ptr
            {
                res[stack_pointer]++;
                if (stack_pointer < k - 1)
                {
                    stack_pointer++;
                    int posx = MAX_DEGREE / 32 - 1;
                    for (int i = offsets[stack[stack_pointer] + 1]; i >= offsets[stack[stack_pointer]]; i--)
                    {
                        actual_stack_values[stack_pointer][posx] = sh_mem[i] & actual_stack_values[stack_pointer - 1][posx];
                        posx--;
                    }
                    // for (int i = stack[stack_pointer] / 32; i < MAX_DEGREE / 32; i++)
                    // {
                    //     actual_stack_values[stack_pointer][i] &= incidence_matrix[stack[stack_pointer]][i];
                    // }
                    stack[stack_pointer] = stack[stack_pointer - 1] + 1;
                }
            }
        }
        for (int i = 0; i < k; i++)
        {
            atomicAdd(&output[i], res[i]); // to optimize
        }
    }
}

int main(int argc, char *argv[])
{
    if (argc != 4)
    {
        std::cout << "Usage ./kcliques <graph input file> <k value> <output file>";
        return 1;
    }
    int k = stoi(std::string(argv[2]));
    std::ifstream input;
    input.open(std::string(argv[1]));
    unsigned int a, b, vertex_num = 0;
    std::vector<std::pair<int, int>> edges_list;
    std::vector<std::vector<int>> edges;
    std::vector<int> edges_sorted, vertex_places;

    while (input >> a >> b)
    {
        if (a > b)
        {
            std::swap(a, b);
        }
        edges_list.emplace_back(a, b);
        vertex_num = std::max(vertex_num, b);
    }
    for (int i = 0; i < vertex_num; i++)
    {
        edges.emplace_back();
    }
    for (int i = 0; i < edges_list.size(); i++)
    {
        edges[edges_list[i].first].push_back(edges_list[i].second);
    }
    int counter = 0;
    for (int i = 0; i < edges.size(); i++)
    {
        vertex_places.push_back(counter);
        for (int z = 0; z < edges[i].size(); z++)
        {
            edges_sorted.push_back(edges[i][z]);
        }
        counter += edges[i].size();
    }
    int output[k];
    for (int i = 0; i < k; i++)
    {
        output[i] = 0;
    }
    vertex_places.push_back(counter);
    int *edges_gpu, *vertex_places_gpu, *output_gpu;
    HANDLE_ERROR(cudaMalloc(&edges_gpu, sizeof(int) * edges_sorted.size()));
    HANDLE_ERROR(cudaMalloc(&vertex_places_gpu, sizeof(int) * (vertex_num + 1)));
    HANDLE_ERROR(cudaMalloc(&output_gpu, sizeof(int) * k));
    HANDLE_ERROR(cudaMemcpy(edges_gpu, edges_sorted.data(), sizeof(int) * edges_sorted.size(), cudaMemcpyHostToDevice));
    HANDLE_ERROR(cudaMemcpy(vertex_places_gpu, vertex_places.data(), sizeof(int) * (vertex_num + 1), cudaMemcpyHostToDevice));
    HANDLE_ERROR(cudaMemcpy(output_gpu, output, sizeof(int) * k, cudaMemcpyHostToDevice));
    // HANDLE_ERROR(cudaDeviceSynchronize());
    int constants_to_gpu[MAX_DEGREE], mem_offset = MAX_DEGREE;
    constants_to_gpu[0] = 0;
    for (int i = 1; i < MAX_DEGREE; i++)
    {
        constants_to_gpu[i] = constants_to_gpu[i - 1] + (mem_offset + 32) / 32;
        mem_offset--;
    }
    cudaMemcpyToSymbol(offsets, constants_to_gpu, sizeof(int) * MAX_DEGREE);
    sort_edges<<<vertex_num, MAX_DEGREE>>>(edges_gpu, vertex_places_gpu); // optimize to infinity
                                                                          // cudaMemcpy(edges_sorted.data(), edges_gpu, sizeof(int) * edges_sorted.size(), cudaMemcpyDeviceToHost);

    // for(int i=0;i<edges_sorted.size();i++)
    //{
    //     output_stream<<edges_sorted[i]<<" ";
    // }
    // output_stream<<"\n";
    // return 0;
    HANDLE_ERROR(cudaDeviceSynchronize());

    printf("%d", vertex_num);
    printf("xddd %ld \n", sizeof(int) * (MAX_DEGREE * MAX_DEGREE / 32 + MAX_DEGREE));
    int maxbytes = 98304;
    cudaFuncSetAttribute(kcliques, cudaFuncAttributeMaxDynamicSharedMemorySize, maxbytes);
    kcliques<<<vertex_num, MAX_DEGREE, maxbytes>>>(k, edges_gpu, vertex_places_gpu, output_gpu);
    // HANDLE_ERROR(cudaMemcpy(output, output_gpu, sizeof(int) * k, cudaMemcpyDeviceToHost));
    HANDLE_ERROR(cudaDeviceSynchronize());
    printf("%d", k);

    std::ofstream output_stream;
    output_stream.open(std::string(argv[3]));
    for (int i = 0; i < k; i++)
    {
        output_stream << output[i] << " ";
    }
}