#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <stdio.h>
#include <unordered_map>
#include <algorithm>
#include <map>
#include "common/errors.h"
#define MAX_DEGREE 1024
#define MAX_DEPTH 12
#define OFFSET (MAX_DEGREE + 2) * MAX_DEGREE / 64
#define MOD 1000000000

__constant__ int offsets[MAX_DEGREE + 1];

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
                if (start < mid && (mid2 >= end || to_sort[start] < to_sort[mid2]))
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
        if(threadIdx.x !=0)
        {
            if(to_sort[threadIdx.x]<to_sort[threadIdx.x-1])
            {
                printf("SEFRR");
            }
        }
    }
    
}

__global__ void kcliques(int k, int *input, int *input_vertex_places, unsigned long long int *output)
{
    extern __shared__ int sh_mem[];
    int main_vertex = blockIdx.x;
    
    int start_main_edge = input_vertex_places[main_vertex];
    int end_main_edge = input_vertex_places[main_vertex + 1];
    int len = end_main_edge - start_main_edge;
    int used = (len - 1) / 32 + 1;
    int not_used = max(MAX_DEGREE / 32 - used, 0);
    if (threadIdx.x < len)
    {
        sh_mem[OFFSET + threadIdx.x] = input[start_main_edge + threadIdx.x];
        for (int i = offsets[threadIdx.x]; i < offsets[threadIdx.x + 1]; i++)
        {
            sh_mem[i] = 0;
        }
    }
  
    __syncthreads();
    int actual_main_edge = threadIdx.x + 1 ;
    unsigned long long int res[MAX_DEPTH];
        
    if (threadIdx.x < end_main_edge - start_main_edge)
    {

        int start_my_edge = input_vertex_places[input[start_main_edge + threadIdx.x]];
        int end_my_edge = input_vertex_places[input[start_main_edge + threadIdx.x] + 1];
        for (int i = 0; i < k; i++)
        {
            res[i] = 0;
        }
        for (int i = start_my_edge; i < end_my_edge; i++)
        {
            while (actual_main_edge < len && sh_mem[OFFSET + actual_main_edge] < input[i])
            {
                actual_main_edge++;
            }
            if (actual_main_edge < len && sh_mem[OFFSET + actual_main_edge] == input[i])
            {
                sh_mem[offsets[threadIdx.x + 1] - MAX_DEGREE / 32 + actual_main_edge / 32] += (1 << (actual_main_edge % 32));
                res[2]++;
            }
        }
    }
    __syncthreads();
    if (threadIdx.x < end_main_edge - start_main_edge)
    {
    
        int stack[MAX_DEPTH];
        int stack_pointer = 2;
        int actual_stack_values[MAX_DEPTH][MAX_DEGREE / 32], pos = MAX_DEGREE / 32 - 1;
        for (int i = offsets[threadIdx.x + 1] - 1; i >= offsets[threadIdx.x]; i--)
        {
            actual_stack_values[1][pos] = sh_mem[i];
            pos--;
        }
        while (pos >= 0)
        {
            actual_stack_values[1][pos] = 0;
            pos--;
        }
        stack[0] = main_vertex; // does not matter
        stack[1] = threadIdx.x; // intersection of first vector and second vector
        stack[2] = threadIdx.x;
        while (stack_pointer >= 2)
        {
            stack[stack_pointer]++;
            if (stack[stack_pointer] >= len)
            {
                stack_pointer--;
            }
            else if (actual_stack_values[stack_pointer - 1][stack[stack_pointer] / 32] & (1 << (stack[stack_pointer] % 32))) // if 1 is on stack_ptr
            {

                int posx = MAX_DEGREE / 32 - 1 - not_used;
                int val = 0;
                int count = 0;

                for (int i = offsets[stack[stack_pointer] + 1] - 1 - not_used; i >= offsets[stack[stack_pointer]]; i--)
                {
                    int x = sh_mem[i] & actual_stack_values[stack_pointer - 1][posx];
                    actual_stack_values[stack_pointer][posx] = x;
                    count += __popc(x);
                    val |= x;
                    posx--;
                }
                res[stack_pointer + 1] += count;

                if (val != 0 && stack_pointer < k - 2)
                {
                    while (posx >= 0)
                    {
                        actual_stack_values[stack_pointer][posx] = 0;
                        posx--;
                    }

                    stack_pointer++;

                    stack[stack_pointer] = stack[stack_pointer - 1];
                }
            }
        }

        for (int i = 0; i < k; i++)
        {
            atomicAdd(&output[i], res[i] % MOD); // to optimize
        }
    }
}
std::unordered_map<int, int> renum_vertex;
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
    std::vector<int> edges_sorted, vertex_places, deg_count;

    while (input >> a >> b)
    {
        if (renum_vertex.find(a) == renum_vertex.end())
        {
            renum_vertex[a] = vertex_num;
            vertex_num++;
        }
        a = renum_vertex[a];
        if (renum_vertex.find(b) == renum_vertex.end())
        {
            renum_vertex[b] = vertex_num;
            vertex_num++;
        }
        b = renum_vertex[b];
        edges_list.emplace_back(a, b);
    }
    std::vector <int> vert_pos;
    for (int i = 0; i < vertex_num; i++)
    {
        edges.emplace_back();
        deg_count.push_back(0);
        vert_pos.push_back(0);
    }
    for (int i = 0; i < edges_list.size(); i++)
    {
        deg_count[edges_list[i].first]++;
        deg_count[edges_list[i].second]++;
    }
    std::vector <std::pair<int,int>> vert_sort{};
    for (int i = 0; i < vertex_num; i++)
    {
        vert_sort.emplace_back(deg_count[i], i);
    }
    std::sort(vert_sort.begin(), vert_sort.end());
    
    
    for(int i=0;i<vertex_num;i++)
    {
        vert_pos[vert_sort[i].second] = i; 
    }

    for (int i = 0; i < edges_list.size(); i++)
    {
        edges_list[i].first = vert_pos[edges_list[i].first];
        edges_list[i].second = vert_pos[edges_list[i].second];
        if (edges_list[i].second < edges_list[i].first)
        {
            std::swap(edges_list[i].first, edges_list[i].second);
        }
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
    vertex_places.push_back(counter);
    unsigned long long int output[MAX_DEGREE];
    for (int i = 0; i < k; i++)
    {
        output[i] = 0;
    }
    int *edges_gpu, *vertex_places_gpu;
    unsigned long long int *output_gpu;
    HANDLE_ERROR(cudaMalloc(&edges_gpu, sizeof(int) * edges_sorted.size()));
    HANDLE_ERROR(cudaMalloc(&vertex_places_gpu, sizeof(int) * (vertex_num+1)));
    HANDLE_ERROR(cudaMalloc(&output_gpu, sizeof(unsigned long long int) * k));
    HANDLE_ERROR(cudaMemcpy(edges_gpu, edges_sorted.data(), sizeof(int) * edges_sorted.size(), cudaMemcpyHostToDevice));
    HANDLE_ERROR(cudaMemcpy(vertex_places_gpu, vertex_places.data(), sizeof(int) * (vertex_num+1), cudaMemcpyHostToDevice));
    HANDLE_ERROR(cudaMemcpy(output_gpu, output, sizeof(unsigned long long int) * k, cudaMemcpyHostToDevice));
    // HANDLE_ERROR(cudaDeviceSynchronize());
    int constants_to_gpu[MAX_DEGREE + 1], mem_offset = MAX_DEGREE;
    constants_to_gpu[0] = 0;
    for (int i = 1; i <= MAX_DEGREE; i++)
    {
        constants_to_gpu[i] = constants_to_gpu[i - 1] + (mem_offset - 1) / 32 + 1;
        mem_offset--;
    }
    cudaMemcpyToSymbol(offsets, constants_to_gpu, sizeof(int) * (MAX_DEGREE + 1));
    sort_edges<<<vertex_num, MAX_DEGREE>>>(edges_gpu, vertex_places_gpu); // optimize to infinity

    // HANDLE_ERROR(cudaDeviceSynchronize());

    //int maxbytes = 65536;
    int maxbytes = 98304;// to restore
    cudaFuncSetAttribute(kcliques, cudaFuncAttributeMaxDynamicSharedMemorySize, maxbytes);
    kcliques<<<vertex_num, MAX_DEGREE, maxbytes>>>(k, edges_gpu, vertex_places_gpu, output_gpu);
    HANDLE_ERROR(cudaMemcpy(output, output_gpu, sizeof(unsigned long long int) * k, cudaMemcpyDeviceToHost));
    // HANDLE_ERROR(cudaDeviceSynchronize());
    output[0] = vertex_num;
    output[1] = edges_sorted.size();
    std::ofstream output_stream;
    output_stream.open(std::string(argv[3]));
    for (int i = 0; i < k; i++)
    {
        output_stream << output[i] % MOD << " ";
    }
}