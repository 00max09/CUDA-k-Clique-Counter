#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#define MAX_DEGREE 1024
#define MAX_DEPTH 12

__global__ sort_edges(int* input, int*input_vertex_places){
    __shared__ int to_sort[MAX_DEGREE];
    __shared__ int to_sort2[MAX_DEGREE];
    
    int main_vertex = blockIdx.x;
    int start_main_edge = input_vertex_places[main_vertex];
    int end_main_edge = input_vertex_places[main_vertex + 1];
    int len = end_main_edge-start_main_edge;
    to_sort[threadIdx.x] = input[start_main_edge + threadIdx.x];
    __syncthreads();
    for(int i = 2; i<=2*end_main_edge-start_main_edge;i*=2)
    {
        if(threadIdx.x * i <len)
        {
           int start = threadIdx * i;
           int end = min(start + i, len);
           int mid = min(start + i/2,len);
           int mid2 = start+i/2;
           int count = threadIdx * i;
           while(mid2 < mid || start<mid){
                if(start<mid && to_sort[start]<to_sort[mid2]){
                    to_sort2[count] = to_sort[start];
                    start++;
                }else{
                    to_sort2[count] = to_sort[end];
                    end++;
                }
                count++;
           }
        }
        __syncthreads();
        to_sort[threadIdx.x] = to_sort2[threadIdx.x];
        __syncthreads();
    }
    input[start_main_edge + threadIdx.x] = to_sort[threadIdx.x];
}

__global__ kcliques(int k, int *input, int *input_vertex_places, int *output) {
    int main_vertex = blockIdx.x;
    __shared__ int incidence_matrix[MAX_DEGREE][MAX_DEGREE / 32];
    __shared__ int main_input[MAX_DEGREE]; //copy there
    int start_main_edge = input_vertex_places[main_vertex];
    int end_main_edge = input_vertex_places[main_vertex + 1];
    int actual_main_edge = 0;
    int start_my_edge = input_vertex_places[start_main_edge + threadIdx.x];   // if x to big don't do
    int end_my_edge = input_vertex_places[start_main_edge + threadIdx.x + 1]; // if x to big don't do
    for (int i = start_my_edge; i < end_my_edge; i++) {
        while(main_input[actual_main_edge] < input[i]){
            actual_main_edge++;
            if(actual_main_edge >= end_main_edge- start_main_edge)
                break;            
        }
        if(main_input[actual_main_edge] == input[i])
        {
            incidence_matrix[main_vertex][actual_main_edge/32] |= 1<<(actual_main_edge%32); //to check move
        }
    }
    int stack[MAX_DEPTH];
    int res[MAX_DEPTH];
    for (int i=0;i<k;i++)
    {
        res[i]=0;
    }
    int stack_pointer = 2;
    __shared__ int actual_stack_values[MAX_DEGREE][MAX_DEGREE/32];
    for(int i=0;i<MAX_DEGREE/32;i++)
    {
        actual_stack_values[threadIdx.x][i] = incidence_matrix[threadIdx.x][i]; 
    }

    stack[0] = main_vertex;
    stack[1] = threadIdx.x;
    stack[2] = threadIdx.x+1;
    while(stack_pointer > 2)
    {
        stack[stack_pointer]++;
        if(stack[stack_pointer] > MAX_DEGREE)
        {
            stack_pointer--;
            for(int i=0;i<stack_pointer;i++)
            {
                for(int z=stack[stack_pointer]/32;i<MAX_DEGREE/32;i++)
                {
                    if(i==0){
                        actual_stack_values[threadIdx.x][z] = incidence_matrix[stack[i]][z];
                    }else
                    {
                        actual_stack_values[threadIdx.x][z] &= incidence_matrix[stack[i]][z];
                    }
                }
            }
        }
        if(actual_stack_values[threadIdx.x][stack[stack_pointer]/32] & 1<<stack[stack_pointer]%32) // if 1 is on stack_ptr
        {
            res[stack_pointer]++;
            if(stack_pointer<k-1)
            {
                for(int i=stack[stack_pointer]/32;i<MAX_DEGREE/32;i++)
                {
                    actual_stack_values[threadIdx.x][i] &= incidence_matrix[stack[stack_pointer]][i];
                }
                stack_pointer++;
                stack[stack_pointer] = stack[stack_pointer-1]+1;
            }
        }
    }
    for(int i=0;i<k;i++)
    {
        atomicAdd(&output[i], res[i]); //to optimize
    } 
}

int main(int argc, char *argv[]) {
    if (argc != 4)
    {
        std::cout << "Usage ./kcliques <graph input file> <k value> <output file>";
        return 1;
    }
    int k = stoi(std::string(argv[2]));
    std::ifstream input;
    input.open(std::string(argv[1]));
    unsigned int a, b, vertex_num;
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
    for(int i=0;i<vertex_num;i++)
    {
        edges.emplace_back();
    }
    for(int i=0;i<edges_list.size();i++)
    {
        edges[edges_list[i].first].push_back(edges_list[i].second);
    }
    int counter = 0;
    for(int i=0;i<edges.size();i++)
    {
        vertex_places.push_back(counter);
        for(int z=0;z<edges[i].size();z++)
        {
            edges_sorted.push_back(edges[i][z]);
        }
        counter+=edges[i].size();
    }
    int output[k];
    for(int i=0;i<k;i++)
    {
        output[i] = 0;
    }
    vertex_places.push_back(counter);
    int * edges_gpu, vertex_places_gpu, output_gpu;
    cudaMalloc((void **)&edges_gpu, sizeof(int) * edges_sorted.size());
    cudaMalloc((void **)&vertex_places_gpu, sizeof(int) * (vertex_num+1));
    cudaMalloc((void **)&output_gpu, sizeof(int) * k);
    cudaMemcpy(edges_gpu, edges_sorted.data(), sizeof(int) * edges_sorted.size(), cudaMemcpyHostToDevice);
    cudaMemcpy(vertex_places_gpu, vertex_places.data(), sizeof(int) * (vertex_num+1), cudaMemcpyHostToDevice);
    cudaMemcpy(output_gpu, output, sizeof(int) * k, cudaMemcpyHostToDevice);
    
    sort_edges<<<vertex_num, MAX_DEGREE>>>(edges_gpu, vertex_places_gpu); //optimize to infinity
    kcliques<<<vertex_num, MAX_DEGREE>>>(edges_gpu, vertex_places_gpu, k, output_gpu);
    cudaMemcpy(output, output_gpu, sizeof(int) * k, cudaMemcpyDeviceToHost);
    

}