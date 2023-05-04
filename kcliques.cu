#include <iostream>
#include <fstream>
#include <vector>
#define MAX_DEGREE 1024
#define MAX_DEPTH 12
__global__ kcliques(int k, int *input, int *input_vertex_places) {
    int main_vertex = blockIdx.x;
    __shared__ int incidence_matrix[MAX_DEGREE][MAX_DEGREE / 32];
    __shared__ int 
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
    int stack_pointer = 2;
    stack[0] = main_vertex;
    stack[1] = threadIdx.x;
    stack[2] = 0;
}

int main(int argc, char *argv[]) {
    if (argc != 4)
    {
        std::cout << "Usage ./kcliques <graph input file> <k value> <output file>";
        return 1;
    }
    std::ifstream input;
    input.open(std::string(argv[1]));
    unsigned int a, b;
    std::vector<std::pair<int, int>> edges;
    while (input >> a >> b)
    {
        if (a > b)
        {
            std::swap(a, b);
        }
        edges.emplace_back(a, b);
    }
}