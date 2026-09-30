import sys

import networkx as nx
import numpy as np
import csv

class Variable:
    def __init__(self, max_val):
        self.max_val = max_val

def write_gph(dag, filename):
    with open(filename, 'w') as f:
        for edge in dag.edges():
            f.write("{}, {}\n".format(edge[0], edge[1]))


def get_file_info(infile):
    # Variables
    header = []
    # (data points x number of variables)
    data = []

    with open(infile, 'r') as f:
        reader = csv.reader(f)
        header = next(reader)

        for row in reader:
            data.append([int(x) for x in row])
    
    # For each variable, store it's maximum instantiation value over all data points
    vars = []
    num_data_points = len(data)
    num_vars = len(data[0])
    for var_idx in range(num_vars):
        col = [data[i][var_idx] for i in range(num_data_points)]
        max_val = max(col)
        vars.append(Variable(max_val))

    # Create an inital graph, start with no edges
    G = nx.DiGraph()
    G.add_nodes_from(header)

    return inital_G, header, data, vars


def sub2ind(parents_max_vals, parents_values):
    k = np.concatenate(([1], np.cumprod(parents_max_vals[:-1])))
    return np.dot(k, parents_values)


# Assume a uniform prior over all graphs
def get_statistics(vars, graph, data):
    # Number of variables
    n = len(data[0])

    # Max value of each variable
    r = [vars[var].max_val for var in range(n)]

    # Number of possible parental instantiations for each variable
    q = [np.prod([r[parent] for parent in graph.predecessors(var)]) for var in range(n)]

    # Matrix with dimensions q_i x r_i per variable i
    M = [np.zeros((q[var], r[var])) for var in range(n)]

    # For each data point & each variable within that data point
    for data_point in data:
        for var in range(n):
            # For 0-indexing
            val_var = data_point[var] - 1
            parents = list(graph.predecessors(var))

            parent_combo = 0
            if len(parents) > 0:
                # Get linear index for combo of insantated parents
                parent_combo = sub2ind([r[parent] for parent in parents], [data_point[parent]-1 for parent in parents])

            # Update counts
            M[var][parent_combo, val_var] += 1.0
    
    return M

def get_score(vars, G, data):
    # TO DO

def rand_graph_neighbor(G):
    

def fit(G, vars, data, max_iters):
    y = get_score(vars, G, data)

    # Repeat for max_iters
    for k in max_iters:
        # Get a random graph neighbor of G
        rand_G_neighbor = rand_graph_neighbor(G)

        # If random graph has a cycle, it is invalid. Assign a score of -infinity
        new_y = -np.inf
        if (not nx.is_directed_acyclic_graph(rand_G_neighbor)):
            new_y = get_score(vars, rand_G_neighbor, data)
        
        # If random graph yields a better score than G, replace G with it
        if new_y > y:
            y = new_y
            G = rand_G_neighbor
            
    # Return best graph found after max_iter iterations
    return G


def compute(infile, outfile, max_iters):
    # Get info about the infile
    inital_G, header, data, vars = get_file_info(infile)

    # Run locally, directed graph search. Oportunistically move to random graph neightbor
    # if it's Baysian Score is greater. 
    G = fit(inital_G, vars, data, max_iters)
       
    return G


def main():
    if len(sys.argv) != 3:
        raise Exception("usage: python project1.py <infile>.csv <outfile>.gph")

    inputfilename = sys.argv[1]
    outputfilename = sys.argv[2]
    G = compute(inputfilename, outputfilename, 1000)

    # Write edges of dag to output file
    write_gph(G, outputfilename + '.gph')


if __name__ == '__main__':
    main()
