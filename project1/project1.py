import sys

import networkx as nx
import numpy as np
import csv
from scipy.special import gammaln
import random

n = 0

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
    
    # Set global number of variables
    n = len(header)

    # For each variable, store it's maximum instantiation value over all data points
    vars = []
    num_data_points = len(data)
    for var_idx in range(n):
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


def get_statistics(vars, graph, data):
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
    
    # Now get the prior, assume uniform (all entries = 1)
    prior = [ones(q[var], r[var]) for var range(n)]
    
    return M, prior


# Assuming uniform prior so we drop log P(G)
def bayesian_score_component(M, prior):
    p = np.sum(gammaln(np.sum(prior, axis=2)))
    p -= np.sum(gammaln(np.sum(prior, axis=2) + np.sum(M, axis=2)))
    p += np.sum(gammaln(prior + M))
    p -= np.sum(gammaln(prior))
    return p

def get_score(vars, G, data):
    M, prior = get_statistics_and_prior(vars, G, data)

    # Sum the bayesian score components of each variable, using that variable's
    # prior & that variable's associated counts
    return sum(bayesian_score_component(M[i], prior[i]) for i in range(n))


def rand_graph_neighbor(G):
    nodes = list(G.nodes)
    node_1 = random.choice(nodes)
    nodes_without_node_1 = [node for node in nodes if node != node_1]
    node_2 = random.choice(nodes_without_node_1)
    
    # Create a copy so we do not mutate original graph 
    G_cpy = G.copy()

    # Preform some graph mutation
    if G.has_edge(node_1, node_2):
        # Remove edge half the time
        if (random.randint(0, 1)):
            G_cpy.remove_edge(node_1, node_2)
        # Reverse edge half the time
        else: 
            G_cpy.remove_edge(node_1, node_2)
            G_cpy.add_edge(node_2, node_1)
    elif G.has_edge(node_2, node_1):
        # Remove edge half the time
        if (random.randint(0, 1)):
            G_cpy.remove_edge(node_2, node_1)
        # Reverse edge half the time
        else: 
            G_cpy.remove_edge(node_2, node_1)
            G_cpy.add_edge(node_1, node_2)
    else:
        G_cpy.add_edge(node_1, node_2)

    return G_cpy



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

    # TO DO - check what other output we might need to create


if __name__ == '__main__':
    main()
