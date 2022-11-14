
#Finding the Highest K-Core Subgraph in a Protein-Protein Interaction Network

import networkx as nx #Python software package for study of complex networks
from sortedcontainers import SortedDict #Sorted dict keys are maintained in sorted order.

ppi_file = 'ppi_file.txt'  #Path to a network file

def create_network(ppi_file) :
    with open(ppi_file, "r") as ppi: # Open the protein- protein interaction file
        graph = nx.Graph() # Creation of a empty network graph
        for pro_int in ppi: # Each protein interaction in the network
            node = pro_int.rstrip("\n") #Strip edge
            node1= node.split("\t") #Split edge
            graph.add_edge(node1[0], node1[1]) # Add nodes to the graph
        #Summary of the graph
        print(nx.info(graph))
    return graph

def kcores(ppi_file):
    k_cores = {}  #Dictionary based adjacency list representation for the PPI network to hold proteins and the highest k-core they belong to
    highest_kcore =0    #Follow the hightest k-core
    network = create_network(ppi_file)         #Create a network via first function
    protein_cores = nx.core_number(network)    #Max k-cores from the network for each protein with networkx package.

    #Group proteins according to the highest number of k-cores a protein has within the network
    for protein, k_core in protein_cores.items():
        if k_core > highest_kcore:  #Find highest k-core by comparing it k-core of each protein.
            highest_kcore = k_core

        #Sort the proteins according to their highest k_core in a dictionary
        if k_core in k_cores:
            k_cores[k_core].append(protein)
        else:
            k_cores[k_core]=[protein]
    return k_cores

k_cores =kcores(ppi_file) #all possible cores
sor_kcores= SortedDict(k_cores) #sort the dictionary with respect to numerical order of keys
#print k-cores from the dictionary and bring the number of proteins specific to that k-core
for key in sor_kcores:
    # with open("output.txt", "a") as k_core  #used for printing output.
    print("For k = " + str(key) + " there are "+ str(len(sor_kcores[key])) + " proteins. \n")
    # add (,file= k_core) while printing the output
