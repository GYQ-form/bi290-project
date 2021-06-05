import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import sys

def main():
    if len(sys.argv) == 1:
        print("please input the name of the file")
        return 1
    fp = open(str(sys.argv[1]), 'r')
    fp = open("slim_extract", 'r')
    node = []
    edge = []
    lines = fp.readlines()
    i = 0
    for line in lines:
        line = line.replace("\n", "")
        line = str(line)
        if i % 2 == 0:
            node.append(line)
            i = i + 1
            continue
        if i % 2 == 1:
            line = line.replace("\t", ",")
            line = line.split(',')
            num = len(line)
            for j in range(num - 1):
                tmp_num = line[j]
                tmp_edge = (node[int(i / 2)], tmp_num)
                edge.append(tmp_edge)
            i = i + 1
            continue
    nx.DiGraph()
    G = nx.DiGraph(name='graph')
    node = set(node)
    edge = set(edge)
    G.add_nodes_from(node)
    G.add_edges_from(edge, color='red', length=5000)
    # print("nodes:", G.nodes())
    # print("edges:", G.edges())
    # print("number of edges:", G.number_of_edges())
    # print(node)
    # print(edge)
    options = {
        "node_color": "#A0CBE2",
        "node_size": 300,
        "font_size": 4,
        "edge_color": "g",
        "width": 0.5,
        # "edge_cmap": plt.cm.Blues,
        "edge_length": 500,
        "with_labels": True,
    }
    plt.figure(3, figsize=(50, 50))
    nx.draw_networkx(G, **options)
    plt.savefig("photo.png", dpi=200)
    # print(len(node))
    # pos=nx.networkx.random_layout(G)
    print("Image generated")
    return 0
main()

