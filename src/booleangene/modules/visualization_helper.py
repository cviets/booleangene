from .BooleanNetwork import BooleanNetwork
import networkx as nx
import matplotlib.pyplot as plt
from .io_utils import read_expressions_from_txt, infer_external_from_expressions
from .regenerate_network import regenerate_network

def make_graph(bnetwork: BooleanNetwork, labels:bool=False, seed=None, pos=None):

    node_size = 3000
    G = nx.DiGraph()
    visual = []
    for node, children in bnetwork.get_out_degs().items():
        for child in children:
            visual.append([node, child])

    G.add_edges_from(visual)
    fig, _ = plt.subplots(1,1,figsize=(15,15))
    if pos is None:
        pos = nx.forceatlas2_layout(G)
    genes = [elt for elt in bnetwork.external_values] + [elt for elt in bnetwork.state]
    nodes = nx.draw_networkx_nodes(G, pos, node_color='white', nodelist=genes, linewidths=3, node_size=node_size)
    nodes.set_edgecolor('black')
    _ = nx.draw_networkx_edges(G, pos, arrows=True, width=3, arrowsize=30, node_size=node_size)

    if labels:
        _ = nx.draw_networkx_labels(G, pos, verticalalignment='center', font_weight='bold', font_family='Arial', font_size=24)

    return fig, pos

def visualize(expressions_file: str, labels:bool=False, seed=None):
    expressions = read_expressions_from_txt(expressions_file)
    external = infer_external_from_expressions(expressions)
    net = BooleanNetwork(expressions, external)
    new = regenerate_network(net, True, True, True, 1)
    graph_original, pos = make_graph(net, labels, seed, None)
    graph_scrambled, _ = make_graph(new, labels, seed, pos)

    return graph_original, graph_scrambled