import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx

# remove all invariant columns from a dataframe
def remove_constant_columns(df):
    return df.loc[:, df.nunique() > 1]

# keep only invariant columns of a dataframe
def keep_constant_columns(df):
    return df.loc[:, df.nunique() == 1]

# create a correlation network from a presence / absence matrix
def create_network(pa_matrix, nodes, thresh=0, pos_only=False, neg_only=False):
    nodes = [x for x in nodes if x in pa_matrix.columns]
    nodes_variable = pa_matrix[nodes]
    corr_mat = nodes_variable.corr().values
    G = nx.Graph()
    G.add_nodes_from(nodes)
    for  i in range(len(nodes)):
        for j in range(i+1, len(nodes)):
            weight = corr_mat[i][j] # extract edge weights from the correlation matrix
            if weight > thresh and not neg_only:
                G.add_edge(nodes[i], nodes[j], weight=weight)
            if weight < -1*thresh and not pos_only:
                G.add_edge(nodes[i], nodes[j], weight=weight)
    return G

# handle negative weights in the correlation network
def handle_neg_only(G):
    G_pos = G.copy()
    if np.array([d['weight'] < 0 for u,v,d in G_pos.edges(data=True)]).all():
        for u, v, d in G_pos.edges(data=True):
            d['weight'] *= -1
    return G_pos

# function for drawing correlation graphs from a set of parameters
def draw_correlation_graph(G, ax, color_palette_dict, sizes_dict, params={}, title="", circled_nodes=None, annotate=False):
    default_params = {'random_seed': 42,
                      'iterations': 50,
                      'edge_alpha':0.5,
                      'edge_width': 0.5,
                      'circled_node_alpha': 0.5,
                      'circled_node_width': 5,
                      #'pos_color':"#b65353",
                      'pos_color':"#808080",
                      'neg_color':"#4b7cbc",
                      'title_font_size':50
                      }
    for k, v in default_params.items():
        if k not in params.keys():
            params[k] = v
    
    # set parameters        
    rs = params["random_seed"]
    it = params["iterations"]
    ea = params["edge_alpha"]
    ew = params["edge_width"]
    cna = params["circled_node_alpha"]
    cnw = params["circled_node_width"]
    pc = params["pos_color"]
    nc = params["neg_color"]
    fs = params["title_font_size"]
    
    G_lcc = G.copy()

    # set node parameters
    node_sizes = [50 * sizes_dict[n] for n in G_lcc.nodes]
    edge_colors = [pc if weight > 0 else nc for u, v, weight in G_lcc.edges(data='weight')]
    node_colors = [color_palette_dict[n] for n in G_lcc.nodes]
    
    G_pos = handle_neg_only(G_lcc)
    pos = nx.spring_layout(G_pos,
                           seed=rs,
                           iterations=it
                           )
    
    # label nodes where the gene groups have sufficiently short names
    labels = {node: node for node in G_pos.nodes() if len(node) <= 5}
    
    density = nx.density(G_pos)
    if density!=0:
        edge_widths = [abs(d['weight']) * abs(d['weight']) * ew / density for u, v, d in G_pos.edges(data=True)]
    else:
        edge_widths = [abs(d['weight']) * abs(d['weight']) * ew for u, v, d in G_pos.edges(data=True)]

    nx.draw_networkx_edges(G_pos, pos,
                           width=edge_widths,
                           edge_color=edge_colors,
                           alpha=ea,
                           ax=ax
                           )
    
    # highlight a specified set of nodes
    if circled_nodes:
        linewidths = [cnw if n in circled_nodes else 0 for n in G_pos.nodes]
        alpha = [1 if n in circled_nodes else cna for n in G_pos.nodes]
        circled_nodes = [n for n in circled_nodes if n in G_pos.nodes]
        non_circled_nodes = [n for n in G_pos.nodes if n not in circled_nodes]

        nx.draw_networkx_nodes(G_pos, pos,
                               nodelist=non_circled_nodes,
                               node_size=[node_sizes[list(G_pos.nodes).index(n)] for n in non_circled_nodes],
                               node_color=[node_colors[list(G_pos.nodes).index(n)] for n in non_circled_nodes],
                               alpha=[alpha[list(G_pos.nodes).index(n)] for n in non_circled_nodes],
                               linewidths=[linewidths[list(G_pos.nodes).index(n)] for n in non_circled_nodes],
                               edgecolors="black",
                               ax=ax
                               )

        nx.draw_networkx_nodes(G_pos, pos,
                               nodelist=circled_nodes,
                               node_size=[node_sizes[list(G_pos.nodes).index(n)] for n in circled_nodes],
                               node_color=[node_colors[list(G_pos.nodes).index(n)] for n in circled_nodes],
                               alpha=[alpha[list(G_pos.nodes).index(n)] for n in circled_nodes],
                               linewidths=[linewidths[list(G_pos.nodes).index(n)] for n in circled_nodes],
                               edgecolors="black",
                               ax=ax
                               )
    else:
        linewidths=0.5
        alpha=1
        nx.draw_networkx_nodes(G_pos, pos,
                               node_size=node_sizes,
                               node_color=node_colors,
                               alpha=alpha,
                               linewidths=linewidths,
                               edgecolors="black",
                               ax=ax
                               )
    if annotate:
        nx.draw_networkx_labels(G_pos, pos,
                                labels,
                                font_size=32,
                                font_color='black',
                                ax=ax)

    for spine in ax.spines:
        ax.spines[spine].set_visible(False)
    ax.set_title(f"{title}",
                 fontsize=fs)

# decompose the graph into k independent sets
def kpartite_decomposition(G):
    coloring = nx.coloring.greedy_color(G, strategy='largest_first')
    coloring_reindexed = {key: value + 1 for key, value in coloring.items()}
    kpartite_sets = {}
    # create a dictionary of k-partite groups to gene group content
    for gene, group_number in coloring_reindexed.items():
        if group_number not in kpartite_sets:
            kpartite_sets[group_number] = []
        kpartite_sets[group_number].append(gene)
    return kpartite_sets

# compute the number of k-partite sets
def k_partite(G):
    try: # try to color the nodes with at most k colors
        coloring = nx.coloring.greedy_color(G, strategy='largest_first')
    except nx.NetworkXError:
        return False
    # check the number of colors used, i.e. the number of k-partite sets
    return len(set(coloring.values()))

# examine various metrics of a given network
def graph_theory_analysis(G):
    print(f"Number of nodes: {nx.number_of_nodes(G)}")
    print(f"Number of edges: {nx.number_of_edges(G)}")
    print(f"Number of connected components: {len(sorted(nx.connected_components(G), key=len, reverse=True))}")
    print(f"Density: {nx.density(G)}")
    print(f"Transitivity: {nx.transitivity(G)}")
    #print(f"Average degree connectivity: {nx.average_degree_connectivity(G)}")
    print(f"Average clustering: {nx.average_clustering(G)}")
    print(f"Maximum clique size: {nx.approximation.large_clique_size(G)}")
    #print(f"Maximum independent set: {nx.approximation.maximum_independent_set(G)}")
    #print(f"Dominating set: {nx.dominating_set(G)}")
    #print(f"Is bipartite? {nx.bipartite.is_bipartite(G)}")
    print(f"k-partite {k_partite(G)}")
    if nx.is_connected(G):
        print(f"Average shortest path length: {nx.average_shortest_path_length(G)}")
        print(f"Diameter: {nx.diameter(G)}")



