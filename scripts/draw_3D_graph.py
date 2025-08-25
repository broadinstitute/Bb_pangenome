#!/usr/bin/env python3

import pandas as pd
import numpy as np
import igraph as ig
import chart_studio.plotly as ply
import plotly.graph_objects as go
import plotly.express as px
import plotly.io as pio
from scipy.spatial import ConvexHull
import random
import os
import pathlib
from tqdm import tqdm
import kaleido
import json
import argparse

def edges_from_matrix(matrix, output, idmap, min_contig_length=750, min_weight=0.0):
    """
    Convert alignment matrix to edge list with weights calculated as Jaccard index.
    """
    fixed_map = idmap[['assembly_id', 'contig_id', 'plasmid_name']]
    fixed_map['name'] = fixed_map.apply(lambda x: x['assembly_id'] + '__' + x['contig_id'], axis=1)
    fixed_map = fixed_map[['name', 'plasmid_name']]
    pmap = {}
    for _, row in fixed_map.iterrows():
        pmap[row['name']] = row['plasmid_name']

    df = pd.read_csv(matrix, sep='\t', index_col=0, header=0)

    source = []
    target = []
    weight = []
    length_i = []
    length_j = []
    align_length = []

    for i in tqdm(range(len(df.index)), desc='Processing matrix rows'):
        for j in range(i+1, len(df.columns)):
            if df.iloc[i,j] > 0:
                if df.iloc[i,i] >= min_contig_length and df.iloc[j,j] >= min_contig_length:
                    score = ((2 * df.iloc[i,j]) / (df.iloc[i,i] + df.iloc[j,j]))
                    
                    if score >= min_weight:
                        source.append(df.index[i])
                        target.append(df.columns[j])
                        weight.append(score)
                        length_i.append(df.iloc[i,i])
                        length_j.append(df.iloc[j,j])
                        align_length.append(df.iloc[i,j])

    edge_df = pd.DataFrame({
        'Source': source,
        'Target': target,
        'weight': weight,
        'source_length': length_i,
        'target_length': length_j,
        'alignment_length': align_length,
        'interaction': ['interacts'] * len(source)
    })
    edge_df['source_plasmid'] = edge_df['Source'].apply(lambda x: pmap.get(x, 'Unclassified'))
    edge_df['target_plasmid'] = edge_df['Target'].apply(lambda x: pmap.get(x, 'Unclassified'))
    print(f"\nEdge Statistics:")
    print(f"Total edges: {len(edge_df)}")
    print(f"Unique source nodes: {len(edge_df['Source'].unique())}")
    print(f"Unique target nodes: {len(edge_df['Target'].unique())}")
    print(f"Weight range: {edge_df['weight'].min():.3f} - {edge_df['weight'].max():.3f}")
    
    edge_df.to_csv(output, index=False)
    print(f'Saved edges to {output}')
    return edge_df

def create_igraph(edge_df):
    nodes = list(set(edge_df['Source'].to_list() + edge_df['Target'].to_list()))
    nodes = [str(node) for node in nodes]
    G = ig.Graph()
    G.add_vertices(nodes)
    edges = [(str(row['Source']), str(row['Target'])) for _, row in edge_df[['Source', 'Target']].iterrows()]
    G.add_edges(edges)
    G.es['weight'] = edge_df['weight'].tolist()
    return G, nodes

def align_labels_and_groups(idmap, vertex_order):
    node_dict = idmap.apply(lambda x: f"{x['assembly_id']}__{x['contig_id'].split(' ')[0]}", axis=1)
    node_dict = dict(zip(node_dict, idmap['plasmid_name']))
    
    aligned_labels = []
    aligned_groups = []

    print("Sample of idmap:")
    print(idmap.head())
    print("\nSample of vertex_order:")
    print(vertex_order[:10])

    plasmid_mappings = {
        'cp32-2': 'cp32-7',
        'cp32-9-4': 'cp32-9',
        'cp32-1+5': 'cp32-1',
        'cp32-3+10': 'cp32-3',
        'cp32-5+1': 'cp32-5',
        'cp32-5-1': 'cp32-5'
    }

    special_cases = {
        'URI88H_contig000014': 'lp21-cp9'
    }

    for vertex in vertex_order:
        if vertex in special_cases:
            aligned_groups.append(special_cases[vertex])
            aligned_labels.append(vertex)
        elif vertex in node_dict:
            group = node_dict[vertex]
            group = plasmid_mappings.get(group, group)
            aligned_groups.append(group)
            aligned_labels.append(vertex)
        else:
            print(f'Fragment found: {vertex} (assigned as unknown)')
            aligned_labels.append(vertex)
            aligned_groups.append('unknown')

    print(f"\nTotal vertices: {len(vertex_order)}")
    print(f"Aligned labels: {len(aligned_labels)}")
    print(f"Aligned groups: {len(aligned_groups)}")
    print(f"Unique groups: {set(aligned_groups)}")

    return aligned_labels, aligned_groups

def create_color_mapping(groups):
    color_map = {
        'cp26': '#d60000', 'lp54': '#018700', 'lp17': '#b500ff', 'lp36': '#05acc6', 
        'lp28-3': '#97ff00', 'lp25': '#ffa52f', 'cp32-7': '#ff8ec8', 'lp28-4': '#79525e', 
        'lp38': '#00fdcf', 'cp32-9': '#afa5ff', 'cp32-6': '#93ac83', 'cp32-11': '#9a6900', 
        'lp28-1': '#366962', 'cp32-3': '#d3008c', 'cp32-1+5': '#fdf490', 'cp32-4': '#c86e66', 
        'cp32-12': '#9ee2ff', 'lp56': '#00c846', 'lp28-6': '#a877ac', 'cp32-8': '#b8ba01', 
        'lp28-2': '#f4bfb1', 'cp32-1': '#ff28fd', 'cp32-2': '#f2cdff', 'lp28-5': '#009e7c', 
        'cp32-5': '#ff6200', 'cp9': '#56642a', 'cp32-10': '#953f1f', 'lp28-7': '#90318e', 
        'lp28-8': '#ff3464', 'cp32-13': '#a0e491', 'lp5': '#8c9ab1', 'lp21': '#829026', 
        'cp9-3': '#ae083f', 'lp21-cp9': '#77c6ba', 'cp32-9-4': '#bc9157', 'lp28-9': '#e48eff', 
        'cp32-3+10': '#72b8ff', 'lp28-11': '#c6a5c1', 'lp32-3': '#ff9070', 'chromosome': '#bceddb',
        'Unclassified': '#d3c37c', "none": "hsla(0, 0.00%, 0.00%, 0.00)"
    }

    unique_groups = set(groups)
    missing = unique_groups - set(color_map.keys())
    if missing:
        print(f"\nWarning: Missing colors for groups: {missing}")
        for group in missing:
            color_map[group] = "#cccccc"
    
    return color_map

def make_3d_plot(G, labels, groups, output_file, layout, include_edges=True):
    """
    Create a 3D visualization of the graph.
    """
    layout = G.layout(layout, dim=3)
    Xn = [layout[k][0] for k in range(len(G.vs))]
    Yn = [layout[k][1] for k in range(len(G.vs))]
    Zn = [layout[k][2] for k in range(len(G.vs))]
    
    color_map = create_color_mapping(groups)
    unique_groups = list(set(groups))
    
    data = []
    
    # Add edges if requested
    if include_edges:
        Xe, Ye, Ze = [], [], []
        for e in G.es:
            Xe += [layout[e.source][0], layout[e.target][0], None]
            Ye += [layout[e.source][1], layout[e.target][1], None]
            Ze += [layout[e.source][2], layout[e.target][2], None]
        
        edge_trace = go.Scatter3d(
            x=Xe, y=Ye, z=Ze,
            mode='lines',
            name='Edges',
            line=dict(color='rgb(125,125,125)', width=0.35),
            hoverinfo='none',
            showlegend=False
        )
        data.append(edge_trace)

    # Add nodes for each group
    for group in unique_groups:
        group_indices = [i for i, g in enumerate(groups) if g == group]
        data.append(
            go.Scatter3d(
                x=[Xn[i] for i in group_indices],
                y=[Yn[i] for i in group_indices],
                z=[Zn[i] for i in group_indices],
                mode='markers',
                name=group,
                marker=dict(
                    symbol='circle',
                    size=6,
                    color=color_map[group],
                    line=dict(color='rgb(38, 38, 38)', width=0.35),
                ),
                text=[f'Label: {labels[i]}<br>Group: {group}' for i in group_indices],
                hoverinfo='text'
            )
        )

    # Create slider for edge opacity
    steps = []
    for step in np.arange(0, 1.1, 0.1):
        step = round(step, 2)
        steps.append(
            dict(
                method="update",
                args=[{"opacity": [step, *[1]*len(unique_groups)]}],  # First trace is edges, rest are nodes
                label=str(step)
            )
        )
    sliders = [dict(
        active=5,
        currentvalue={"prefix": "Edge Opacity: "},
        pad={"t": 5, "b": 10},
        steps=steps
    )]

    layout3d = go.Layout(
        title="Sequence Homology between contigs (3D)",
        scene=dict(
            xaxis=dict(title=''),
            yaxis=dict(title=''),
            zaxis=dict(title=''),
        ),
        margin=dict(r=0, l=0, b=0, t=100),
        hovermode='closest',
        legend=dict(
            itemsizing='constant',
            title_text='Plasmids',
            bgcolor='rgba(255,255,255,0.5)',
            bordercolor='rgba(0,0,0,0)',
            borderwidth=2
        ),
        sliders=sliders
    )

    fig = go.Figure(data=data, layout=layout3d)
    fig.update_layout(
        scene_aspectmode='data',
        autosize=True,
        uirevision=True
    )

    # Save as PNG
    img = pio.to_image(fig, 'png', width=5000, height=5000, scale=1)
    png_file = f'{pathlib.Path(output_file).parent}/{pathlib.Path(output_file).stem}.png'
    with open(png_file, 'wb') as outfile:
        outfile.write(img)
    print(f"PNG written to {png_file}")

    # Save as HTML if edges are included
    if include_edges:
        config = {
            'responsive': True,
            'scrollZoom': True,
        }
        pio.write_html(fig, file=output_file, full_html=False, include_plotlyjs='cdn', config=config)
        print(f"HTML written to {output_file}")

def igraph_to_coords_file(G, aligned_labels, groups):
    """
    Convert igraph object to coordinates file format.
    """
    color_map = create_color_mapping(groups)
    color_values = [color_map[group] for group in groups]
    layout = G.layout('kk', dim=3)
    
    Xn = [layout[k][0] for k in range(len(G.vs))]
    Yn = [layout[k][1] for k in range(len(G.vs))]
    Zn = [layout[k][2] for k in range(len(G.vs))]
    
    nodes = {
        i: {
            "location": [Xn[i], Yn[i], Zn[i]],
            "name": aligned_labels[i],
            "color": color_values[i],
            "plasmid": groups[i],
        } for i in range(G.vcount())
    }
    
    edges = [{"source": e.source, "target": e.target} for e in G.es]
    return {"nodes": nodes, "edges": edges}

def main():
    parser = argparse.ArgumentParser(description='Create 3D network visualization from alignment matrix')
    parser.add_argument('--matrix', required=True, help='Path to alignment matrix TSV')
    parser.add_argument('--edge-matrix', required=True, help='Path to save/load edge matrix CSV')
    parser.add_argument('--id-mapping', required=True, help='Path to ID mapping CSV')
    parser.add_argument('--output-dir', required=True, help='Directory for output files')
    parser.add_argument('--layout', default='kk3d', help='Layout algorithm (default: kk3d)')
    parser.add_argument('--min-contig-length', type=int, default=750, help='Minimum contig length')
    parser.add_argument('--min-weight', type=float, default=0.0, help='Minimum edge weight')
    
    args = parser.parse_args()
    
    # Create output directory if it doesn't exist
    output_dir = pathlib.Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

        # Process ID mapping
    idmap = pd.read_csv(args.id_mapping, sep=',', header=0)
    idmap['contig_id'] = idmap['contig_id'].apply(lambda x: x.split('[')[0])
    idmap['contig_id'] = idmap['contig_id'].apply(lambda x: x.split(' ')[0])

    # Process edge matrix
    edge_df = edges_from_matrix(args.matrix, args.edge_matrix, idmap, args.min_contig_length, args.min_weight)
    
    # Create graph
    G, vertex_order = create_igraph(edge_df)
    # Generate aligned labels and groups given our mapping.
    aligned_labels, aligned_groups = align_labels_and_groups(idmap, vertex_order)
    
    # Create output files
    base_name = f'igraph_asm_ava_homology_nucl_{args.layout}'
    output_html = output_dir / f'{base_name}.html'
    
    # Generate visualizations
    make_3d_plot(G, aligned_labels, aligned_groups, output_html, args.layout, include_edges=True)
    make_3d_plot(G, aligned_labels, aligned_groups, 
                 output_dir / f'{base_name}_no_edges.html', args.layout, include_edges=False)
    
    # Save network data
    graph_data = igraph_to_coords_file(G, aligned_labels, aligned_groups)
    with open(output_dir / 'network.json', 'w') as f:
        json.dump(graph_data, f, indent=4)
    print(f'Network data saved to {output_dir / "network.json"}')

if __name__ == '__main__':
    main()
