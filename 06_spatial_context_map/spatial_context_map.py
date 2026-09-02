"""
spatial_context_map.py

Spatial context map of ecotype adjacency (Fig. 1f), adapting the approach
of Hickey et al. for Visium ST data. For each spot, spatial neighborhoods
are defined by a 100-unit-radius Euclidean search. Within each neighborhood,
the relative frequency of each ecotype is computed and the dominant
combination contributing >= 90% of local composition is kept. Combinations
present in > 0.5% of neighborhoods are retained, and a directed graph is
built connecting combinations that differ by exactly one ecotype.

Input (place in ./data/):
    <section>_Coordinates.csv - per-spot x, y coordinates and ecotype ('cluster') label

NOTE: confirmed -- this script intentionally runs on a single
representative section (S20_3289_A3) for Fig. 1f, not pooled across
all 24 profiled sections.
"""

import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx
from collections import Counter
import itertools
import numpy as np

# Function to load real-world sample data with x, y coordinates and cluster labels
def load_sample_data(file_path):
    # Read the CSV file into a pandas DataFrame
    data = pd.read_csv(file_path, index_col=0)
    
    # Return the DataFrame
    return data


# Function to create combinations of cluster neighborhoods based on proximity
def create_combinations_from_coordinates(data, window_size=100, contribution_threshold=0.90):
    combinations = []
    
    # Iterate through all spots in the dataset
    for index, spot in data.iterrows():
        # Get the coordinates of the current spot
        x, y, cluster = spot['x'], spot['y'], spot['cluster']
        
        # Find neighbors within a 100-unit Euclidean radius (not k-nearest-neighbors)
        neighbors = data[((data['x'] - x)**2 + (data['y'] - y)**2) <= window_size**2]
        
        # Count occurrences of each cluster in the neighbors
        cluster_counts = neighbors['cluster'].value_counts()
        total_neighbors = cluster_counts.sum()
        print(cluster_counts)
        print(total_neighbors)

        # Keep adding clusters to the combination until the cumulative percentage exceeds the threshold
        cumulative_percentage = 0
        selected_clusters = []
        
        for cluster, count in cluster_counts.items():
            cumulative_percentage += count / total_neighbors
            selected_clusters.append(cluster)
            
            if cumulative_percentage >= contribution_threshold:
                break
        
        # Store this combination of clusters if it meets the 85% rule
        if selected_clusters:
            combinations.append(tuple(sorted(selected_clusters)))
            #print(combinations)
    
    return combinations

# Function to filter combinations based on frequency threshold
def filter_combinations(combinations, threshold=0.005):
    combination_counts = Counter(combinations)
    total_combinations = sum(combination_counts.values())
    
    # Filter combinations based on threshold
    filtered_combinations = {
        combo: count for combo, count in combination_counts.items()
        if count / total_combinations >= threshold
    }
    
    return filtered_combinations

# Function to build the network graph from combinations
def build_layered_network(combinations):
    G = nx.DiGraph()
    
    combination_counts = Counter(combinations)
    for combination, count in combination_counts.items():
        G.add_node(combination, size=count)
    
    for combo1 in combination_counts:
        for combo2 in combination_counts:
            if len(combo2) == len(combo1) + 1 and set(combo1).issubset(set(combo2)):
                G.add_edge(combo1, combo2)
    
    return G

# Function to plot the spatial context network and save as PDF
custom_cc_order = ['CC9', 'CC4', 'CC5', 'CC2', 'CC3', 'CC6', 'CC7', 'CC8']

def custom_sort_key(combination):
    # Prioritize combinations with CC9 by checking if 'CC9' is in the combination
    if 'CC9' in combination:
        return (0, [custom_cc_order.index(cc) if cc in custom_cc_order else len(custom_cc_order) for cc in combination])
    else:
        return (1, [custom_cc_order.index(cc) if cc in custom_cc_order else len(custom_cc_order) for cc in combination])

# Function to assign layers and control horizontal positioning based on custom CC order
def assign_layers_and_positions(G):
    pos = {}  # Position dictionary for the layout
    layer_positions = {}  # Track horizontal positioning for each layer

    # Group nodes by layer (based on the number of elements in the combination)
    layers = {}
    for node in G.nodes():
        layer = len(node)
        if layer not in layers:
            layers[layer] = []
        layers[layer].append(node)

    # Sort combinations within each layer based on the custom CC order
    for layer, nodes in layers.items():
        nodes.sort(key=custom_sort_key)  # Sort nodes in the layer using the custom sorting key

        # Use np.linspace to evenly space nodes in this layer along the x-axis
        num_nodes = len(nodes)
        x_positions = np.linspace(-num_nodes / 2, num_nodes / 2, num_nodes)
        
        # Assign each node in this layer to a horizontal position
        for i, node in enumerate(nodes):
            pos[node] = (x_positions[i], -layer)  # Negative layer to place nodes lower for higher layers

    return pos

# Function to plot the spatial context network with hierarchical layout and custom sorting
def plot_spatial_context_network_layers_save(G, community_colors, filename):
    # Assign layers and horizontal positions to each node
    plt.figure(figsize=(10, 7))
    pos = assign_layers_and_positions(G)

    # Draw edges (connections) between circles, no nodes
    nx.draw_networkx_edges(G, pos, edge_color='gray', width=1.0, alpha=0.6, arrows=False)

    # Draw each combination as colored triangles and black circle for frequency
    total_frequency = sum(nx.get_node_attributes(G, 'size').values())  # Calculate total frequency for normalization
    for node, (x, y) in pos.items():
        combination = node
        frequency = G.nodes[node]['size']
        percentage = frequency / total_frequency  # Normalize the frequency to get a percentage
        size = percentage * 500  # Scale the size of the circle based on the percentage

        # Stack triangles vertically above the circle
        for i, community in enumerate(combination):
            dy = i * 0.15  # Vertical offset for stacking triangles
            plt.scatter(x, y + dy, marker='^', s=150, color=community_colors[community], zorder=2)

        # Draw the black circle representing frequency BELOW the triangles
        plt.scatter(x, y - 0.2, s=size, color='black', zorder=1)

    # Legend for triangle colors
    legend_elements = [plt.Line2D([0], [0], marker='^', color='w', label=community, markerfacecolor=color, markersize=10)
                       for community, color in community_colors.items()]
    color_legend = plt.legend(handles=legend_elements, title='CCs', loc='center left', bbox_to_anchor=(1, 0.5))

    # Add the first legend (community colors)
    plt.gca().add_artist(color_legend)

    # Legend for black circle sizes (frequency representation)
    for size in [50, 100, 200]:  # Example sizes (adjust these based on your actual data scale)
        plt.scatter([], [], c='black', alpha=0.6, s=size, label=f' {round(size / 500 * 100, 1)}%')
    plt.legend(title='Frequency', loc='center left', bbox_to_anchor=(1, 0.1))

    plt.title('Spatial Context Network')

    # Save as PDF
    plt.savefig(filename, format='png', bbox_inches='tight', dpi=400)
    plt.show()



# Main function to load data, create combinations, and plot the network
def main(file_path, community_colors, output_path):
    # Step 1: Load the real-world data
    data = load_sample_data(file_path)
    
    # Step 2: Create combinations of clusters based on spatial proximity (100 nearest neighbors)
    combinations = create_combinations_from_coordinates(data, window_size=100)
    
    # Step 3: Filter combinations by frequency
    filtered_combinations = filter_combinations(combinations, threshold=0.005)
    
    # Step 4: Build the network from combinations
    G = build_layered_network(filtered_combinations)
    
    # Step 5: Plot the network and save the figure
    plot_spatial_context_network_layers_save(G, community_colors, output_path)



import os
os.makedirs("results", exist_ok=True)

if __name__ == "__main__":
    # Define the color map for clusters (adapt this to your data)
    community_colors = {

 'CC9': '#8F30A1',  
'CC4' : '#F0027F', 
'CC5' : '#BF5B17', 

'CC2' :'#7FC97F', 
'CC3' : '#FDC086', 
'CC1' : '#666666', 

'CC6' : "#046307",
'CC7' : '#17DEEE', 
'CC8' :  '#FFD700', 

'CC10' :'#6495ED'
        
    }

    # Call the main function with the file path to your real-world data
    main("data/S20_3289_A3_Coordinates.csv", community_colors, "results/S20_3289_A3_spatial_context_network.png")
