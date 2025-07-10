        
    
import networkx as nx
import json

from mapping import  get_a_matrix_threshold, get_master_regulator_list, get_TF_lists, generate_gene_id_name_map, get_sorted_gene_order
from genexpression_utils import fetch_fibroblast_data




class GeneUtils():
    def __init__(self, thresh=300):
        self.net_effect_thresh = -5
        self.a_matrix_adata, self.b_matrix = get_a_matrix_threshold(thresh)

        self.master_regulator_list = get_master_regulator_list()
        self.repressorlist, self.activatorlist, self.conflictedlist, tf_list = get_TF_lists()

        self.gene_id_name_map, self.gene_name_id_map = generate_gene_id_name_map()
        
        self.genes = get_sorted_gene_order()

    def print_network_metrics(self, G):
        print(f"Number of nodes: {G.number_of_nodes()}")
        print(f"Number of edges: {G.number_of_edges()}")
        print(f"Is directed: {G.is_directed()}")
        print(f"Is weighted: {'weight' in next(iter(G.edges(data=True)))[-1] if G.number_of_edges() > 0 else False}")
        print(f"Density: {nx.density(G):.4f}")
        print(f"Average degree: {sum(dict(G.degree()).values()) / G.number_of_nodes():.2f}")
        print(f"Average clustering coefficient: {nx.average_clustering(G):.4f}")
        print(f"Connected components: {nx.number_connected_components(G) if not G.is_directed() else 'N/A for directed graphs'}")
        print(f"Diameter: {nx.diameter(G) if nx.is_connected(G) else 'Graph not connected'}")
        print(f"Average shortest path length: {nx.average_shortest_path_length(G) if nx.is_connected(G) else 'Graph not connected'}")

    
    def get_all_targets(self, matrix, perturbation):
        gene_row_index = matrix[matrix.obs_names == perturbation , :]
        gene_row = gene_row_index.X

        if hasattr(gene_row, "toarray"):
            gene_row = gene_row.toarray().ravel()
        else:
            gene_row = gene_row.ravel()
            
        matching_mask = (gene_row == 1) 
        targets = matrix.var_names[matching_mask]

        return targets

    def get_gene_targets(self, matrix, perturbation, nodelist):

        # Genes infulenced by a transcription factor

        gene_row_index = matrix[matrix.obs_names == perturbation , :].copy()
        
        var_node_mask = matrix.var_names.isin(nodelist)
        
        gene_row = gene_row_index.X
        
        if hasattr(gene_row, "toarray"):
            gene_row = gene_row.toarray().ravel()
        else:
            gene_row = gene_row.ravel()
        
        matching_var_mask = var_node_mask & (gene_row == 1) 
        targets = matrix.var_names[matching_var_mask]
        print(targets)
        
        return targets

    def get_gene_targets_multi(self, perturblist, nodelist):

        # Genes infulenced by multiple factors
        var_node_mask = self.b_matrix.var['GeneStableID'].isin(nodelist)

        gene_row_indices = self.b_matrix[self.b_matrix.obs['TFStableID'].isin(perturblist)].copy()


        gene_rows = gene_row_indices.X
        source_order = gene_row_indices.obs['TFStableID']

        if hasattr(gene_rows, "toarray"):
            gene_rows = gene_rows.toarray()
            
        targetlist = []
        sourcelist = []
        for gene_row, source in zip(gene_rows, source_order):
            matching_var_mask = var_node_mask & (gene_row == 1) 
            targets = self.b_matrix.var[matching_var_mask]['GeneStableID']
            targetlist.append(targets)
            sourcelist.append(source)

        return targetlist, sourcelist

    def get_influencers(self, gene_target, subnet=[]):
        
        # Transcription factors influencing a Gene 
        gene_index = self.b_matrix.var_names == gene_target
        gene_cols = self.b_matrix.X[:, gene_index].copy()
        
        if hasattr(gene_cols, "toarray"):
            gene_cols = gene_cols.toarray().ravel()
        else:
            gene_cols = gene_cols.ravel()
        
        matching_obs_mask = gene_cols == 1
        
        gene_influencers = self.b_matrix.obs_names[matching_obs_mask]

        if len(subnet) > 0:
            subset_influencers = set(gene_influencers).intersection(set(subnet))
            gene_influencers = list(subset_influencers)

        
        net_effect = 0
        for gene_influencer in gene_influencers:
            if gene_influencer in self.activatorlist:
                net_effect += 1
            elif gene_influencer in self.repressorlist:
                net_effect -= 1
        
        # activations.X
        return gene_influencers, net_effect
    

    def calculate_net_effect(self, genelist):
        
        net_effect = 0
        for gene in genelist:
            if gene in self.activatorlist:
                net_effect += 1
            elif gene in self.repressorlist:
                net_effect -= 1
        return net_effect

    # get_all_targets(a_matrix_adata, 'ENSG00000129152')
    # targetset = get_gene_targets(a_matrix_adata, 'ENSG00000129152', fibroblast_expressed_genes)


    def construct_network_2d(self, subnet_x, subnet_y, HWG_subnetwork):
        '''
        a_matrix is an anndata object which we will filter the nodes and construct network based on
        '''
        a_matrix = self.a_matrix_adata
        if HWG_subnetwork:
            subnetwork_nodes = list(set(subnet_x).union(subnet_y).intersection(set(HWG_subnetwork)))
        else:
            subnetwork_nodes = subnet_x + subnet_y

            
        subnetwork_adata = a_matrix[a_matrix.obs_names.isin(subnet_x), a_matrix.var_names.isin(subnet_y)].copy()
        matrix = subnetwork_adata.X
        G = nx.from_numpy_array(matrix, create_using=nx.Graph())
        G.add_nodes_from(subnetwork_nodes)
        G.add_nodes_from(HWG_subnetwork)

        
        mapping = dict(zip(range(len(subnetwork_adata.obs_names)), subnetwork_adata.obs_names))
        G = nx.relabel_nodes(G, mapping)
        return G

    def construct_network_x_y(self, specific_subnet_x, specific_subnet_y, HWG_subnetwork):
        '''
        a_matrix is an anndata object which we will filter the nodes and construct network based on
        '''
        a_matrix = self.a_matrix_adata
            
        subnetwork_adata = a_matrix[a_matrix.obs_names.isin(specific_subnet_x), a_matrix.var_names.isin(specific_subnet_y)].copy()
        G = nx.Graph()
        matrix = subnetwork_adata.X.tocoo()
        for i, j in zip(matrix.row, matrix.col):
            source = subnetwork_adata.obs_names[i]
            target = subnetwork_adata.var_names[j]
            G.add_edge(source, target)

        
        # mapping = dict(zip(range(len(subnetwork_adata.obs_names)), subnetwork_adata.obs_names))
        # G = nx.relabel_nodes(G, mapping)
        return G

    def construct_network(self, specific_subnet, HWG_subnetwork):
        '''
        a_matrix is an anndata object which we will filter the nodes and construct network based on
        '''
        a_matrix = self.a_matrix_adata
        if len(HWG_subnetwork) > 0:
            subnetwork_nodes = list(set(specific_subnet).intersection(set(HWG_subnetwork)))
            print(len(subnetwork_nodes))
        else:
            subnetwork_nodes = specific_subnet

            
        subnetwork_adata = a_matrix[a_matrix.obs_names.isin(subnetwork_nodes), a_matrix.var_names.isin(subnetwork_nodes)].copy()
        matrix = subnetwork_adata.X
        G = nx.from_numpy_array(matrix, create_using=nx.Graph())
        G.add_nodes_from(subnetwork_nodes)
        G.add_nodes_from(HWG_subnetwork)

        
        mapping = dict(zip(range(len(subnetwork_adata.obs_names)), subnetwork_adata.obs_names))
        G = nx.relabel_nodes(G, mapping)
        return G

    def add_network_props(self, G, ntype):
        
        repressorlist = self.repressorlist
        activatorlist = self.activatorlist
        conflictedlist = self.conflictedlist
        master_regulator_list = self.master_regulator_list
        gene_id_name_map = self.gene_id_name_map
        

        highlight_color = "#a345f0"
        repressor_color = "#d77132"
        initiator_color = "#f83e56"

        

        # g1_color = '#0096c7'
        # g2_color = '#00b4d8'
        # g3_color = "#90e0ef"
        g1_color = "#6b7fa6"
        g2_color = "#86c9d7"
        g3_color = "#c6d6d8"

        

        graph_nodes = G.nodes
        graph_degree = G.degree
        eigenvector_centrality = nx.eigenvector_centrality(G)
        degree_centrality = nx.degree_centrality(G)
        for node in graph_nodes:
            # Basic Info
            if ntype == 'size':
                G.nodes[node]['size'] = 5 + 1 * G.degree[node] 
            elif ntype == 'base':
                G.nodes[node]['size'] = 1 + 0.5 * G.degree[node]
            G.nodes[node]['genename'] = gene_id_name_map.get(node, "NA")
            G.nodes[node]['degree'] = G.degree[node]
            G.nodes[node]['eigenvector'] = eigenvector_centrality.get(node, 'NA')
            G.nodes[node]['degree_centrality'] = degree_centrality.get(node, 'NA')

            
            # classify node
            if node in master_regulator_list:
                G.nodes[node]['type'] = 'g1'
                G.nodes[node]['color'] = g1_color
            elif node in repressorlist or node in activatorlist or node in conflictedlist:
                G.nodes[node]['type'] = 'g2'
                G.nodes[node]['color'] = g2_color
            else:
                G.nodes[node]['type'] = 'g3'
                G.nodes[node]['color'] = g3_color

            if node in repressorlist:
                G.nodes[node]['tftype'] = 'Repressor'
            elif node in activatorlist:
                G.nodes[node]['tftype'] = 'Activator'
            elif node in conflictedlist:
                G.nodes[node]['tftype'] = 'Conflicted'
            else:
                G.nodes[node]['tftype'] = 'Normal'
            
            # activity of node
            if graph_degree[node] > 0:
                G.nodes[node]['active'] = True
            else:
                G.nodes[node]['active'] = False


            # Activated, Repressed, Neutral
            if not(G.nodes[node].get('state')):
                G.nodes[node]['state'] = 'Neutral'
        return G

    def insert_perturbation(self, G, activation_nodes, ref_matrix):
        # initiator_color = "#f83e56"
        initiator_color = '#ff1760'
        # first_order_color = "#fc7ebd"
        first_order_color = "#ff59c5"


        for node in activation_nodes:
            G.add_node(node)

            
            
            # G.nodes[node]['size'] = 15 + 1 * G.degree[node] 
            G.nodes[node]['genename'] = self.gene_id_name_map.get(node, "NA")
            G.nodes[node]['color'] = initiator_color
            targets = get_all_targets(ref_matrix, node)
            print(f'Total number of targets{len(targets)}')
                    
            for target in targets:
                if target not in G.nodes:
                    continue

                eigenvector_centrality = nx.eigenvector_centrality(G)
                degree_centrality = nx.degree_centrality(G)
                    
                # G.nodes[target]['size'] = 40 + 10 * G.degree[target] 
                G.nodes[target]['size'] = 5 + 2 * G.degree[target] 
                G.nodes[target]['color'] = first_order_color
                G.nodes[target]['degree'] = G.degree[target]
                # G.nodes[target]['state'] = 'Activated'
                G.nodes[target]['eigenvector'] = eigenvector_centrality.get(target, 'NA')
                G.nodes[target]['degree_centrality'] = degree_centrality.get(target, 'NA')
                G.nodes[target]['active'] = True

                G.add_edge(node, target)

            eigenvector_centrality = nx.eigenvector_centrality(G)
            degree_centrality = nx.degree_centrality(G)
                
            G.nodes[node]['state'] = 'Activated'
            G.nodes[node]['size'] = 5 + 2 * G.degree[node] 
            G.nodes[node]['degree'] = G.degree[node]
            G.nodes[node]['eigenvector'] = eigenvector_centrality.get(node, 'NA')
            G.nodes[node]['degree_centrality'] = degree_centrality.get(node, 'NA')
            G.nodes[node]['active'] = True
        return G
        
    def update_perturbation_network(self, G, subnet):
        '''
        This function updates the network based on the status of each node.
        Activated nodes propagate the changes further
        Repressed nodes stop the flow
        inactive nodes that have net positive get activated.
        '''

        nodelist = list(G.nodes)
        
        for node in nodelist:
            state = G.nodes[node]['state']
            active = G.nodes[node]['active']
            # pure = G.nodes[node]['pure']
            tftype = G.nodes[node].get('tftype', 'Normal')

            connected_edges = list(nx.node_connected_component(G, node))
            
            influencers, net_effect = self.get_influencers(node, subnet=connected_edges)

            
            if net_effect > self.net_effect_thresh:
                if active:
                    continue
                    
                if state == 'Activated':
                    G.nodes[node]['active'] = True
                elif (state == 'Neutral' or state == 'Repressed') and not active:
                    G.nodes[node]['active'] = True
                    G.nodes[node]['state'] = 'Activated'

                gene_targets = self.get_all_targets(self.a_matrix_adata, node)
                for target in gene_targets:
                    if target not in G.nodes:
                        G.add_node(target)
                        G.nodes[target]['state'] = 'Activated'
                        G.nodes[target]['active'] = True
                    G.add_edge(node, target)
                # Here we assume that activate edges for anything that is not already active indicating newly activated nodes.
                # We don't increase the transcription of existing activated nodes.
        
            elif net_effect <= self.net_effect_thresh:
                G.nodes[node]['state'] = 'Repressed'
                G.nodes[node]['active'] = False
                edges_to_remove = list(G.edges(node))
                
                G.remove_edges_from(edges_to_remove)
        return G
                    
                

            
            # if net_effect <= self.net_effect_thresh:
                
            #     G.nodes[node]['state'] = 'Repressed'
            #     G.nodes[node]['active'] = False
            #     edges_to_remove = list(G.edges(node)) 
                
            #     G.remove_edges_from(edges_to_remove)
                
            # elif net_effect > self.net_effect_thresh:
            #     if not active and (tftype in ['Activator', 'Repressor']):
            #         G.nodes[node]['active'] = False
            #         G.nodes[node]['state'] = 'Activated'
                    
                    
            #     if tftype == 'Activator' and not pure and active:
            #         targetset = get_gene_targets(self.b_matrix, node, a_matrix_adata.obs_names)
            #         for target in targetset:
            #             G.add_edge(node, target)
            #             G.nodes[target]['state'] = 'Activated'
                
                    
            #     if pure and state == 'Activated':
                    
                    
            #     if state == 'Repressed' and active:
            #         edges_to_remove = list(G.edges(node)) 
            #         G.remove_edges_from(edges_to_remove)
                    
            #         G.nodes[node]['active'] = False
        return G        

        
    def activate_perturbation(self, G, activation_nodes, ref_matrix, initial=True):
        initiator_color = "#f83e56"

        for node in activation_nodes:
            G.add_node(node)
            
            targets = get_all_targets(ref_matrix, node)
            print(f'Total number of targets{len(targets)}')
                    
            for target in targets:
                if target not in G.nodes:
                    continue
                # TO-DO: FIX THIS FOR CELL TYPE SPECIFIC B MATRIX
                # influencers, net_effect = get_influencers(b_matrix, target)

                
                current_influencers = G.edges(target)
                
                net_effect = self.calculate_net_effect(current_influencers)
                print(net_effect)
                
                if net_effect <= 0:
                    G.nodes[target]['state'] = 'Repressed'
                else:
                    G.nodes[target]['state'] = 'Activated'
            
            edges_to_add = [(node, target) for target in targets]
            G.nodes[node]['state'] = 'Activated'
            G.nodes[node]['size'] = 15 + 1 * G.degree[node] 
            G.nodes[node]['genename'] = self.gene_id_name_map.get(node, "NA")
            G.nodes[node]['color'] = initiator_color


            if initial:
                G.nodes[node]['tftype'] = 'Perturbation'
            
            G.add_edges_from(edges_to_add)
        return G

    def update_network(self, G):
        new_activations = []
        for node in G.nodes:
            if G.nodes[node]['tftype'] == 'Perturbation':
                continue

            if G.nodes[node]['state'] == 'Repressed':
                G.remove_edges_from(list(G.edges(node)))
            elif G.nodes[node]['state'] == 'Activated':
                new_activations.append(node)

        return G, new_activations
                
            
    # def perturb_network(G, perturbed_nodes, b_matrix):
    #     # set of new perturbed nodes. This will be used as the new list for perturbation
    #     for node in perturbed_nodes:
    #         if not(G.nodes[node]['active']):
    #             G.nodes[node]['active'] = True
                
    #             if G.nodes[node]['tftype'] == 'Repressor':
    #                 G.remove_edges_from(list(G.edges(node)))
                    
    #             elif G.nodes[node]['tftype'] == 'Activator':
    #                 targets = get_gene_targets_multi(b_matrix, node, list(G.nodes))
    #                 # activations = b_matrix[b_matrix.obs_names == node, b_matrix.var_names.isin(list(G.nodes))].copy()
                    
    #             else:
    #                 # both unknown and normal nodes
    #                 pass
                
    #         # else:
    #         #     if G.nodes[node]['type'] == 'g1' or G.nodes[node]['type'] == 'g2':
    #         #         pass

    def export_graph_as_json(self, G, file_path=None):
        ll = [ ]
        for n in G.nodes():
            nodeattrs = G.nodes[n]
            genename = str(nodeattrs.get('genename', "NA"))
            if 'genename' in nodeattrs:
                del nodeattrs['genename']
            if genename == "NA" or None:
                genename = self.gene_id_name_map.get(str(n), str(n))
            ll.append({"data": {"id": str(n), "genename": genename, **nodeattrs}})

        # Convert to Cytoscape.js format
        cy_data = {
            "nodes": ll,
            "edges": [
                {"data": {"id": f"{u}-{v}", "source": str(u), "target": str(v), **G.edges[u, v]}} for u, v in G.edges()
            ]
        }

        # Save to JSON
        if file_path:
            with open(file_path, "w") as f:
                json.dump(cy_data, f, indent=2)
        return cy_data


### BACKUP FUNCS FOR JUPYTER NOTEBOOK
def get_all_targets(matrix, perturbation):
    gene_row_index = matrix[matrix.obs_names == perturbation , :]
    gene_row = gene_row_index.X

    if hasattr(gene_row, "toarray"):
        gene_row = gene_row.toarray().ravel()
    else:
        gene_row = gene_row.ravel()
        
    matching_mask = (gene_row == 1) 
    targets = matrix.var_names[matching_mask]

    return targets

def get_gene_targets(matrix, perturbation, nodelist):

    # Genes infulenced by a transcription factor

    gene_row_index = matrix[matrix.obs_names == perturbation , :].copy()
    
    var_node_mask = matrix.var_names.isin(nodelist)
    
    gene_row = gene_row_index.X
    
    if hasattr(gene_row, "toarray"):
        gene_row = gene_row.toarray().ravel()
    else:
        gene_row = gene_row.ravel()
    
    matching_var_mask = var_node_mask & (gene_row == 1) 
    targets = matrix.var_names[matching_var_mask]
    
    return targets

def get_gene_targets_multi(b_matrix, perturblist, nodelist):

    # Genes infulenced by multiple factors
    var_node_mask = b_matrix.var['GeneStableID'].isin(nodelist)

    gene_row_indices = b_matrix[b_matrix.obs['TFStableID'].isin(perturblist)].copy()


    gene_rows = gene_row_indices.X
    source_order = gene_row_indices.obs['TFStableID']

    if hasattr(gene_rows, "toarray"):
        gene_rows = gene_rows.toarray()
        
    targetlist = []
    sourcelist = []
    for gene_row, source in zip(gene_rows, source_order):
        matching_var_mask = var_node_mask & (gene_row == 1) 
        targets = b_matrix.var[matching_var_mask]['GeneStableID']
        targetlist.append(targets)
        sourcelist.append(source)

    return targetlist, sourcelist

def get_influencers(b_matrix, gene_target):    
    repressorlist, activatorlist, conflictedlist, tf_list = get_TF_lists(log=False)
    
    # Transcription factors influencing a Gene 
    gene_index = b_matrix.var_names == gene_target
    gene_cols = b_matrix.X[:, gene_index].copy()
    
    if hasattr(gene_cols, "toarray"):
        gene_cols = gene_cols.toarray().ravel()
    else:
        gene_cols = gene_cols.ravel()
    
    matching_obs_mask = gene_cols == 1
    
    gene_influencers = b_matrix.obs_names[matching_obs_mask]

    net_effect = 0
    for gene_influencer in gene_influencers:
        if gene_influencer in activatorlist:
            net_effect += 1
        elif gene_influencer in repressorlist:
            net_effect -= 1
    
    # activations.X
    return gene_influencers, net_effect

def calculate_net_effect(genelist):
    repressorlist, activatorlist, conflictedlist, tf_list = get_TF_lists(log=False)
    
    net_effect = 0
    for gene in genelist:
        if gene in activatorlist:
            net_effect += 1
        elif gene in repressorlist:
            net_effect -= 1
    return net_effect

# get_all_targets(a_matrix_adata, 'ENSG00000129152')
# targetset = get_gene_targets(a_matrix_adata, 'ENSG00000129152', fibroblast_expressed_genes)




def construct_network_2d(a_matrix, subnet_x, subnet_y, HWG_subnetwork):
    '''
    a_matrix is an anndata object which we will filter the nodes and construct network based on
    '''
    if HWG_subnetwork:
        subnetwork_nodes = list(set(subnet_x).union(subnet_y).intersection(set(HWG_subnetwork)))
    else:
        subnetwork_nodes = subnet_x + subnet_y

        
    subnetwork_adata = a_matrix[a_matrix.obs_names.isin(subnet_x), a_matrix.var_names.isin(subnet_y)].copy()
    matrix = subnetwork_adata.X
    G = nx.from_numpy_array(matrix, create_using=nx.Graph())
    G.add_nodes_from(subnetwork_nodes)
    G.add_nodes_from(HWG_subnetwork)

    
    mapping = dict(zip(range(len(subnetwork_adata.obs_names)), subnetwork_adata.obs_names))
    G = nx.relabel_nodes(G, mapping)
    return G

def construct_network(a_matrix, specific_subnet, HWG_subnetwork):
    '''
    a_matrix is an anndata object which we will filter the nodes and construct network based on
    '''
    if HWG_subnetwork:
        subnetwork_nodes = list(set(specific_subnet).intersection(set(HWG_subnetwork)))
    else:
        subnetwork_nodes = specific_subnet

        
    subnetwork_adata = a_matrix[a_matrix.obs_names.isin(subnetwork_nodes), a_matrix.var_names.isin(subnetwork_nodes)].copy()
    matrix = subnetwork_adata.X
    G = nx.from_numpy_array(matrix, create_using=nx.Graph())
    G.add_nodes_from(subnetwork_nodes)
    G.add_nodes_from(HWG_subnetwork)

    
    mapping = dict(zip(range(len(subnetwork_adata.obs_names)), subnetwork_adata.obs_names))
    G = nx.relabel_nodes(G, mapping)
    return G

def add_network_props(G, ntype):
    
    global repressorlist, activatorlist, conflictedlist, master_regulator_list, gene_id_name_map
    

    highlight_color = "#a345f0"
    repressor_color = "#d77132"
    initiator_color = "#f83e56"

    

    # g1_color = '#0096c7'
    # g2_color = '#00b4d8'
    # g3_color = "#90e0ef"
    g1_color = "#6b7fa6"
    g2_color = "#86c9d7"
    g3_color = "#c6d6d8"

    

    graph_nodes = G.nodes
    graph_degree = G.degree
    for node in graph_nodes:
        # Basic Info
        if ntype == 'size':
            G.nodes[node]['size'] = 50 + 7.5 * G.degree[node] 
        elif ntype == 'base':
            G.nodes[node]['size'] = 5 + 1 * G.degree[node]
        G.nodes[node]['genename'] = gene_id_name_map.get(node, "NA")

        
        # classify node
        if node in master_regulator_list:
            G.nodes[node]['type'] = 'g1'
            G.nodes[node]['color'] = g1_color
        elif node in repressorlist or node in activatorlist or node in conflictedlist:
            G.nodes[node]['type'] = 'g2'
            G.nodes[node]['color'] = g2_color
        else:
            G.nodes[node]['type'] = 'g3'
            G.nodes[node]['color'] = g3_color

        if node in repressorlist:
            G.nodes[node]['tftype'] = 'Repressor'
        elif node in activatorlist:
            G.nodes[node]['tftype'] = 'Activator'
        elif node in conflictedlist:
            G.nodes[node]['tftype'] = 'Conflicted'
        else:
            G.nodes[node]['tftype'] = 'Normal'
        
        # activity of node
        if graph_degree[node] > 0:
            G.nodes[node]['active'] = True
        else:
            G.nodes[node]['active'] = False


        # Activated, Repressed, Neutral
        if not(G.nodes[node].get('state')):
            G.nodes[node]['state'] = 'Neutral'
    return G

def insert_perturbation(G, activation_nodes, ref_matrix):
    # initiator_color = "#f83e56"
    initiator_color = '#ff1760'
    # first_order_color = "#fc7ebd"
    first_order_color = "#ff59c5"

    for node in activation_nodes:
        G.add_node(node)
        
        
        # G.nodes[node]['size'] = 15 + 1 * G.degree[node] 
        G.nodes[node]['genename'] = gene_id_name_map.get(node, "NA")
        G.nodes[node]['color'] = initiator_color
        targets = get_all_targets(ref_matrix, node)
        print(f'Total number of targets{len(targets)}')
                
        for target in targets:
            if target not in G.nodes:
                continue
                
            G.nodes[target]['size'] = 40 + 10 * G.degree[target] 
            G.nodes[target]['color'] = first_order_color
            # G.nodes[target]['state'] = 'Activated'

            G.add_edge(node, target)
            
        G.nodes[node]['state'] = 'Activated'
        G.nodes[node]['size'] = 900

    return G
        
    
def activate_perturbation(G, activation_nodes, ref_matrix, initial=True):
    initiator_color = "#f83e56"

    for node in activation_nodes:
        G.add_node(node)
        
        targets = get_all_targets(ref_matrix, node)
        print(f'Total number of targets{len(targets)}')
                
        for target in targets:
            if target not in G.nodes:
                continue
            # TO-DO: FIX THIS FOR CELL TYPE SPECIFIC B MATRIX
            # influencers, net_effect = get_influencers(b_matrix, target)

            
            current_influencers = G.edges(target)
            
            net_effect = calculate_net_effect(current_influencers)
            print(net_effect)
            
            if net_effect <= 0:
                G.nodes[target]['state'] = 'Repressed'
            else:
                G.nodes[target]['state'] = 'Activated'
        
        edges_to_add = [(node, target) for target in targets]
        G.nodes[node]['state'] = 'Activated'
        G.nodes[node]['size'] = 15 + 1 * G.degree[node] 
        G.nodes[node]['genename'] = gene_id_name_map.get(node, "NA")
        G.nodes[node]['color'] = initiator_color


        if initial:
            G.nodes[node]['tftype'] = 'Perturbation'
        
        G.add_edges_from(edges_to_add)
    return G

def update_network(G):
    new_activations = []
    for node in G.nodes:
        if G.nodes[node]['tftype'] == 'Perturbation':
            continue

        if G.nodes[node][sate] == 'Repressed':
            G.remove_edges_from(list(G.edges(node)))
        elif G.nodes[node][sate] == 'Activated':
            new_activations.append(node)

    return G, new_activations
            
        
# def perturb_network(G, perturbed_nodes, b_matrix):
#     # set of new perturbed nodes. This will be used as the new list for perturbation
#     for node in perturbed_nodes:
#         if not(G.nodes[node]['active']):
#             G.nodes[node]['active'] = True
            
#             if G.nodes[node]['tftype'] == 'Repressor':
#                 G.remove_edges_from(list(G.edges(node)))
                
#             elif G.nodes[node]['tftype'] == 'Activator':
#                 targets = get_gene_targets_multi(b_matrix, node, list(G.nodes))
#                 # activations = b_matrix[b_matrix.obs_names == node, b_matrix.var_names.isin(list(G.nodes))].copy()
                
#             else:
#                 # both unknown and normal nodes
#                 pass
            
#         # else:
#         #     if G.nodes[node]['type'] == 'g1' or G.nodes[node]['type'] == 'g2':
#         #         pass
            
def export_graph_as_json(G, file_path=None):
    ll = [ ]
    for n in G.nodes():
        nodeattrs = G.nodes[n]
        genename = str(nodeattrs.get('genename', "NA"))
        if 'genename' in nodeattrs:
            del nodeattrs['genename']
        if genename == "NA" or None:
            genename = gene_id_name_map.get(str(n), str(n))
        ll.append({"data": {"id": str(n), "genename": genename, **nodeattrs}})

    # Convert to Cytoscape.js format
    cy_data = {
        "nodes": ll,
        "edges": [
            {"data": {"id": f"{u}-{v}", "source": str(u), "target": str(v), **G.edges[u, v]}} for u, v in G.edges()
        ]
    }

    # Save to JSON
    if file_path:
        with open(file_path, "w") as f:
            json.dump(cy_data, f, indent=2)
    return cy_data
