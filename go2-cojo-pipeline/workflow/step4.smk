import os
import operator
import itertools
from pathlib import Path

import scipy # This isn't used, but we add it here because it seems to be an unlisted dependency of networkx
import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt



container: 'worker_3.2.3'


OUTDIR = config['outdir']




checkpoint create_network_graphs:
    input:
        comb=f'{{output}}/{{ethnicity}}/across-pheno/{OUTDIR}/comb.csv',
        to_run=f'{{output}}/{{ethnicity}}/across-pheno/{OUTDIR}/to_run_cond.csv',
        same_var=f'{{output}}/{{ethnicity}}/across-pheno/{OUTDIR}/same_var.csv',
        all_cma=f'{{output}}/{{ethnicity}}/across-pheno/{OUTDIR}/conditional/all.cma.cojo'
    output:
        edges=f'{{output}}/{{ethnicity}}/across-pheno/{OUTDIR}/edges.csv',
        graphs_json=f'{{output}}/{{ethnicity}}/across-pheno/{OUTDIR}/graphs.json',
        index_data=f'{{output}}/{{ethnicity}}/across-pheno/{OUTDIR}/index_data.csv',
        comb_index_data=f'{{output}}/{{ethnicity}}/across-pheno/{OUTDIR}/comb_index_data.csv'
    run:
        comb = pd.read_csv(input.comb)
        to_run = pd.read_csv(input.to_run)
        same_var = pd.read_csv(input.same_var)
        all_cma = pd.read_csv(input.all_cma, sep = '\t')

        # Find signals which are significant for multiple phenotypes
        # and change the variant ID
        for group in set(same_var['groups']):
            group_signals = same_var.loc[same_var['groups']==group]

            for var in set(group_signals['rightVar']):
                all_phenos = '_'.join(group_signals.loc[group_signals['rightVar']==var, 'rightPheno'].drop_duplicates().to_list())
                same_var.loc[(same_var['groups']==group) & (same_var['rightVar']==var), 'rightID'] = all_phenos + '_' + var

        comb2 = comb.copy()
        comb2.insert(1, 'display_ID', comb2['pheno_ID'])
        for _, (group, left, right) in same_var[['groups', 'leftID', 'rightID']].iterrows():
            comb2.loc[(comb2['groups']==group) & (comb2['pheno_ID']==left), 'display_ID'] = right

        nodes_all = comb2[['groups', 'pheno_ID', 'display_ID', 'p']].copy()

        # For variants with the same display_ID, keep the row with the lowest p-value
        nodes_all = nodes_all.sort_values(['groups', 'display_ID', 'p'])
        nodes = nodes_all[~nodes_all['display_ID'].duplicated(keep='first')].reset_index(drop=True)

        okay_all_cma = all_cma[all_cma['cond-flag'].isna()].drop(columns='cond-flag')
        m = to_run.merge(okay_all_cma, on = 'pairID')
        edges_all = m[['groups', 'pairID', 'leftID', 'rightID', 'to-condition', 'p', 'pC', 'cond-flag']].copy()

        # Add reason for dependency and flagging info to edges
        edges_all['reason'] = [list() for i in range(edges_all.shape[0])]
        for i in edges_all.loc[edges_all['pC'].isna(), 'reason']:
            i.append('collinear')

        # conditioned p-value above study-wide significance threshold
        for i in edges_all.loc[(1.3e-8<edges_all['pC']), 'reason']:
            i.append('insignif')

        edges_all['flag'] = [list() for i in range(edges_all.shape[0])]
        # For pairwise conditional result with no dependency,
        # flag pairs which have it's p-value attenuated by 4-orders of magnitude.
        for i in edges_all.loc[(4<((-np.log10(edges_all['p'])) - (-np.log10(edges_all['pC'])))) & (edges_all['reason'].str.len()==0), 'flag']:
            i.append('4-orders')

        edges_all['reason'] = edges_all['reason'].str.join(',')
        edges_all['flag'] = edges_all['flag'].str.join(',')

        edges_all.insert(edges_all.columns.to_list().index('leftID') + 1, 'left_display_ID', np.nan)
        edges_all.insert(edges_all.columns.to_list().index('rightID') + 1, 'right_display_ID', np.nan)

        for _id in set(pd.concat([edges_all['leftID'], edges_all['rightID']])):
            edges_all.loc[edges_all['leftID']==_id, 'left_display_ID'] = nodes_all.loc[nodes_all['pheno_ID']==_id, 'display_ID'].iloc[0]
            edges_all.loc[edges_all['rightID']==_id, 'right_display_ID'] = nodes_all.loc[nodes_all['pheno_ID']==_id, 'display_ID'].iloc[0]
            
        edges = edges_all[['groups', 'left_display_ID', 'right_display_ID', 'to-condition', 'pC', 'cond-flag', 'reason', 'flag']].copy()
        edges = edges[edges['reason']!=''].sort_values(['groups', 'left_display_ID', 'right_display_ID', 'pC'])
        edges = edges[~edges[['left_display_ID', 'right_display_ID']].duplicated(keep='first')].reset_index(drop=True)


        # Create the network graphs
        signals = {g: nx.DiGraph() for g in set(nodes['groups'])}
        for _, (g, pheno_id, display_id, pval) in nodes.iterrows():
            signals[g].add_node(display_id, p=pval, pheno_ID=pheno_id)

        for _, (group, left, right, to_cond, pval, cond_flag, reason, flag) in edges.iterrows():
            signals[group].add_edge(left, right, to_cond=to_cond, pC=pval, cond_flag=cond_flag, reason=reason, flag=flag)
            
            
        for graph in signals.values():
            subgraph = dict()
            traversed = set()
            num = 1

            def register_node(node, num):
                if node not in traversed:
                    if num not in subgraph:
                        subgraph[num] = {
                            'nodes': list(),
                            'index': None,
                            'p': None
                        }
                    subgraph[num]['nodes'].append(node)

                    traversed.add(node)

            for node in graph.nodes:
                register_node(node, num)
                connected_nodes = nx.descendants(graph.to_undirected(), node)
                # connected_nodes.update(nx.ancestors(graph.to_undirected(), node))
                for node2 in connected_nodes:
                    register_node(node2, num)
                num = max(subgraph.keys()) + 1

            for num, info in subgraph.items():
                node_pvals = sorted([(node, graph.nodes[node]['p']) for node in info['nodes']],
                                    key=operator.itemgetter(1))
                info['index'], info['p'] = node_pvals[0]

            graph.graph['subgraph'] = subgraph

        # Create the index table
        data = list()
        for group, graph in signals.items():
            subgraph = graph.graph['subgraph']
            for subgroup, v in subgraph.items():
                for node in v['nodes']:
                    index = node==v['index']
                    # p = graph.nodes[node]['p']
                    pheno_id = graph.nodes[node]['pheno_ID']
                    data.append([group, subgroup, pheno_id, node, index])

        # Use nx.node_link_graph to read the individual graphs        
        with open(output.graphs_json, 'w') as f:
            json.dump({group: nx.node_link_data(graph) for group, graph in signals.items()}, f)

        no_proxy = to_run.loc[to_run['to-condition'].isna(), ['groups', 'pairID', 'leftID', 'rightID', 'to-condition', 'cond-flag']].copy()
        no_proxy['p'] = np.nan

        cond_failed = all_cma.loc[all_cma['cond-flag'].notna(), ['pairID', 'p', 'cond-flag']].copy()
        cond_failed = to_run[['groups', 'pairID', 'leftID', 'rightID', 'to-condition']].merge(cond_failed, on = 'pairID')

        edges_all2 = pd.concat([edges_all, no_proxy, cond_failed])\
                            .sort_values(['groups', 'pairID', 'leftID', 'rightID', 'pC'])\
                            .reset_index(drop=True)

        edges_all2.loc[edges_all2['reason']=='', 'reason'] = np.nan
        edges_all2.loc[edges_all2['flag']=='', 'flag'] = np.nan

        assert to_run.shape[0] == edges_all2.shape[0], 'We have missing edge info!'
        edges_all2.to_csv(output.edges, index=False)


        index_data = pd.DataFrame(data, columns = ['groups', 'subgroup', 'pheno_ID', 'display_ID', 'is_index'])\
                                 .sort_values(['groups', 'subgroup', 'display_ID', 'pheno_ID'])
        index_data = index_data[['groups','subgroup','display_ID','pheno_ID','is_index']]
        index_data.to_csv(output.index_data, index=False)
        comb2_index = comb2.merge(index_data, how = 'left')
        comb2_index = comb2_index[['groups', 'subgroup', 'display_ID', 'pheno_ID', 'phenotype', 'ID', 'chrom', 'pos', 'b', 'se',
                                   'p', 'bC', 'bC_se', 'pC', 'type', 'is_index']]

        # Assign missing subgroup and is_index info to variants
        # with signif signal for multiple phenotypes but those that were not index variant because
        # it wasn't the lowest p-value variant.
        for display_id in set(same_var['rightID']):
            subgroup = index_data.loc[index_data['display_ID']==display_id, 'subgroup'].iloc[0]
            comb2_index.loc[comb2_index['display_ID']==display_id, 'subgroup'] = subgroup
        comb2_index.loc[comb2_index['is_index'].isna(), 'is_index'] = False
        comb2_index['subgroup'] = comb2_index['subgroup'].astype(int)
        comb2_index.sort_values(['groups', 'subgroup', 'display_ID', 'pheno_ID'], inplace=True)
        comb2_index.to_csv(output.comb_index_data, index=False)


rule create_graph_plot:
    input:
        rules.create_network_graphs.output.graphs_json
    output:
        f'{{output}}/{{ethnicity}}/across-pheno/{OUTDIR}/graph-plots/{{group}}.png'
    run:
        with open(input[0]) as f:
            signals = json.load(f)
        group = nx.node_link_graph(signals[wildcards.group])
        plt.figure(figsize=(10, 10), dpi=120)
        pos = nx.kamada_kawai_layout(group)
        nx.draw(group, pos, with_labels=True, arrowstyle="->", arrowsize=15, width = 1)
        plt.savefig(output[0])
        plt.close()


def input_create_all_graph_plots(w):
    graphs_json_file = checkpoints.create_network_graphs.get(output=w.output, ethnicity=w.ethnicity).output.graphs_json
    with open(graphs_json_file) as f:
        signals = json.load(f)
    return expand('{output}/{ethnicity}/across-pheno/{OUTDIR}/graph-plots/{group}.png',
                    output=w.output, OUTDIR=OUTDIR, ethnicity=w.ethnicity, group=signals.keys())

rule create_all_graph_plots:
    input:
        input_create_all_graph_plots
    output:
        touch(f'{{output}}/{{ethnicity}}/across-pheno/{OUTDIR}/graph-plots/.done')


# rule create_graph_plots:
#     input:
#         rules.create_network_graphs.output.graphs_json
#     output:
#         directory(f'output/{{ethnicity}}/across-pheno/{OUTDIR}/graph-plots')
#     run:
#         with open(input[0]) as f:
#             signals = json.load(f)
#         signals = {group: nx.node_link_graph(graph) for group, graph in signals.items()}
        
#         output_dir = Path(output[0])
#         output_dir.mkdir()
#         for group, graph in signals.items():
#             # https://stackoverflow.com/a/50459950
#             plt.figure(figsize=(10, 10), dpi=120)
#             pos = nx.kamada_kawai_layout(graph)
#             # pos = nx.shell_layout(graph)
#             # x_values, y_values = zip(*pos.values())
#             # x_max = max(x_values)
#             # x_min = min(x_values)
#             # x_margin = (x_max - x_min) * 0.25
#             # plt.xlim(x_min - x_margin, x_max + x_margin)
#             nx.draw(graph, pos, with_labels=True, arrowstyle="->", arrowsize=15, width = 1)
#             plt.savefig(output_dir.joinpath(f'{group}.png'))
#             plt.close()
