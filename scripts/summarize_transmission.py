#!/usr/bin/python3

import numpy as np
import pandas as pd
import networkx as nx
import argparse
import arviz as az

parser = argparse.ArgumentParser()
parser.add_argument('--inputFile', "-i",
                    help='input file is a TnT transmission analyser output log obtained on full trees. '
                         'Full trees are obtained either by assuming complete sampling or running stochastic mapping '
                         'on inferred trees.',
                    default="")
parser.add_argument('--truthFile', "-t", help='summary file of simulated truth if applicable.', default="")
parser.add_argument('--transmissionTimes', "-tt",
                    help='summary file of inferred transmission times as outputed by the mapper', default="")
parser.add_argument('--outputDir', "-o",
                    help='name for the output directory containing the summary of inferred transmissions.', default="")
parser.add_argument('--precision', "-p", help='precision in decimal points (default 2 results in format x.xx).',
                    type=float, default=2)
parser.add_argument('--probThreshold', "-pt", help='probability threshold to include edge in graph (default 0.1).',
                    type=float, default=0.1)
parser.add_argument('--hpdProb', "-hpd", help='probability for HPD calculation (default 0.95).',
                    type=float, default=0.95)
parser.add_argument('--noPlotDirect', '-nPD', help='do not plot the direct transmission graph', type=bool,
                    default=False)
parser.add_argument('--noPlotIndirect', '-nPI', help='do not plot the direct transmission graph', type=bool,
                    default=False)

args = parser.parse_args()
with_truth = True
with_times = True

if args.inputFile == "" or args.outputDir == "":
    print("Error, input and output files must be specified with -i and -o options.")
    exit()
if args.truthFile == "":
    with_truth = False
if args.transmissionTimes == "":
    with_times = False

tnt_analyser = pd.read_csv(args.inputFile, sep="\t")
tnt_analyser = tnt_analyser.loc[:, ~tnt_analyser.columns.str.match("Unnamed")]
hosts = np.unique([i.split('_', 1)[0] for i in tnt_analyser.columns])

if with_times:
    tnt_tr_times = pd.read_csv(args.transmissionTimes, sep="\t")
    tnt_tr_times = tnt_tr_times.loc[:, ~tnt_tr_times.columns.str.match("Unnamed")]

if with_truth:
    sim_truth_file = args.truthFile
    lookup_start = 'transmission times'
    lookup_end = 'SA count'
    start = False
    end = False
    from_ = []
    to_ = []
    time_ = []
    with open(sim_truth_file) as simFile:
        for line in simFile:
            if lookup_end in line:
                end = True
            if start and not end:
                pair = line.rstrip().split(':')[0]
                from_.append(pair.split('_')[0])
                to_.append(pair.split('_')[1])
                time_.append(line.rstrip().split(':')[1])
            if lookup_start in line:
                start = True

    df_sim = pd.DataFrame({'from': from_,
                           'to': to_,
                           'time': time_})

    # get direct transmission pairs that were sampled
    from_ = []
    to_ = []
    tr_times_in_sampled_hosts = []
    tr_times_in_sampled = []
    for index in df_sim.index:
        from_val = df_sim.loc[index, 'from']
        to_val = df_sim.loc[index, 'to']
        if to_val in hosts:
            tr_times_in_sampled_hosts.append(to_val)
            tr_times_in_sampled.append(df_sim.loc[index, 'time'])
            if from_val in hosts:
                from_.append(from_val)
                to_.append(to_val)

    df_direct_tr_in_sampled_truth = pd.DataFrame({'from': from_,
                                                  'to': to_})
    df_tr_times_in_sampled_truth = pd.DataFrame({'host': tr_times_in_sampled_hosts,
                                                 'time': tr_times_in_sampled})

    # get indirect transmission pairs that were sampled
    dict_to = dict(zip(df_sim['to'], df_sim['from']))
    from_in = []
    to_in = []
    n_intermediate_ = []
    intermediate_ = []
    source = ''
    # root = ''
    for h in hosts:
        n_int = 0
        if h == 'unsampled':
            continue
        sink = h
        tmp_intermediate = []

        if df_direct_tr_in_sampled_truth.empty or h not in df_direct_tr_in_sampled_truth['to'].values:
            source = dict_to.get(sink)
            while df_direct_tr_in_sampled_truth.empty or source not in df_direct_tr_in_sampled_truth['to'].values:
                if source not in dict_to.keys():  # or n_int >= hosts.size:
                    source = 'unsampled'
                    #               root=source
                    break
                sink = source
                source = dict_to.get(sink)
                n_int += 1
                tmp_intermediate.append(sink)

            if sink != '1':
                if n_int == 0:
                    from_.append(source)
                    to_.append(sink)
                else:
                    from_in.append(source)
                    to_in.append(h)
                    n_intermediate_.append(n_int)
                    intermediate_.append(tmp_intermediate)

    df_direct_tr_in_sampled_truth = pd.DataFrame({'from': from_,
                                                  'to': to_})

    intermediate_ = [' '.join(str(y) for y in x) for x in intermediate_]
    df_indirect_tr_in_sampled_truth = pd.DataFrame({'from': from_in,
                                                    'to': to_in,
                                                    'intermediate': intermediate_,
                                                    'n_intermediate': n_intermediate_})

    ##################################
    # save true transmission history #
    ##################################
    df_indirect_tr_in_sampled_truth.to_csv(args.outputDir + 'true_indirect_sampled_transmission.csv')
    df_direct_tr_in_sampled_truth.to_csv(args.outputDir + 'true_direct_sampled_transmission.csv')
    df_tr_times_in_sampled_truth.to_csv(args.outputDir + 'true_sampled_transmission_times.csv')

################################################################
# get inferred probabilities of different kinds of transmission #
################################################################


# prob of direct transmission:
dd = tnt_analyser.loc[:, :] == 1
prob_direct = dd.sum() / dd.shape[0]

# prob of indirect or direct transmission:
dd = tnt_analyser.loc[:, :] >= 1
prob_indirectAndDirect = dd.sum() / dd.shape[0]

# prob of indirect transmission:
dd = tnt_analyser.loc[:, :] > 1
prob_indirect = dd.sum() / dd.shape[0]

# prob of no transmission:
dd = tnt_analyser.loc[:, :] == 0
prob_noTr = dd.sum() / dd.shape[0]

# count intermediate unobserved transmissions
n_unobserved_transmissions = tnt_analyser[tnt_analyser.loc[:, :] >= 1] - 1

#######################
# indirect dataframes #
#######################
names_indirect = tnt_analyser.columns[prob_indirect > 0]
from_indirect = [i.split('_', 1)[0] for i in names_indirect]
to_indirect = [i.split('_', 1)[1] for i in names_indirect]

edges_indirect = pd.DataFrame({'from': from_indirect,
                               'to': to_indirect,
                               'probability': prob_indirect[prob_indirect > 0],
                               'median_unobserved': n_unobserved_transmissions.loc[:, names_indirect].median()})
nodes_indirect = pd.DataFrame({'id': np.unique([i.split('_', 1)[0] for i in from_indirect + to_indirect])})

# add edge labels, set precision
edges_indirect.loc[:, 'label'] = np.around(edges_indirect.loc[:, 'probability'],
                                           decimals=args.precision).astype(str)


###############################################
# direct and indirect transmission dataframes #
###############################################
names_indirectAndDirect = tnt_analyser.columns[prob_indirectAndDirect > 0]
from_indirectAndDirect = [i.split('_', 1)[0] for i in names_indirectAndDirect]
to_indirectAndDirect = [i.split('_', 1)[1] for i in names_indirectAndDirect]

inf_intermediates = []
hpd_lower_unobserved = []
hpd_upper_unobserved = []
for name in names_indirectAndDirect:
    inf_intermediates.append(n_unobserved_transmissions[n_unobserved_transmissions[name].notnull()][name].values)
    hpd = az.hdi(n_unobserved_transmissions[n_unobserved_transmissions[name].notnull()][name].values,
                                 hdi_prob=args.hpdProb)
    hpd_lower_unobserved.append(hpd[0])
    hpd_upper_unobserved.append(hpd[1])
# print(inf_intermediates)
# inf_intermediates = [' '.join(str(y) for y in x) for x in inf_intermediates]

edges_indirectAndDirect = pd.DataFrame({'from': from_indirectAndDirect,
                                        'to': to_indirectAndDirect,
                                        'probability': prob_indirectAndDirect[prob_indirectAndDirect > 0],
                                        # 'n_unobserved': inf_intermediates,
                                        'median_unobserved': np.nan_to_num(
                                            n_unobserved_transmissions.loc[:, names_indirectAndDirect].median()),
                                        'hpd_lower_unobserved': hpd_lower_unobserved,
                                        'hpd_upper_unobserved': hpd_upper_unobserved})
nodes_indirectAndDirect = pd.DataFrame(
    {'id': np.unique([i.split('_', 1)[0] for i in from_indirectAndDirect + to_indirectAndDirect])})

# add edge labels, set precision
edges_indirectAndDirect.loc[:, 'label'] = np.around(edges_indirectAndDirect.loc[:, 'probability'],
                                                    decimals=args.precision).astype(str)

# calculate root probability
root = []
tr_time = []
tr_time_hpd_lower = []
tr_time_hpd_upper = []
unsampled_exist = False
for h in nodes_indirectAndDirect['id']:
    if h != 'unsampled':
        root.append(1 - edges_indirectAndDirect.loc[edges_indirectAndDirect.loc[:, 'to'] == h, 'probability'].sum())
        if with_times:
            tr_time.append(tnt_tr_times[h].median())
            hpd = az.hdi(tnt_tr_times[h].values, hdi_prob=args.hpdProb, skipna=True)
            tr_time_hpd_lower.append(hpd[0])
            tr_time_hpd_upper.append(hpd[1])
    else:
        unsampled_exist = True
        root.append(0)
        if with_times:
            tr_time.append(np.NaN)
            tr_time_hpd_lower.append(np.NaN)
            tr_time_hpd_upper.append(np.NaN)
nodes_indirectAndDirect.loc[:, 'root_probability'] = np.around(root, decimals=3)
if with_times:
    nodes_indirectAndDirect.loc[:, 'median_infection_times'] = np.around(tr_time, decimals=3)
    nodes_indirectAndDirect.loc[:, 'hpd_lower_infection_time'] = tr_time_hpd_lower
    nodes_indirectAndDirect.loc[:, 'hpd_upper_infection_time'] = tr_time_hpd_upper
if unsampled_exist:
    nodes_indirectAndDirect.loc[nodes_indirectAndDirect.loc[:, 'id'] == 'unsampled', 'root_probability'] = \
        np.around(1 - sum(root), decimals=3)


#####################
# direct dataframes #
#####################
names_direct = tnt_analyser.columns[prob_direct > 0]
from_direct = [i.split('_', 1)[0] for i in names_direct]
to_direct = [i.split('_', 1)[1] for i in names_direct]
p = prob_direct[prob_direct>0]

# if it is from 'unsampled' the intermediates don't matter since they are also unsampled. So there is a direct
# transmission form an unsampled host.
for i in list(range(0,len(from_direct))):
    if from_direct[i] == 'unsampled':
        p[i] = \
            edges_indirectAndDirect.loc[(edges_indirectAndDirect['from'] == 'unsampled') & (edges_indirectAndDirect['to'] == to_direct[i]), 'probability']

edges_direct = pd.DataFrame({'from': from_direct,
                             'to': to_direct,
                             'probability': p})
nodes_direct = pd.DataFrame({'id': np.unique([i.split('_', 1)[0] for i in from_direct + to_direct])})

# add edge labels, set precision
edges_direct.loc[:, 'label'] = np.around(edges_direct.loc[:, 'probability'],
                                         decimals=args.precision).astype(str)

######################################
# save inferred transmission history #
######################################
edges_indirectAndDirect.to_csv(args.outputDir + 'inferred_transmission.csv')
edges_direct.to_csv(args.outputDir + 'inferred_direct_transmission.csv')
nodes_indirectAndDirect.to_csv(args.outputDir + 'inferred_host_summary.csv')

###############################
# plot the transmission graph #
###############################
# TODO : if no plotting of indirect - remove edges
# TODO: Add truth indicator to edges
if not args.noPlotDirect or not args.noPlotIndirect:
    import matplotlib as mpl

    # Apply probability threshold fo inclusion in the graph
    edges_indirectAndDirect_trh = edges_indirectAndDirect[edges_indirectAndDirect.probability >= args.probThreshold]

    # Build your graph. Note that we use the DiGraph function to create the graph!
    G = nx.from_pandas_edgelist(edges_indirectAndDirect_trh, 'from', 'to', ['probability', 'median_unobserved'],
                                create_using=nx.DiGraph())

    # Node colors by root probability
    cmap = mpl.cm.get_cmap("Greens")
    node_col = dict(zip(nodes_indirectAndDirect.id, nodes_indirectAndDirect.root_probability))
    for key, value in node_col.items():
        rgba = cmap(value)
        node_col[key] = mpl.colors.rgb2hex(rgba)

    nx.set_node_attributes(G,
                           node_col,
                           "fillcolor")

    # set edge width by transmission probability
    nx.set_edge_attributes(G,
                           nx.get_edge_attributes(G, 'probability'),
                           "penwidth")

    # merge transmission probability and mean intermediate transmissions for edge labels
    dct_rounded_probs = {k: round(v, args.precision) for k, v in nx.get_edge_attributes(G, 'probability').items()}
    dct_rounded_unobserved_count = {k: round(v, args.precision) for k, v in
                                    nx.get_edge_attributes(G, 'median_unobserved').items()}
    dct_label = {k: str(dct_rounded_probs[k]) + ', ' + str(dct_rounded_unobserved_count[k]) for k in
                 nx.get_edge_attributes(G, 'probability').keys()}

    nx.set_edge_attributes(G,
                           dct_label,
                           "label")

    A = nx.nx_agraph.to_agraph(G)
    A.layout(prog='dot')
    A.draw(args.outputDir + 'inferred_transmission.svg',
           args='-Grankdir=LR -Gsplines=true -Goverlap="false" -Nshape=hexagon -Nstyle=filled',
           prog='dot')

exit()
