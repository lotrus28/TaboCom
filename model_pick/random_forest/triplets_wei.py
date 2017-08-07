import copy
import itertools
import re
import sys
from difflib import SequenceMatcher

import numpy as np
import pandas as pd
import scipy
import statsmodels.formula.api as sm

global tax_thr
tax_thr = 0

class Graph:
    def __init__(self):
        self.nodes = []
        self.len = 0
        self.wh_len = 0

    def add_node(self, n):
        present = self.find(n)
        if present is None:
            self.nodes.append(Node(n))
            self.nodes[-1].index = len(self.nodes) - 1
            self.len += 1
            self.wh_len += 1
            return (self.nodes[-1])
        else:
            return (present)

    def find(self, name):
        for n in self.nodes:
            if n.name == name:
                return (n)
        return (None)

    def add_edge(self, name1, name2):
        N = [name1, name2]
        present = [self.find(el) for el in N]
        for i in range(2):
            if present[i] is None:
                new_node = self.add_node(N[i])
                present[i] = new_node
        a = self.find(name1)
        b = self.find(name2)
        a.neigh.append(b)
        b.neigh.append(a)

    def add_weighted_edge(self, i, j, pairG):
        if (pairG.nodes[i] in pairG.nodes[j].neigh):
            self.nodes[i].neigh += [self.nodes[j]]
            self.nodes[j].neigh += [self.nodes[i]]

            self.nodes[j].neigh = list(set(self.nodes[j].neigh))
            self.nodes[i].neigh = list(set(self.nodes[i].neigh))

            self.nodes[j].weight[self.nodes[i]] += 1
            self.nodes[i].weight[self.nodes[j]] += 1

class Node:
    def __init__(self, name):
        self.name = name
        self.neigh = []
        self.color = 'white'
        self.index = None
        self.weight = {}

# No dots in dimnames!
# in: data with meta_data on pair graph (practically, an edge list), path to counts table
# out: edges list for truplets, table with triple models parameters
def reduce_to_triplets(fstats, counts):
    ##################################
    # Takes in a meta-df for double models!
    def table_to_graph(df):
        G = Graph()
        for i in range(1, len(df)):
            edge = df.iloc[i]
            G.add_edge(edge[0], edge[1])
        return (G)

    def tr_neighs(root):
        closest = [list(x) for x in itertools.combinations(root.neigh, 2)]
        far = {k: k.neigh for k in root.neigh if k.neigh != root}
        temp = []
        for el in far:
            for i in far[el]:
                temp.append([el] + [i])
        all = closest + temp
        all = [el for el in all if (root not in el)]
        all = [el + [root] for el in all]
        return (all)

    def remove_zeroes(df):
        # ((df.T == 2).all() == False).all()
        # is true, if no lines have zero values
        if not ((df.T == 0).all() == False).all():
            # Remember one zero line
            zero_line = df[(df.T == 0).all()].iloc[0]
            # Now remove all zeros
            df = df[(df.T != 0).any()]
            # And add our zero-line
            df.append(zero_line)
        return (df)

    def md_remove_outliers(df):
        inv_cov = df.cov().as_matrix()
        means = df.mean().as_matrix()
        md = df.apply((lambda x: scipy.spatial.distance.mahalanobis(x, means, inv_cov)), axis=1)
        # Q =scipy.stats.chi2.cdf(md, df.shape[1])
        Q = scipy.stats.chi2.ppf(0.975, df.shape[1])
        df = df[md < Q]
        return (df)

    def remove_outliers(df):
        q = df.quantile([0.025, 0.975])
        filt_df = df[(df.iloc[:, 0] > q.iloc[0, 0]) &
                     (df.iloc[:, 1] > q.iloc[0, 1]) &
                     (df.iloc[:, 2] > q.iloc[0, 2]) &
                     (df.iloc[:, 0] < q.iloc[1, 0]) &
                     (df.iloc[:, 1] < q.iloc[1, 1]) &
                     (df.iloc[:, 2] < q.iloc[1, 2])]
        return (filt_df)

    all_models = pd.DataFrame([np.nan] * 6,
                              index=['Response', 'Predictor1', 'Predictor2', 'Coef_p1', 'Coef_p2', 'Intercept']).T
    all_models.drop(0, 0)

    # Col: species; Row: samples
    do = table_to_graph(fstats)

    # Recursion depth issue
    # tr = copy.deepcopy(do)
    tr = table_to_graph(fstats)

    # Add zero weight to all edges
    # Later this weight will be showing the number of 3-models
    # with this particular pair of taxons

    # Then erase all edges, as we are making a new graph,
    # although containing all the vertices for double graph
    for n in tr.nodes:
        n.weight = {x: 0 for x in n.neigh}
        n.neigh = []
    print('Nodes in pair graph: ' + str(do.len))

    # For each node see all possible triplets that contain it
    # Then check if corresponding triplet linear models are better than pair-models
    # If they are --

    computed = []

    for i in range(do.len):
        sets = tr_neighs(do.nodes[i])
        # Get IDs of all vertices in all triplets to quickly look them up in our graph

        temp = copy.copy(sets)
        # Remove set if it has 1 black vertex
        # Black vertex means all triplets containing have been accounted for previously
        for el in temp:
            colors = [j.color for j in el]
            if 'black' in colors:
                sets.remove(el)

        # zip(set_indices, set_names)
        indices_and_names = [([el.index for el in j],[el.name for el in j]) for j in sets]

        # Now calculate models for the sets
        # for ind, na in zip(set_indices, set_names):
        for ind, na in indices_and_names:

            pairs = list(itertools.permutations([na[0],na[1],na[2]], 2))

            # Make sure there are no dublicates
            test = any([x[0] == x[1] for x in pairs])

            taxons = ['s__','g__', 'f__', 'o__', 'c__', 'p__']
            tax_name = taxons[tax_thr]

            if not(test):
                for pair in pairs:
                    temp = SequenceMatcher(None, pair[0],pair[1])
                    temp = temp.find_longest_match(0, len(pair[0]), 0, len(pair[1]))
                    lcs = pair[0][temp[0]:(temp[0] + temp[2])]
                    if (tax_name in lcs):
                        test = True
            if test:
                continue

            # How this line even works?
            # Is counts a df? Or what?

            temp = counts[[na[0], na[1], na[2]]]
            temp = remove_zeroes(temp)
            temp = remove_outliers(temp)

            orders = [[na[i]] + [p for p in na if p != na[i]] for i in range(len(na))]
            for order in orders:

                resp = order[0]
                pred1 = order[1]
                pred2 = order[2]

                text = resp + ' ~ ' + pred1 + ' + ' + pred2

                if ([resp,pred1,pred2] in computed) or ([resp,pred2,pred1] in computed):
                    continue

                computed.append(order)
                if not (temp.empty):
                    model = sm.ols(formula=text, data=temp).fit()

                    # Pick a threshold
                    # First get all F-stats for all 6 possible pair models within a triplet
                    sub_meta = fstats[
                        ((fstats['Response'] == na[0]) & (fstats['Predictor'] == na[1])) |
                        ((fstats['Response'] == na[0]) & (fstats['Predictor'] == na[2])) |
                        ((fstats['Response'] == na[1]) & (fstats['Predictor'] == na[2])) |
                        ((fstats['Response'] == na[1]) & (fstats['Predictor'] == na[0])) |
                        ((fstats['Response'] == na[2]) & (fstats['Predictor'] == na[0])) |
                        ((fstats['Response'] == na[2]) & (fstats['Predictor'] == na[1]))
                        ]

                    # Now pick the smallest one
                    # It is now your threshold for letting the triplet model in
                    cut = min(sub_meta.ix[:, 'F_Stat'])
                    if model.f_pvalue < cut:
                        new_row = pd.DataFrame(
                            [resp, pred1, pred2, model.f_pvalue, model.params[pred1], model.params[pred2], model.params['Intercept']],
                            index=['Response', 'Predictor1', 'Predictor2', 'F_Stat', 'Coef_p1', 'Coef_p2', 'Intercept']).T
                        all_models = all_models.append(new_row)
                        poss_edges = itertools.combinations(ind, 2)
                        # Add weights to our graph
                        for el in poss_edges:
                            tr.add_weighted_edge(el[0], el[1], do)
            # Finally, paint the seed vertex black: we are not coming back here
        do.nodes[i].color = 'black'

    # Remove zero-neighbor vertices
    tr.nodes = [el for el in tr.nodes if len(el.neigh) > 1]
    # Remove zero-weight edges from graph
    for el in tr.nodes:
        el.weight = {x: y for x, y in el.weight.items() if y > 0}
    return (tr, all_models)

# in: graph for triplets
# out: adjacency table w/out weights
def turn_3_graph_to_adj_table(graph):
    adj = []
    for n in graph.nodes:
        adj.append([n.name])
        for el in n.neigh:
            adj[-1] += [el.name]
    return (adj)


# in: graph for triplets
# out: list of edges with weights
def turn_3_graph_to_edge_list(graph):
    edges = []
    for el in graph.nodes:
        edges += [sorted([el.name, i.name]) + [el.weight[i]] for i in el.neigh]
    temp = [i[0] + '\t' + i[1] + '\t' + str(i[2]) for i in edges]
    temp = sorted(list(set(temp)))
    return (temp)


def return_tax_names(ID_to_taxon, profiles):
    tax_code = {}
    profiles = pd.read_table(profiles, header=0, index_col=0, sep='\t', engine = 'python')
    with open(ID_to_taxon) as f:
        for line in f.read().splitlines():
            temp = line.split()
            # tax ID is in the 1st column tax name -- 2nd
            tax_code[temp[0]] = temp[1]
    profiles.index = [tax_code[x] for x in profiles.index]
    return (profiles)


def code_tax_names(taxon_to_ID, profiles):
    tax_code = {}
    profiles = pd.read_table(profiles, header=0, index_col=0, sep='\t', engine = 'python')
    with open(taxon_to_ID) as f:
        for line in f.read().splitlines():
            temp = line.split()
            # tax ID is in the 1st column tax name -- 2nd
            tax_code[temp[1]] = temp[0]
    profiles.index = [tax_code[x] for x in profiles.index]
    return (profiles)


def fstat_to_triplet_edges(pair_mod, counts, path_out):
    # Make replacemens in data to avoid further confusion
    # E.g. KEGG considers Ruminococcus to be in Ruminococcacea
    # While Greengenes -- in Lachnospiraceae
    def prepair_counts_and_edges(counts, pair_mod):
        # counts.index = [re.sub(r'_x(\w+)x$', r'_\1', x) for x in counts.index]
        # counts.index = [re.sub(r'_x(\w+)x(?=,)', r'_\1', x) for x in counts.index]
        # counts.index = [re.sub(r'c__Erysipelotrichi$|c__Erysipelotrichi(?=,)', 'c__Erysipelotrichia', x) for x in
        #                 counts.index]
        # counts.index = [re.sub(r'f__Lachnospiraceae(?=,g__Ruminococcus)', 'f__Ruminococcaceae', x) for x in
        #                 counts.index]
        # pair_mod = pair_mod.replace({r'_x(\w+)x$': r'_\1'}, regex=True)
        # pair_mod = pair_mod.replace({r'_x(\w+)x(?=,)': r'_\1'}, regex=True)
        # pair_mod = pair_mod.replace({r'c__Erysipelotrichi$|c__Erysipelotrichi(?=,)': 'c__Erysipelotrichia'}, regex=True)
        # pair_mod = pair_mod.replace({r'f__Lachnospiraceae(?=,g__Ruminococcus)': 'f__Ruminococcaceae'}, regex=True)

        # Model calling can't work with commas or brackets in variable names
        counts.index = [x.replace(',', '00') for x in counts.index]
        counts = counts.T
        pair_mod['Response'] = pair_mod['Response'].str.replace(',', '00')
        pair_mod['Predictor'] = pair_mod['Predictor'].str.replace(',', '00')

        return (counts, pair_mod)

    with open(pair_mod) as f:
        if f.readline() == 'No significant models\n':
            fail = 'No pair models to base triple models on'
            outfile = path_out + 'all_triplet_models.txt'
            with open(outfile, 'w') as out:
                out.write(fail)
            return ()

    a = pd.read_table(pair_mod, header=0, index_col=None, engine = 'python')

    # NEW LINE!
    counts, a = prepair_counts_and_edges(counts, a)

    tr = reduce_to_triplets(a, counts)
    out = turn_3_graph_to_edge_list(tr[0])
    out = [i.replace('00', ',') for i in out]
    model_table = tr[1]

    model_table.index = range(0, model_table.shape[0])
    model_table.drop(0, 0, inplace=True)

    with open(path_out + 'triplet_edges.txt', 'a') as f:
        for line in out:
            f.write(line + '\n')

    out_model = path_out + 'all_triplet_models.txt'
    model_table = model_table[['Response', 'Predictor1', 'Predictor2', 'F_Stat', 'Coef_p1', 'Coef_p2', 'Intercept', ]]
    model_table = model_table.replace(re.compile(r'(?<=\w)00(?=\w)'), ',')
    model_table.to_csv(out_model, sep='\t', index=False)

    return(model_table)

p_pair_edges = sys.argv[1]
p_counts = sys.argv[2]
tax_code = sys.argv[3]
p_out = sys.argv[4]
if p_out[-1] != '/':
    p_out += '/'

# counts = return_taxonomies(p_counts, tax_code)
counts = return_tax_names(tax_code, p_counts)
x = fstat_to_triplet_edges(p_pair_edges, counts, p_out)
