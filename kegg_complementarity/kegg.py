import copy
import re
import pandas as pd
import requests
import itertools
import os.path
import numpy as np
import logging
from scipy.spatial import distance
import scipy.stats as sc
import plotly as py
import plotly.graph_objs as go

def parse_edges(filename):
    edges = pd.read_table(filename, header=0, index_col=None, sep='\t')
    tax_cols = []
    for el in list(edges.columns):
        try:
            if edges.ix[0, el][:3] == 'k__':
                tax_cols.append(el)
        except:
            pass
    for c in tax_cols:
        for i in edges.index:
            edges = edges.set_value(i, c, re.sub(r'_x(\w+)x$', r'_\1', edges.ix[i, c]))
            edges = edges.set_value(i, c, re.sub(r'_x(\w+)x(?=,)', r'_\1', edges[c][i]))
            edges = edges.set_value(i, c, re.sub(r'(?<=_)BRA(?=\w)', '', edges[c][i]))
            edges = edges.set_value(i, c, re.sub(r'KET$|KET(?=,)', '', edges[c][i]))
            edges = edges.set_value(i, c, re.sub(r'DASH', '-', edges[c][i]))
            edges = edges.set_value(i, c,
                                    re.sub(r'c__Erysipelotrichi$|c__Erysipelotrichi(?=,)', 'c__Erysipelotrichia',
                                           edges[c][i]))
            edges = edges.set_value(i, c, re.sub(r'f__Lachnospiraceae(?=,g__Ruminococcus)', 'f__Ruminococcaceae',
                                                 edges[c][i]))
    return (edges)
def pull_from_kegg(identifiers, ):
    print('Downloading data from KEGG')
    with open('kos.txt', 'a') as kos:
        with open('gene_links.txt', 'a') as gl:
            for id in identifiers:
                perc_now = round((float(identifiers.index(id)) / len(identifiers)) * 100, 2)
                perc_prev = round((float(identifiers.index(id) - 1) / len(identifiers)) * 100, 2)
                if (perc_now % 10 < 5) and (perc_prev % 10 > 5):
                    print(str(perc) + '%')
                req = 'http://rest.kegg.jp//link/ko/%s' % id
                temp = requests.get(req).text
                kos.write(temp + '\n')
                temp = temp.replace('\n', '').split('path:ko%s\tko:' % id)[1:]
                gene_links = [requests.get('http://rest.kegg.jp/link/genes/%s' % k).text for k in temp]
                # path_to_ko_counts[id] = {k: v.replace('\n', '').split('ko:%s\t' % k)[1:]
                #                          for k, v in zip(temp, gene_links)}
                gl.write('===' + id + '\n')
                for el in gene_links:
                    gl.write(el)
                gl.write('\n')
                # path_to_ko_counts[id] = {k:
                #                              [
                #                                  i.split(':')[0] for i in path_to_ko_counts[id][k]
                #                                  ]
                #                          for k in path_to_ko_counts[id].keys()}
    return ()
def parse_kegg_data():
    path_to_ko_counts = {}
    with open('kos.txt') as kos_txt:
        with open('gene_links.txt') as gl_txt:
            kos = {}
            new_path = True
            for i in kos_txt.read().splitlines():
                if new_path:
                    pathname = re.search(r'^path:(ko\d{5})\tko:(K\d{5})$', i)
                    kos[pathname.group(1)] = [pathname.group(2)]
                    new_path = False
                    continue
                if i == '':
                    new_path = True
                else:
                    kos[pathname.group(1)].append(re.search(r'^path:ko(\d{5})\tko:(K\d{5})$', i).group(2))
            path_to_ko_counts = {k[0]: {v: [] for v in k[1]} for k in list(kos.items())}
            for j in gl_txt.read().splitlines():
                if j == '':
                    continue
                if j[0] == '=':
                    path = j[3:]
                    continue
                new_entry = re.search(r'^ko:(K\d+)\t(\w+):.+$', j)
                path_to_ko_counts[path][new_entry.group(1)].append(new_entry.group(2))
    return (path_to_ko_counts)
	
############
# Turn KEGG-taxonomy to a acr->taxon dictionary
def tax_tree_to_acr_dict(tax_file):
    def get_full_taxonomy(index, taxons, number):
        perc = (float(number) / len(taxons)) * 100.
        perc_prev = (float(number - 1) / len(taxons)) * 100.
        if ((perc % 10) < 1) and ((perc_prev % 10) > 1):
            print(str(round(perc, 0)) + '%')
        acronym = taxons[number]
        dist_to_spec = min([number - sp for sp in index['s__'] if (number - sp) > 0])
        species = taxons[number - dist_to_spec]
        dist_to_genus = min([number - ge for ge in index['g__'] if (number - ge) > 0])
        genus = taxons[number - dist_to_genus]
        dist_to_family = min([number - fa for fa in index['f__'] if (number - fa) > 0])
        family = taxons[number - dist_to_family]
        dist_to_order = min([number - ord for ord in index['o__'] if (number - ord) > 0])
        order = taxons[number - dist_to_order]
        dist_to_cla = min([number - cl for cl in index['c__'] if (number - cl) > 0])
        cla = taxons[number - dist_to_cla]
        dist_to_phylum = min([number - ph for ph in index['p__'] if (number - ph) > 0])
        phylum = taxons[number - dist_to_phylum]
        dist_to_kingdom = min([number - ki for ki in index['k__'] if (number - ki) > 0])
        kingdom = taxons[number - dist_to_kingdom]
        taxonomy = '{},{},{},{},{},{},{}'.format(kingdom, phylum, cla, order, family, genus, species)
        out = {acronym[3:]: taxonomy}
        return (out)
    def preproc(lines):
        def test_if_ordered(spisok):
            levs = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J']
            ind_prev = 0
            ordered = True
            for i in range(len(spisok)):
                ind_now = levs.index(spisok[i][0])
                if (ind_prev - ind_now) < (-1):
                    ordered = False
                    return (ordered)
                ind_prev = ind_now
            return (ordered)
        temp = []
        levs = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J']
        # Make whole groups disappear based on words in their names
        omit_words = False
        for i in range(len(lines)):
            if omit_words:
                if any(x == lines[i][0] for x in remain_lev):
                    omit_words = False
                else:
                    continue
            if any(x in lines[i] for x in ['unclassified ', ' culture', ' symbiont', ' samples', 'Candidatus ']):
                remain_lev = levs[:levs.index(lines[i - 1][0]) + 1]
                omit_words = True
                continue
            temp.append(lines[i])
        # Skip non-taxon lines
        skip = []
        for i in range(len(temp)):
            fam_end = re.search(r'aceae$', temp[i])
            group_end = (re.search(r'eae$', temp[i]) and not (fam_end))
            if any(x in temp[i] for x in [' group', ' cluster', ' complex', 'subgen.']) or group_end:
                # prefix =
                # out.append(prefix + temp[i+1][1:])
                skip.append(i)
        # out = []
        for i in skip:
            prefix = temp[i][0]
            temp[i] = 'REMOVE THIS'
            j = i + 1
            while levs.index(temp[j][0]) > levs.index(prefix):
                temp[j] = levs[levs.index(temp[j][0]) - 1] + temp[j][1:]
                j += 1
        temp = [x for x in temp if x != 'REMOVE THIS']
        # Rename prefixes until they are in order
        ordered = False
        while not (ordered):
            move = False
            move_levs = []
            # temp = copy.copy(out)
            out = []
            for i in range(len(temp)):
                if temp[i + 1] == temp[-1]:
                    break
                if move:
                    if temp[i][0] in move_levs:
                        out.append(levs[levs.index(temp[i][0]) - dist + 1] + temp[i][1:])
                        continue
                    else:
                        move = False
                now_ind = levs.index(temp[i][0])
                next_ind = levs.index(temp[i + 1][0])
                dist = next_ind - now_ind
                out.append(temp[i])
                if not (dist <= 1):
                    move = dist
                    move_levs = levs[levs.index(temp[i][0]) + dist:]
                    continue
            ordered = test_if_ordered(out)
        out = ['_'.join(x.split()) for x in out]
        taxes = {'A': 'k__', 'B_': 'p__', 'C_': 'c__', 'D_': 'o__', 'E_': 'f__', 'F_': 'g__', 'G_': 's__', 'H_': 'a__',
                 'I_': 'a__', 'J_': 'a__'}
        final_out = []
        for el in out:
            temp = [re.sub(r'^%s' % k, v, el) for k, v in taxes.items() if v in re.sub(r'^%s' % k, v, el)][0]
            if re.search(r'^a__', temp):
                try:
                    temp = re.search(r'^a__\w+?(?=_)', temp).group(0)
                except:
                    print(temp)
                    print(el)
            if re.search(r'^s__', temp):
                temp = re.sub(r'(?<=__).+?_', '', temp)
                temp = re.sub(r'_\[TAX.+$', '', temp)
            final_out.append(temp)
        taxes = ['a__', 's__', 'g__', 'f__', 'o__', 'c__', 'p__', 'k__']
        for i in range(len(taxes) - 1):
            miss = taxes[i]
            pres = taxes[i + 1]
            to_erase = []
            for j in range(len(final_out) - 1):
                if final_out[j][:3] == pres and final_out[j + 1][:3] != miss:
                    to_erase.append(j)
            final_out = [final_out[k] for k in range(len(final_out)) if k not in to_erase]
        # Remove unmatched taxons from end of tree
        last_acr = 0
        for i in range(len(final_out)):
            if final_out[i][:3] == 'a__':
                last_acr = i
        if last_acr != (len(final_out) - 1):
            final_out = final_out[:(last_acr + 1)]
        return (final_out)
    print('Turning taxon tree to acronym dict')
    with open(tax_file) as f:
        lines = f.read().splitlines()
        euk = lines.index("AEukaryota")
        bact = lines.index("ABacteria")
        arch = lines.index("AArchaea")
        end = max(euk, arch)
        lines = lines[bact:end]
    index = {'k__': [], 'p__': [], 'c__': [], 'o__': [], 'f__': [], 'g__': [], 's__': [], 'a__': []}
    lines = preproc(lines)
    for i in range(len(lines)):
        prefix = lines[i][:3]
        index[prefix].append(i)
    kegg_tax = {}
    for el in index['a__']:
        kegg_tax.update(get_full_taxonomy(index, lines, el))
    return (kegg_tax)
	
# Turn model edges file into a list of all present taxons
def edges_to_tax_list(edges_file):
    edges = pd.read_table(edges_file, header=0, index_col=None, sep='\t')
    tax_cols = []
    for el in list(edges.columns):
        try:
            if edges.ix[0, el][:3] == 'k__':
                tax_cols.append(el)
        except:
            pass
    all_tax = [list(edges.ix[:, i]) for i in tax_cols]
    all_tax = [i for x in all_tax for i in x]
    all_tax = list(set(all_tax))
    all_tax = [re.sub(r'_x(\w+)x$', r'_\1', el) for el in all_tax]
    all_tax = [re.sub(r'_x(\w+)x(?=,)', r'_\1', el) for el in all_tax]
    all_tax = [re.sub(r'(?<=_)BRA(?=\w)', '', el) for el in all_tax]
    all_tax = [re.sub(r'KET$|KET(?=,)', '', el) for el in all_tax]
    all_tax = [re.sub(r'DASH', '-', el) for el in all_tax]
    all_tax = [re.sub(r'c__Erysipelotrichi$|c__Erysipelotrichi(?=,)', 'c__Erysipelotrichia', el) for el in all_tax]
    all_tax = [re.sub(r'f__Lachnospiraceae(?=,g__Ruminococcus)', 'f__Ruminococcaceae', el) for el in all_tax]
    return (all_tax)
	
# Create a dict of acr->taxon associations for a custom taxon list
def tax_to_acr_dict(tax_list, acr_dict):
    assoc = {}
    for tax in tax_list:
        for acr in acr_dict:
            if tax in acr_dict[acr]:
                try:
                    assoc[acr]
                except KeyError:
                    assoc[acr] = [tax]
                    continue
                assoc[acr].append(tax)
    return (assoc)
	
# Create a file with bulk characteristics of taxon metabolism
def pathway_taxon_chars(tax_tree, edges_files, identifiers):
    acr_dict = tax_tree_to_acr_dict(tax_tree)
    tax_list = []
    for el in edges_files:
        tax_list += edges_to_tax_list(el)
    tax_list = list(set(tax_list))
    assoc = tax_to_acr_dict(tax_list, acr_dict)
    if not (os.path.exists('kos.txt')):
        pull_from_kegg(identifiers)
    path_to_ko_counts = parse_kegg_data()
    for path in path_to_ko_counts:
        for ko in path_to_ko_counts[path]:
            path_to_ko_counts[path][ko] = [i for i in path_to_ko_counts[path][ko] if i in assoc]
    columns = ['Taxon', 'Pathway', 'KOs present (%)', 'Species within present KOs (% all species in this taxon)',
               'KOs total', 'Species total']
    counts_per_tax = pd.DataFrame(index=None, columns=columns)
    for tax in tax_list:
        acr_of_this_tax = [i for i in assoc if tax in assoc[i]]
        for path in identifiers:
            kos = path_to_ko_counts[path].keys()
            N_tot = len(kos)
            acr_by_ko = {k: list(set([v for v in path_to_ko_counts[path][k] if v in acr_of_this_tax])) for k in kos}
            N_non_emp = len([k for k in kos if len(acr_by_ko[k]) > 0])
            try:
                av_acr_non_emp = sum([len(acr_by_ko[k]) for k in kos if len(acr_by_ko[k]) > 0]) / float(N_non_emp)
                norm_av_acr = (float(av_acr_non_emp) / len(acr_of_this_tax)) * 100
            except:
                av_acr_non_emp = 0
                norm_av_acr = 0
            counts_per_tax.loc[len(counts_per_tax.index)] = [tax, path, round((N_non_emp / float(N_tot)) * 100, 1),
                                                             round(norm_av_acr, 2), N_tot, len(acr_of_this_tax)]
    counts_per_tax.to_csv('taxon_kegg.csv', index=None, sep='\t')
    return ()
	
# Take bulk metabolism info and add it to models
# as columns describing pathway fullness
def add_path_info_to_edges(edges_files, kegg_info_file):
    kegg = pd.read_table(kegg_info_file, header=0, index_col=None, sep='\t')
    for filename in edges_files:
        edges = parse_edges(filename)
        all_paths = set(list(kegg['Pathway']))
        if not(edges.empty):
            tax_cols = [x for x in edges.columns if ( (type(edges.ix[0, x]) is str) and ('k__' in edges.ix[0, x]) )]
            for el in all_paths:
                for i in range(len(tax_cols)):
                    edges[el + '_' + str(i)] = pd.Series(np.nan, index=edges.index)
            for i in edges.index:
                for path in all_paths:
                    for tax in tax_cols:
                        temp = kegg.ix[(kegg['Taxon'] == edges.ix[i, tax]) & (kegg['Pathway'] == path)].ix[:, 2:4]
                        temp = temp.product(axis=1)[temp.index[0]] / 100
                        edges.ix[i, path + '_' + str(tax_cols.index(tax))] = temp
        else:
            print('You got an empty edges list: ' + filename)
        new_filename = filename.replace('.txt', '_PATH_ADDED.txt')
        edges.to_csv(new_filename, sep='\t', index=None)
    return()
	
# Create tables with metabolic profiles of taxons
# for each pathway
def get_kegg_profiles(tax_tree, identifiers, edges_files =  None, tax_list = None):
    acr_dict = tax_tree_to_acr_dict(tax_tree)
    if not(edges_files is None):
        tax_list = []
        for el in edges_files:
            tax_list += edges_to_tax_list(el)
    tax_list = list(set(tax_list))
    assoc = tax_to_acr_dict(tax_list, acr_dict)
    if not (os.path.exists('kos.txt')):
        pull_from_kegg(identifiers)
    path_to_ko_counts = parse_kegg_data()
    pathway_profiles = {}
    for path in identifiers:
        columns = [x for x in path_to_ko_counts[path]] + ['Total']
        pathway_profiles[path] = pd.DataFrame(index=tax_list, columns=columns)
    for path in pathway_profiles:
        print('Starting calculations for: ' + path)
        path_table = pathway_profiles[path]
        for tax in path_table.index:
            acr_of_this_tax = [i for i in assoc if tax in assoc[i]]
            for KO in path_table.columns:
                if KO == 'Total':
                    L = len(acr_of_this_tax)
                    if L != 0:
                        # Get percentage
                        path_table.loc[tax, :] = (path_table.loc[tax, :]/float(L))
                    path_table.loc[tax, 'Total'] = L
                else:
                    temp = list(set([x for x in path_to_ko_counts[path][KO] if x in acr_of_this_tax]))
                    path_table.loc[tax,KO] = len(temp)
    return(pathway_profiles)
	
def compare_model_dist_to_global(edges_files_and_tax_code_and_used_taxons,pathways,all_dists_or_average,euc_or_jac):
    if all_dists_or_average == 'average':
        output = {k[0] : None for k in edges_files_and_tax_code_and_used_taxons}
    if all_dists_or_average == 'all':
        output = {k[0] : {} for k in edges_files_and_tax_code_and_used_taxons}
    for tup in edges_files_and_tax_code_and_used_taxons:
        print('=== Starting calculations for: ' + tup[0])
        edges = pd.read_table(tup[0],header=0,index_col=None, sep = '\t')
        coded_taxons = list(pd.read_table(tup[2],header=0,index_col=0, sep = '\t').index)
        code_to_tax = pd.read_table(tup[1], header=None, index_col=0, sep='\t')
        taxons = [code_to_tax.ix[id, 1] for id in coded_taxons]
        global_profiles = get_kegg_profiles('br08610.keg',pathways,tax_list=taxons)
        sample_pairs = [(edges.ix[x,'Predictor'],edges.ix[x,'Response']) for x in list(edges.index)]
        all_pairs = [x for x in itertools.combinations(taxons,2) if not( (x[0],x[1]) in sample_pairs or (x[1],x[0]) in sample_pairs)]
        av_dists_by_path = {k:{'global':None, 'edges':None, 'MW_pval': 1} for k in pathways}
        all_dists_by_path = {k:{'global':None, 'edges':None, 'MW_pval' : 1} for k in pathways}
        for path in pathways:
            gl_dists = []
            edg_dists = []
            for pairs in (all_pairs, sample_pairs):
                for el in pairs:
                    total1 = global_profiles[path].ix[el[0], -1]
                    total2 = global_profiles[path].ix[el[1], -1]
                    if not ((total1 == 0) or (total2 == 0)):
                        profile1 = global_profiles[path].ix[el[0], :-1]
                        profile2 = global_profiles[path].ix[el[1], :-1]
                        if euc_or_jac == 'euc':
                            temp = distance.euclidean(profile1,profile2)
                        if euc_or_jac == 'jac':
                            min_prof = [min(profile1[i], profile2[i]) for i in range(len(profile1))]
                            max_prof = [max(profile1[i], profile2[i]) for i in range(len(profile1))]
                            if (sum(max_prof) == 0):
                                temp = 0
                            else:
                                temp = 1 - (sum(min_prof) / sum(max_prof))
                        if pairs == all_pairs:
                            gl_dists.append(temp)
                        if pairs == sample_pairs:
                            edg_dists.append(temp)
            if euc_or_jac == ['euc']:
                L = len(list(global_profiles[path].ix[0, :-1]))
                max_dist = distance.euclidean([0]*L,[1]*L)
                gl_dists = [x / float(max_dist) for x in gl_dists]
            MW = sc.mannwhitneyu(gl_dists,edg_dists).pvalue
            if MW < 0.05:
                glob_av = sum(gl_dists) / len(gl_dists)
                edg_av = sum(edg_dists) / len(edg_dists)
                print("Grats! You've found something interesting!")
                print(path + ':  ' + str(MW))
                print("Non-model average distance: " + str(glob_av))
                print("Inter-model average distance: " + str(edg_av))
                if all_dists_or_average == 'average':
                    av_dists_by_path[path]['global'] = glob_av
                    av_dists_by_path[path]['edges'] = edg_av
                    av_dists_by_path[path]['MW_pval'] = MW
                if all_dists_or_average == 'all':
                    all_dists_by_path[path]['global'] = gl_dists
                    all_dists_by_path[path]['edges'] = edg_dists
                    all_dists_by_path[path]['MW_pval'] = MW
        if all_dists_or_average == 'all':
            output[tup[0]] = all_dists_by_path
        if all_dists_or_average == 'average':
            output[tup[0]] = av_dists_by_path
    return(output)
	
def generate_profiles_by_pathway(models_to_check,profiles_to_check):
    for model_file in models_to_check:
        model_to_check = pd.read_table(model_file, sep='\t', header=0, index_col=None)
        for c in ('Response', 'Predictor'):
            for i in model_to_check.index:
                model_to_check = model_to_check.set_value(i, c, re.sub(r'_x(\w+)x$', r'_\1', model_to_check.ix[i, c]))
                model_to_check = model_to_check.set_value(i, c, re.sub(r'_x(\w+)x(?=,)', r'_\1', model_to_check[c][i]))
                model_to_check = model_to_check.set_value(i, c, re.sub(r'(?<=_)BRA(?=\w)', '', model_to_check[c][i]))
                model_to_check = model_to_check.set_value(i, c, re.sub(r'KET$|KET(?=,)', '', model_to_check[c][i]))
                model_to_check = model_to_check.set_value(i, c, re.sub(r'DASH', '-', model_to_check[c][i]))
                model_to_check = model_to_check.set_value(i, c, re.sub(r'c__Erysipelotrichi', 'c__Erysipelotrichia', model_to_check[c][i]))
                model_to_check = model_to_check.set_value(i, c, re.sub(r'f__Lachnospiraceae,g__Ruminococcus', 'f__Ruminococcaceae,g__Ruminococcus', model_to_check[c][i]))
        for profile_file in profiles_to_check:
            profile_to_check = pd.read_table(profile_file, sep='\t', header=0, index_col=0)
            taxon_pairs = [(model_to_check.ix[i,'Predictor'],model_to_check.ix[i,'Response']) for i in model_to_check.index]
            interest_profiles = {k:(profile_to_check.ix[k[0],:],profile_to_check.ix[k[1],:]) for k in taxon_pairs if not ((profile_to_check.ix[k[1], -1] == 0) or (profile_to_check.ix[k[0], -1] == 0))}
            out = pd.DataFrame(data=None, columns = ['Taxon']+ list(interest_profiles[list(interest_profiles.keys())[0]][1].index))
            for pair in interest_profiles:
                out.loc[len(out)] = [pair[0]] + list(interest_profiles[pair][0])
                out.loc[len(out)] = [pair[1]] + list(interest_profiles[pair][1])
            filename = re.search(r'^[A-Z]+', model_file).group(0) + '_' + re.search(r'^ko\d+', profile_file).group(0) + '.txt'
            out.to_csv(filename,sep='\t',index=None)
    return()
	
	
vitamin_pathways = {
    'Thiamin': 'ko00730',
    'Riboflavin': 'ko00740',
    'Pyridoxine': 'ko00750',
    'Nicotinate': 'ko00760',
    'Pantothenate': 'ko00770',
    'Biotin': 'ko00780',
    'Lipoate': 'ko00785',
    'Folate': 'ko00790',
    'Retinol': 'ko00830',
    'Ubiquinone': 'ko00130'
}


identifiers = [y for x, y in vitamin_pathways.items()]
edges_files = ['HEALTHY_pair_models_0.001.txt', 'HEALTHY_triplet_models_0.001.txt', 'IBD_pair_models_0.001.txt', 'IBD_triplet_models_0.001.txt']

edges_files_and_tax_code_and_used_taxons = [
    ('IBD_0.3_0.001_3_0.3_0.001_pair_models.txt',
     'IBD_0.3_tax_code.txt',
     'IBD_0.3_train.txt'),
    ('HEAL_0.3_0.001_2_0.3_0.001_pair_models.txt',
     'HEAL_0.3_tax_code.txt',
     'HEAL_0.3_train.txt'),
]

euc = compare_model_dist_to_global(edges_files_and_tax_code_and_used_taxons,
                                 identifiers, 'all', 'euc')

temp = get_kegg_profiles('br08610.keg', ['ko00780'], edges_files =  ['IBD_0.3_0.001_3_0.3_0.001_pair_models.txt'])
temp['ko00780'].to_csv('x.txt', sep ='\t')
generate_profiles_by_pathway(['IBD_0.3_0.001_3_0.3_0.001_pair_models.txt'],['ko00780.txt'])