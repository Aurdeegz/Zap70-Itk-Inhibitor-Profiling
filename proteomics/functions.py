## Importables
from platform import python_version
print(f"Python version {python_version()}")

import sys
import os

#sys.path.append("/windir/c/Users/redas/Desktop/jupyter_directory/helpers/src/helpers/")

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from sklearn.decomposition import PCA, NMF
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import StandardScaler, normalize
from math import log10, log2, ceil, floor, sqrt

# And grab the helpers
import sys
sys.path.append("/windir/c/Users/redas/Desktop/jupyter_directory/helpers/src/helpers/")
from helpers import general_helpers as gh
from helpers import stats_helpers as sh
from helpers import mpl_plotting_helpers as mph
from helpers.proteomics_helpers import Peptide
from helpers import argcheck_helpers as ah

# for the enrichment dotplots, will be moved eventually
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.patches import Circle

# For GCT file writing
## Need to keep only the input values and flanks

def ptmsea_transform(a_qvalue, a_foldchange, head = "No"):
    """
    This function takes a single q-value (or p-value) and the foldchange, and performs the following transformation:
    -log10(q)*(-1) if the foldchange is between 0 and 1
    -log10(q) if the foldchange is greater than 1
    Assumes the foldchanges are not log transformed
    """
    if type(a_qvalue) == str or type(a_foldchange) == str:
        return head
    if a_qvalue != a_qvalue or a_foldchange != a_foldchange:
        return float("nan")
    if 0 < a_foldchange < 1:
        return -10 * log10(a_qvalue) *-1
    else:
        return -10 * log10(a_qvalue) * 1

def duplicate_flanks(a_row,
                     flank_locs,   # list of indices for the places with flanking sequences
                     ):
    """
    This function takes the n-columns defining flanking sequences in a row and returns new rows where
    all the values are duplicated but contain a single flanking sequence. Builds new lists iteratively.

    Ex input: [Zap70, S491Y492Y493, KALGADDSYYTARSA , ALGADDSYYTARSAG, LGADDSYYTARSAGK] 
    Ex output: [[Zap70, S491Y492Y493, KALGADDSYYTARSA] ,
                [Zap70, S491Y492Y493, ALGADDSYYTARSAG] ,
                [Zap70, S491Y492Y493, LGADDSYYTARSAGK] ]
    """
    flank_num = len(flank_locs)
    new_rows = [[] for _ in range(flank_num)] # Make the new rows
    flanks = []
    f_ind = -1
    seen = False
    for i in range(len(a_row)):               # loop over the size of the row
        if i not in flank_locs:               # If this is not a flank
            for j in range(flank_num):
                new_rows[j].append(a_row[i])  # Add this info to each sublist
        elif i in flank_locs:                 # If you hit a flank index
            flanks.append(a_row[i])           # Add the flank to the flanks list
            if not seen:                      # if it's new
                for j in range(flank_num):    # Then add a placeholder 
                    new_rows[j].append("Flankity Flank Flank")
                seen = True                   # and mark as seen
            if f_ind == -1:                   # and update the flank index
                f_ind = i
    for i in range(flank_num):
        new_rows[i][f_ind] = flanks[i]        # Then add the flanks back in
    return new_rows

def dup_flanks_matrix(a_matrix, flank_locs):
    """
    Wraps duplicate_flanks for every row in a matrix, building a new matrix with the single flank columns
    """
    new_matr = []
    for row in a_matrix:
        exp = duplicate_flanks(row, flank_locs)
        new_matr += exp
    return new_matr

def write_gct_file(index,    # Flanking sequences with -p 
                    values,   # values per flanking sequence, len(values[i]) == len(headers), len(values) == len(index)
                    headers, # Headers for the columns
                    outfile = "bullshit.gct" # File to write, include path and .gct
                    ):
    """
    Writes a GCT file given the index of properly formatted flanking sequences (index), the values for each flanking sequence,
    and the headers. len(header) == len(values[i])

    returns None, writes a file to the given path instead. 
    """
    # GCT file has some info up top, we take it
    gct = [[fr"#1.3"], # Tells programs its a GCT 1.3 file
           [len(index), len(headers), 0, 0], # number of rows, number of cols, and metadata shit
           ["flanking_seq"] + headers] 
    for i in range(len(index)):
        gct.append([index[i]] + values[i])
    gct = [gh.list_to_str(row) for row in gct]
    gh.write_outfile(gct, outfile, writestyle = "w")
    return None







def filter_df_rows(a_df, column = "U_ID"):
    """
    Takes the rows of a PeptideDepot xls file (as a Pandas DataFrame) and a unique id, then filters the dataframe
    to include only one row with that unique id. The dataframe should be pre-sorted.
    """
    cols = list(a_df.columns)
    keycol = cols.index(column)
    a_list = list(a_df.to_numpy())
    first_occ = []
    first_key = []
    for row in a_list:
        if row[keycol] not in first_key:
            first_occ.append(row)
            first_key.append(row[keycol])
        else:
            pass
    return pd.DataFrame(first_occ, columns = cols)

def read_and_filter(a_file, columns):
    """
    Given an xls file (from peptide depot) and a list of columns for that file, filter the dataframe for 
    further processing. 
    """
    file = pd.read_excel(a_file)
    file = file[columns]
    file["U_ID"] = file[columns[0]] + file[columns[1]]
    fcs = list(file[columns[14:16]].astype(float).to_numpy())
    fcs = [[abs(value) for value in row] for row in fcs]
    qs = list(file[columns[16:18]].astype(float).to_numpy())
    # Originally replaced NANs with 1 in qs, but that's stupid I think.
    #qs = [replace_nans(row) for row in qs]
    #qs = [replcae_value(row, float('nan'), )]
    file["q_num"] = [sum([1 for _ in row if _ < 0.05]) for row in qs]
    file["min_q"] = [min(row) for row in qs]
    file["max_fc"] = [max(row) for row in fcs]
    file = file.sort_values(["U_ID", "q_num", "min_q", "max_fc"],
                            ascending = [True, False, True, False])
    file = filter_df_rows(file)
    return file 

def find_lims(all_qs,all_fcs):
    """
    Given a list of all q-values and all fold changes, return the max -log10(q) and
    foldchange (assumes foldchanges are log transformed)
    """
    qs = gh.unpack_list(all_qs)
    fc = [abs(item) for item in gh.unpack_list(all_fcs) if item == item]
    return round(-log10(min(qs)))+1, round(max(fc))+1






# PCA function to make my life less miserable
def square_bounds(mpl_axes, inplace = False):
    """
    This just finds the biggest difference and creates square axes, which is just
    nicer sometimes.
    """
    ticks = list(mpl_axes.get_xticks()) + list(mpl_axes.get_yticks())
    if inplace:
        mpl_axes.set_xlim(min(ticks), max(ticks))
        mpl_axes.set_ylim(min(ticks), max(ticks))
    else:
        return min(ticks), max(ticks)

def pca_analysis(a_df, pca_kwargs = dict(n_components = 2,
                                         whiten = False,
                                         svd_solver = "full",
                                         tol = 0)):
    """
    Given a dataframe and a dictionary with keyword arguments for scikit learns PCA,
    perform dimenstionality reduction and reutrn the components and the PCA object.

    Note that the components contain the coordinates for the dimensionally reduced
    points.
    """
    std_scalar = StandardScaler()
    scaled = std_scalar.fit_transform(a_df.transpose())
    pca_analysis = PCA(**pca_kwargs)
    components = pca_analysis.fit_transform(scaled)
    components = gh.transpose(*[list(point) for point in list(components)])
    return components, pca_analysis

def nmf_analysis(a_df, nmf_kwargs = dict(n_components = 2, 
                                     init = "nndsvd", # preserves sparseness
                                     solver = "mu", # multiplicative update
                                     beta_loss = "frobenius", # stable, but slow
                                     alpha_W = 0,  # default
                                     alpha_H = 0,  # default
                                     l1_ratio = 0  # default
                                    )):
    """
    Given a dataframe and a dictionary with keyword arguments for scikit learns NMF,
    perform dimenstionality reduction and reutrn the components and the NMF object.

    Note that the components contain the coordinates for the dimensionally reduced
    points.
    """
    #std_scalar = StandardScaler()
    #scaled = std_scalar.fit_transform(a_df.transpose())
    norms = normalize(a_df)
    nmf_analysis = NMF(**nmf_kwargs)
    W = nmf_analysis.fit_transform(norms)
    H = nmf_analysis.components_
    return H, W
    
def cluster_plotting(dataframes, # list with minimum 1 df
                 groups,     
                 expnames, # should match len(dataframes)
                 filenames,# should match len(dataframes)
                 group_slices, #assumes reps are clustered in list
                 labels,       # should correspons to len(group_slices)
                 colours,      # list of lists, each sublist should correspond to len(group_slices)
                 markers = ["o","^", "s"],
                 cluster = 'PCA', # other option is NNMF
                 markersize=100, 
                 square = True,
                     forced_axes = False,
                 textdict = dict(fontfamily = "sans-serif",
                 font = "Arial",
                 fontweight = "bold",
                 fontsize = 10),
                 pca_kwargs = dict(n_components = 2,
                                   whiten = False,
                                   svd_solver = "full",
                                   tol = 0),
                 nmf_kwargs = dict(n_components = 2, 
                                   init = "nndsvd", # preserves sparseness
                                   solver = "mu", # multiplicative update
                                   beta_loss = "frobenius", # stable, but slow
                                   alpha_W = 0,  # default
                                   alpha_H = 0,  # default
                                   l1_ratio = 0  # default
                                    )):
    """
    Function that wraps PCA or NMF analysis and makes a 2D plot with my specific graphing style
    """
    # Get the data columns from the dataframes, remove
    # missing values, run PCA analysis with sklearn,
    # scatter
    
    # Grab the columns corresponding to the groups of data. Assumes the 'groups' strings
    # are a substring of the column headers
    dfs = [df[[name for name in list(df.columns) if any(s in name for s in groups)]] for df in dataframes]
    # Remove any column with all missing values 
    dfs = [df.dropna(axis=1, how="all") for df in dfs]
    # then drop any row with missing values, as PCA doesn't tolerate MVs
    dfs = [df.dropna() for df in dfs]
    axes = []
    i = 0
    for df in dfs:
        if cluster == "PCA":
            components,pca = pca_analysis(df, pca_kwargs = pca_kwargs)
        else:
            # will add nmf soon
            components,nmf = nmf_analysis(df, nmf_kwargs = nmf_kwargs)
        fig, ax = plt.subplots(figsize = (6,6))
        # Next, loop over the slices and scatter
        j = 0
        for g in group_slices:
            if type(g) == list:
                newslice = slice(sorted(g)[0], sorted(g)[-1]+1)
                ax.scatter(components[0][newslice], components[1][newslice], 
                       color = colours[i][j],
                       marker = markers[j], 
                       s = markersize, 
                       alpha = 0.50,          # my preference
                       label = labels[j],
                       edgecolor = "black",   # my preference
                      )
            else:
                ax.scatter(components[0][g], components[1][g], 
                       color = colours[i][j],
                       marker = markers[j], 
                       s = markersize, 
                       alpha = 0.50,          # my preference
                       label = labels[j],
                       edgecolor = "black",   # my preference
                      )
            j+=1
        ax.set_title(expnames[i], **textdict)
        if cluster == "PCA":
            ax.set_xlabel(f"PC1 ({100*pca.explained_variance_ratio_[0]:.2f}%)",**textdict)
            ax.set_ylabel(f"PC2 ({100*pca.explained_variance_ratio_[1]:.2f}%)", **textdict)
            #square_bounds(ax, inplace = True)
        else:
            # will add nmf soon
            ax.set_xlabel(f"Component 1", **textdict)
            ax.set_ylabel(f"Component 2", **textdict)
        if square and not forced_axes:
            square_bounds(ax, inplace = True)
        elif not square and bool(forced_axes):
            print(forced_axes[i])
            ax.set_xlim(forced_axes[i][0][0], forced_axes[i][0][1])
            ax.set_ylim(forced_axes[i][1][0], forced_axes[i][1][1])
        mph.update_ticks(ax, which = "x")
        mph.update_ticks(ax, which ="y")
        ax.legend()
        plt.tight_layout()
        plt.savefig(filenames[i])
        axes.append(ax)
        plt.close()
        i+=1
    return axes





## Multiple Regression, because fuck doing the pairwise bullshit
reg_gs = ["thresholded timepoint1", "thresholded timepoint2", "thresholded timepoint3",
         "thresholded timepoint4", "thresholded timepoint5", "thresholded timepoint6",
         "thresholded timepoint7", "thresholded timepoint8", "thresholded timepoint9",]

def format_linreg_strs(coeffs, r2, intercept, label = "bd"):
    """
    Function that creates the equation of a line given n coefficients and an intercept.
    """
    outstr = f"{label}\n$y={intercept:.3f}"
    for i in range(len(coeffs)):
        outstr += fr"+{coeffs[i]:.3f}x_{{{i+1}}}"
    outstr += f"$\n$R={r2:.3f}$"
    return outstr

def plot_linreg_strs(strs, save = "test.pdf"):
    """
    Takes the output of format_linreg_strs and plots using matplotlib and LaTeX for formatting
    """
    fig, ax = plt.subplots()
    # Check how many strings there are, and adjust the axes accordingly
    num = len(strs)
    ax.set_xlim(-1,1)
    ax.set_ylim(-1,num)
    ax.set_yticks(list(range(num)))
    # turn off the bounding box
    ax.axis("off")
    # plot the strings
    for i in range(num):
        ax.text(0, i, strs[i], fontsize = 4, ha = "center", va = "center")
    plt.savefig(save)
    plt.close()

def _handle_missing_regression_reps(a_split_file):
    # First, detect a missing replicate
    new_splitf = []
    test = [[row for row in g if all([item == item for item in row])] for g in a_split_file]
    test = [i for i in range(len(test)) if len(test[i]) == 0]
    print(test)
    if test == []:
        return a_split_file
    else: # there's a missing replicate, so find it
        for i in range(len(a_split_file)):
            if i not in test:
                new_splitf.append(a_split_file[i])
            else:
                # Transpose and find the index of all missing values
                tmp = gh.transpose(*a_split_file[i])
                m_rep = [j for j in range(len(tmp)) if all([item != item for item in tmp[j]])]
                print(m_rep)
                tmp = [tmp[j] for j in range(len(tmp)) if j not in m_rep]
                new_splitf.append(gh.transpose(*tmp))
    return new_splitf
    
def multi_reg_lineplot(file, #dataframe
                       groups = ["t1", "t2", "t3"], #substring of header, group indicator
                       labels = ["0 min", "2 min", "5 min"], # Goes above the strings
                       log2_trans = True,
                       savefile = "test.pdf", # path/to/file.pdf
                       ):
    """
    Takes in a Pandas dataframe as well as some groups and labels, then plots the multiple regression
    output with LaTeX formatting. 
    """
    g_num = len(groups)
    split_f = [file[[c for c in list(file.columns) if groups[i] in c]] for i in range(g_num)]
    if log2_trans:
        split_f = [[list(row) for row in list(g.astype(float).transform(np.log2).to_numpy())] for g in split_f]
    else:
        split_f = [[list(row) for row in list(g.astype(float).to_numpy())] for g in split_f]
    # LinearRegression can't take missing values
    split_f = _handle_missing_regression_reps(split_f)
    split_f = [[row for row in g if all([item == item for item in row])] for g in split_f]
    print([len(row) for row in split_f])
    xs = [[row[:-1] for row in g] for g in split_f]
    ys = [gh.transpose(*[[row[-1]] for row in g])[0] for g in split_f]
    # Set up the model
    linmods = [LinearRegression() for _ in range(g_num)]
    # and fit it, always assume y is the last replicate in a group
    regs = [linmods[i].fit(xs[i], ys[i]) for i in range(g_num)]
    scores = [sqrt(regs[i].score(xs[i],ys[i])) for i in range(g_num)]
    # Now we make the strings
    strs = [format_linreg_strs(regs[i].coef_, scores[i], regs[i].intercept_, labels[i]) for i in range(g_num)]
    # And pass them to the plotter
    plot_linreg_strs(strs, save=savefile)
    return None




def pepdep_volcano_arr(file_dfs,
                       q_cols,
                       fc_cols,
                       wide = False,
                       **volcano_kwargs):
    # q-value column headers are the same from the standard peptide depot dump
    # so get a list of all the q-values
    all_qs = [[list(item) for item in list(df[q_cols].astype(float).to_numpy())] for df in file_dfs]
    all_fc = [[list(item) for item in list(df[fc_cols].astype(float).to_numpy())] for df in file_dfs]
    if wide:
        all_qs = gh.transpose(*[gh.transpose(*row) for row in all_qs])
        all_fc = gh.transpose(*[gh.transpose(*row) for row in all_fc])
    mph.volcano_array(all_qs, all_fc, **volcano_kwargs)
    return None

def make_all_pepplots(peptide_list, path = "outputs/graphics",
                      subset = ["all"], exclude = [], comparisons = [],
                      foldchange_group = None, heatname_tag = "",
                      global_max = 6,
                      heatmap_kwargs = {'aspect': 'equal', 'remove_spines': False, 'subplot_args': {'figsize': (14, 1)}, 
                                        'colorbar_args': {'orientation': 'vertical', 'location': 'right', 'shrink': 2}, 
                                        'textdict': {'fontfamily': 'sans-serif', 'font': 'Arial', 'fontweight': 'bold'}, 
                                        'img_name': 'figs/pep_plots/tbc1d5/heatmaps/tbc1d5_2_foldchange_all', 
                                        'heat_title': '', 
                                        'clb_label': 'log$_{2}$(FC)',
                                        'maxs': [-6, 6], 
                                        "sig_bounds" : [0.1,0.05,0.01],
                                        'cmap': mph.trans}):
    """
    """
    if global_max == None:
        max_val = math.ceil(max([max(p.vals) for p in peptide_list]))
        min_val = math.floor(min([min(p.vals) for p in peptide_list]))
        hm_max = math.ceil(max([max(p.means) for p in peptide_list]))
        hm_min = math.floor(min([min(p.means) for p in peptide_list]))
        hm = [-max([abs(hm_min), abs(hm_max)]), max([abs(hm_min), abs(hm_max)])]
        fc_max = math.ceil(max([max(p.fc) for p in peptide_list]))
        fc_min = math.floor(min([min(p.fc) for p in peptide_list]))
        fc = [-max([abs(fc_min), abs(fc_max)]), max([abs(fc_min), abs(fc_max)])]
    else:
        fc = [-global_max, global_max]
    for p in peptide_list:
        if not os.path.exists(f"{path}/{p.gene.lower()}/heatmaps"):
            os.makedirs(f"{path}/{p.gene.lower()}/heatmaps")
        p.heatmap(d_type = "foldchange",
                  path = f"{path}/{p.gene.lower()}/heatmaps/", 
                  maxs = fc,
                  heatname_tag = heatname_tag,
                  subset = subset, exclude = exclude,
                  heatmap_args = heatmap_kwargs)
        plt.close()











database_decoding = {
    "PERT-PSP" : "PhosphoSitePlus Perturbation\nSignatures",
    "DISEASE-PSP" : "PhosphoSitePlus Disease\nSignatures",
    "KINASE-PSP" : "PhosphoSitePlus Kinase\nSignatures",
    "KINASE-iKiP" : "In Vitro Kinase\nSubstrate Signatures",
    "PATH-BI" : "BI Pathway Signatures",
    "PATH-NP" : "NetPath Pathway\nSignatures",
    "PATH-WP" : "WikiPathway Signatures",
    "PERT-P100-DIA2" : "P100-DIA2\nPertubation Signatures",
}

def find_ids(str_list, delim = "_", position = 0, obj = []):
    new_ids = {}
    for s in str_list:
        newid = s.split(delim)[position]
        if newid not in list(new_ids.keys()):
            new_ids[newid] = obj
    return new_ids

def split_by_dbtype(sea_hm_file, delim = "_", id_col = 0, position = 0):
    # read the file, bin by id_col split on delim, return dict of sublists
    parsed = find_ids(gh.transpose(*sea_hm_file)[0][1:], delim = delim, position = position)
    parsed = {key : [sea_hm_file[0]] for key, value in parsed.items()}
    for row in sea_hm_file[1:]:
        newrow = [gh.list_to_str(row[0].split(delim)[1:], delimiter = " ", newline = False)] + [float(item) for item in gh.replace_value(row[1:],"NA", float("nan"))]
        parsed[row[id_col].split(delim)[position]].append(newrow)
    return parsed

def find_radius(sig, sig_mapper = {1 : 0.1,
                             0.05 : 0.2,
                             0.005 : 0.3,
                             0.0005 : 0.49}):
    if float(sig) != float(sig):
        return 0
    holder = 0.1
    for thresh in list(sig_mapper.keys()):
        if sig <= thresh:
            holder = sig_mapper[thresh]
    return holder

def score_to_colour(score, low = 0, high = 0.99999999, max_score = 1, 
                    cmap = cm.get_cmap("cool")):
    if float(score) != float(score):
        return "grey"
    mid = (low+high)/2
    moved_low = low - mid
    moved_high = high - mid
    scalar = moved_high/max_score
    color = (score*scalar) + mid
    return cmap(color)

def all_scores_to_colours(score_matrix, **stc_kwargs):
    # This should be a n (rows) x m (cols) matrix of scores (numbers)
    return [[score_to_colour(score, **stc_kwargs) for score in row] for row in score_matrix]

def find_all_radii(sig_matrix, sig_mapper = {1 : 0.1,
                             0.05 : 0.2,
                             0.005 : 0.3,
                             0.0005 : 0.49}):
    return [[find_radius(sig, sig_mapper = sig_mapper) for sig in row] for row in sig_matrix]

def sea_arr_poss(w,h, colours = None, radii = None, labels = None):
    if colours == None:
        colours = [["white" for j in range(w)] for i in range(h)]
    if radii == None:
        radii = [[0.1 for j in range(w)] for i in range(h)]
    arr = [[Circle((i,j), radii[i][j],ec="black", lw=0.5, color=colours[i][j]) 
            for j in range(w)] for i in range(h)]
    return arr

def legend_points(axis, sig_mapper, leg_scale, 
                  fontdict = dict(fontfamily="sans-serif",
                                      font = "Arial",
                                      fontweight= "bold",
                                      fontsize = 6)
                 ):
    # Make circles that are 0.5 apart with text next to them
    circs = []
    index = 0
    for key, value in sig_mapper.items():
        circs.append(Circle((1,index), radius = value, color = "white", ec = 'black' ) )
        index += 1
    keys = list(sig_mapper.keys())
    for i in range(len(circs)):
        axis.add_patch(circs[i])
        axis.text(1.5, i, f"$q < {keys[i]}$",**fontdict)
    return None

def add_legend(ax1,ax2, cmap, vmin = -10, vmax = 10, 
               leg_scale = 69,
               sig_mapper = {1 : 0.1,
                             0.05 : 0.2,
                             0.005 : 0.3,
                             0.0005 : 0.49},
               fontdict = dict(fontfamily="sans-serif",
                                      font = "Arial",
                                      fontweight= "bold",
                                      fontsize = 6),
               **stc_kwargs):
    # First, make fake circles that are white
    sig_matr = [[value for key, value in sig_mapper.items()]]
    rads= find_all_radii(sig_matr, sig_mapper)
    #points = sea_arr_poss(len(rads[0]), len(rads), labels = [[str(item) for item in row] for row in sig_matr])
    img = plt.imshow([[-10,10]], cmap = cmap, vmin = vmin, vmax = vmax)
    img.set_visible(False)
    cb = plt.colorbar(ax = ax1, fraction = 0.046, pad = 0.04)
    cb.set_ticks([item for item in cb.get_ticks() if item == int(item)])
    cb.set_ticklabels([int(item) for item in cb.get_ticks()], **fontdict)
    #cb.set_label(label, fontfamily = "sans-serif",
    #                  font = "Arial", fontweight= "bold", loc = "top")
    legend_points(ax2, sig_mapper, leg_scale, fontdict= fontdict)
    return ax1,ax2

def enrich_bubbleplot(enriched_dict,
                      savefile, # Should equal len(keys), be a directory to put files
                      filetype = "pdf",
                      group_heads = ["Large dong" for _ in range(20)],
                      bubblenum = 20, 
                      colourmap = mph.trans,
                      max_score = 10,
                      fontdict = dict(fontfamily="sans-serif",
                                      font = "Arial",
                                      fontweight= "bold",
                                      fontsize = 6),
                      vmin = -10, vmax = 10,
                      ):
    # Make an enrichment bubbleplot for every grouping, containing only
    # bubblenum groups.
    # Any filtering should be done preemptively.
    cmap = cm.get_cmap(colourmap)
    rgba = cmap(0.4999999995)
    # Loop over the keys and values of the dictionary
    saver = 0
    for key, value in enriched_dict.items():
        # The 0th row is headers, so ignore them.We will use group_heads for this
        # Grab the q-values and the enrichment scores. the number of
        # q-cols and nes-cols = len(group_heads)
        # Also, the 0th column is the labels
        qs = [row[1:len(group_heads)+1] for row in value[1:]]
        #qs = gh.transpose(*qs)
        qs = [qs[bubblenum*i:bubblenum*(i+1)] for i in range(len(qs)//bubblenum + 1)]
        qs = [q_list for q_list in qs if q_list != []]
        # For some reason, the above listcomp can produce empty lists at the end,
        # which breaks things downstream. So we'll just remove them
        #qs = [gh.transpose(*c) for c in qs]
        nes = [row[len(group_heads)+1:] for row in value[1:]]
        #nes = gh.transpose(*nes)
        nes = [nes[bubblenum*i:bubblenum*(i+1)] for i in range(len(nes)//bubblenum + 1)]
        # For some reason, the above listcomp can produce empty lists at the end,
        # which breaks things downstream. So we'll just remove them
        nes = [nes_list for nes_list in nes if nes_list != []]
        
        
        #nes = [gh.transpose(*c) for c in nes]
        labels = [row[0] for row in value[1:]]
        labels = [labels[bubblenum*i:bubblenum*(i+1)] for i in range(len(labels)//bubblenum + 1)]
        # The sea_arr_poss fails if there are no entries
        points = []
        #print(len(qs))
        #for i in range(len(qs)):
        #    print(qs[i])
        #    points.append(sea_arr_poss(len(qs[i]),    # Rows
        #                          len(qs[i][0]), # Cols
        #                          colours = all_scores_to_colours(gh.transpose(*nes[i]),
        #                                                          max_score=max_score,
        #                                                          cmap=cmap),
        #                          radii= find_all_radii(gh.transpose(*qs[i]))) )
        try:
            for i in range(len(qs)):
                points.append(sea_arr_poss(len(qs[i]),    # Rows
                                  len(qs[i][0]), # Cols
                                  colours = all_scores_to_colours(gh.transpose(*nes[i]),
                                                                  max_score=max_score,
                                                                  cmap=cmap),
                                  radii= find_all_radii(gh.transpose(*qs[i]))) )
        except:
            print(f"Skipping category {database_decoding[key]}: No groups were enriched.\n")
        if False:
            continue
        else:
            index = 0
            # Begin plotting each cluster
            for cluster in points:
                fig, ax = plt.subplots(1,2, figsize = (6,6))
                for row in cluster:
                    for p in row:
                        ax[0].add_artist(p)
                # Set some of the parameters
                ax[0].set_xticks(list(range(len(group_heads))))
                ax[0].set_xticklabels(group_heads, rotation = 90, **fontdict)
                ax[0].set_yticks(list(range(len(labels[index]))))
                ax[0].set_yticklabels(labels[index], **fontdict)
                ax[0].set_xlim(-1, bubblenum+1)
                ax[0].set_ylim(-1, bubblenum+1)
                ax[1].set_xlim(-1, bubblenum+1)
                ax[1].set_ylim(-1, bubblenum+1)
                ax[0].set_aspect("equal")
                ax[1].set_aspect("equal")
                ax[1].axis('off')
                ax[0].set_title(database_decoding[key],**fontdict)
                add_legend(ax[0], ax[1], cmap, leg_scale = bubblenum,
                          vmin = -max_score, vmax = max_score)
                ax[0].spines[:].set_visible(False)
                plt.tight_layout()
                plt.savefig(f"{savefile}/{key}_{index}.{filetype}")
                plt.close()
                index +=1
        saver += 1
    return None

def enrich_bubbleplot_list(ptmsea_outfiles, 
                           savefiles,
                           sig_exception = ["filenames"],
                           significance = 1, # 0 < significance < 1
                           filetype = "pdf",
                           group_heads = ["Large dong" for _ in range(20)],
                           bubblenum = 15, 
                           colourmap = mph.trans,
                           max_score = 10,
                           fontdict = dict(fontfamily="sans-serif",
                                      font = "Arial",
                                      fontweight= "bold",
                                      fontsize = 6),
                           vmin = -10, vmax = 10,):
    # Wraps bubbleplot to do a list of them and send them where they need to go
    for i in range(len(ptmsea_outfiles)):
        print(ptmsea_outfiles[i])
        file = gh.read_file(ptmsea_outfiles[i])
        file = split_by_dbtype(file)
        if ptmsea_outfiles[i] in sig_exception:
            print("exception")
            file = {key : [value[0]] + [row for row in value[1:] if any([item < 1 for item in row[1:len(group_heads)+1]])] 
               for key, value in file.items()}
        else:
            file = {key : [value[0]] + [row for row in value[1:] if any([item < significance for item in row[1:len(group_heads)+1]])]
               for key, value in file.items()}
            #print(file["PATH-NP"])
        # Remove any empties
        file = {key : value for key, value in file.items() if len(value) >1}
        enrich_bubbleplot(file, savefiles[i],
                         filetype = filetype,
                         group_heads = group_heads,
                         bubblenum = bubblenum,
                         colourmap = colourmap,
                          max_score = max_score,
                          fontdict = fontdict,
                          vmin = vmin, vmax = vmax)
        #break
    return None












