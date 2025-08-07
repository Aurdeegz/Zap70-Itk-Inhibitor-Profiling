###################################################################
#
# Making line plots and doing the statistics
#  #  #  # Add to mpl_plotting_helpers at some point

def _logical_ignore_comps(labelled_line_groups,
                          group_strs,
                          xgroup_strs):
    """
    Only want to compare along a line group (e.g. timecourse) or
    down an x-column (e.g. JE6 DMSO 0m vs JE6 U0126 0m), but not
    all the random other comparisons because statistically they're
    kind of useless
    
    So this function will find all of the pairs that are useless
    """
    groups_unpacked = []
    for group in labelled_line_groups:
        groups_unpacked += group
    # This will hold the ignored pairs
    ignore_me_senpai = []
    # First, get all pairs
    paired = gh.make_pairs(groups_unpacked,
                           dupes = False,
                           reverse = False)
    # Then iterate over and check the labels
    for p in paired:
        gs_check = 0
        xs_check = 0
        # Check all the group strings
        for gs in group_strs:
            if gs_check == 1:
                pass
            elif gs in p[0][0] and gs in p[1][0]:
                gs_check = 1
        # Check all the xgroup strings
        for xs in xgroup_strs:
            if xs_check == 1:
                pass
            elif xs in p[0][0] and xs in p[1][0]:
                xs_check = 1
        # If there isn't a match, in either, ignore
        if gs_check == 0 and xs_check == 0:
            ignore_me_senpai.append(p)
    # Return the ignored pairs at the end
    return ignore_me_senpai

def perform_line_statistics(labelled_line_groups,
                            ignore_comps,
                            comp_type,
                            statsfile):
    """
    labelled_line_groups -> data with labels
                            list of lists of [label, [d1,d2,...,dn]]
    ignore_comps -> list of pairs ("group 1", "group 2") to not be
                    compared
    comp_type -> statistics to use, currently only
                 ["HolmSidak", "TukeyHSD"] are supported
                 (both do an ANOVA first by default)
    statsfile -> a string to the output path and filename
                 for the statistics file output
    #####
    Returns None, just dumps the statsfile
    """
    assert comp_type in ["HolmSidak", "TukeyHSD"], f"Invalid comparison type: {comp_type}"
    groups_unpacked = []
    for group in labelled_line_groups:
        groups_unpacked += group
    if comp_type == "HolmSidak":
        comparison = sh.HolmSidak(*groups_unpacked,
                                  labels = True,
                                  override = True,
                                  alpha = 0.05,
                                  no_comp = ignore_comps)
    elif comp_type == "TukeyHSD":
        comparison = sh.TukeyHSD(*groups_unpacked,
                                  labels = True,
                                  override = True,
                                  alpha = 0.05,
                                  no_comp = ignore_comps)
    comparison.write_output(filename = statsfile,
                            file_type = "csv")
    return None

def find_centres(plotting_info):
    """
    plotting_info -> output from get_data_info, a list of
                     data info and the raw data
                     
    goal: grab the centres for xticks
    """
    centres = []
    for group in plotting_info:
        if len(centres) <= len(group[0]["centers"]):
            centres = group[0]["centers"]
    return centres

def line_plot(labelled_line_groups,
              show_points = False,
              show_legend = False,
              colours = ["grey" for _ in range(20)],
              group_labs = [f"Thing {i}" for i in range(20)],
              markers = ["s" for _ in range(20)],
              linestyles = ["solid" for _ in range(20)],
              xlabels = [f"Time {i}" for i in range(20)],
              ylabel = ["Fold change"],
              ylims = None,
              ignore_comps = [],
              statsfile = None,
              comp_type = "HolmSidak",
              figfile = None):
    """
    labelled_line_groups -> list of lists, where each sublist contains labelled groups
    """
    # First, get some basic plotting information
    plotting_info = [get_data_info(line) for line in labelled_line_groups]
    # Then manage the statistics
    if statsfile != None:
        perform_line_statistics(labelled_line_groups, 
                                ignore_comps, 
                                comp_type, 
                                statsfile)
    # Begin plotting c::
    if ylims == None:
        ylims = floor(min([item for item in gh.unpack_list(labelled_line_groups) if type(item) in [int, float]])), ceil(max([item for item in gh.unpack_list(labelled_line_groups) if type(item) in [int, float]]))
    # 
    fig, ax = plt.subplots(figsize = (6,6))
    # 
    for i in range(len(labelled_line_groups)):
        #
        ax.plot(plotting_info[i][0]["centers"],
                plotting_info[i][0]["means"],
                color = colours[i],
                label = group_labs[i],
                linestyle = linestyles[i])
        #
        for j in range(len(labelled_line_groups[i])):
            add_errorbar(ax, 
                         plotting_info[i][0]["centers"][j],
                         plotting_info[i][0]["means"][j],
                         plotting_info[i][0]["sems"][j],
                         color = colours[i])
            if show_points:
            #
                ax.scatter(plotting_info[i][0]["xs"][j],
                           plotting_info[i][1][j][1],
                           color = colours[i],
                           edgecolor = "black", alpha = 0.3,
                           marker = markers[i],
                           s = 10)
            else:
            #
                ax.scatter(plotting_info[i][0]["centers"],
                           plotting_info[i][0]["means"],
                           color = colours[i],
                           edgecolor = "black", alpha = 0.3,
                           marker = markers[i],
                           s = 30)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    xticks = find_centres(plotting_info)
    ax.set_xticks(xticks)
    ax.set_xticklabels(xlabels[:len(xticks)],
                       fontfamily = "sans-serif", 
                       font = "Arial", 
                       fontweight = "bold", 
                       fontsize = 12,
                       rotation = 45,
                       ha = "center")
    ax.set_ylim(*ylims)
    mph.update_ticks(ax, which = "y")
    ax.set_ylabel(ylabel, fontfamily = "sans-serif",
                  font = "Arial", fontweight = "bold",
                  fontsize = 14)
    if show_legend:
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5),
                  prop = mpl_fm.FontProperties(family = "sans-serif",
                                               weight = "bold"))
    if figfile == None:
        plt.show()
    else:
        plt.savefig(figfile)
    plt.close()
    return None

def replace_neg(a_list, value = float("nan")):
    """
    replace any value <0 with 0
    """
    newlist = []
    for item in a_list:
        try:
            truth = item < 0
        except:
            newlist.append(item)
        else:
            if truth:
                newlist.append(value)
            else:
                newlist.append(item)
    return newlist

def safe_log2(number):
    try:
        log2(number)
    except:
        return float("nan")
    else:
        return log2(number)
    
#
#
###################################################################