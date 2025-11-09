"""
=================================================================================================
Kenneth P. Callahan

9 July 2021
                  
=================================================================================================
Python >= 3.8.5

homebrew_stats.py

This module is meant to help with general statistical functions. Currently, there is only
a small number of statistics options supported, but I suspect this will grow in the future.

Currently supported:

    FDR Estimation (Storey)
    T-tests

=================================================================================================
Dependencies:

 PACKAGE     VERSION
  Pandas  ->  1.2.3
  Numpy   ->  1.20.1
  SciPy   ->  1.6.1

=================================================================================================
"""
print(f"Loading the module: helpers.homebrew_stats\n")
#####################################################################################################################
#
#   Importables

import fnmatch                                        # Unix-like string searching

import pandas as pd                                   # General use for data
import numpy as np                                    # General use for data

import scipy
import scipy.stats as ss                              # General use for Statisitcs
 
from scipy.interpolate import splrep, splev           # Used for Storey Q-value estimation, fitting cubic spline
from scipy.interpolate import UnivariateSpline

from statsmodels.stats.multicomp import pairwise_tukeyhsd # Used for Tukey-HSD test
from statsmodels.stats.multitest import multipletests     # Not currently in use, but I don't want to forget it's there

print(f"nummpy        {np.__version__}")
print(f"scipy         {scipy.__version__}")
print(f"pandas        {pd.__version__}\n")

#
#
#####################################################################################################################
#
#   Q-Value Estimation Algorithms

#### Storey

def storey_check_groups(groups):
    
    """
    =================================================================================================
    storey_check_groups(groups)
    
    This function is meant to check the groups argument input into the function storey()
    
    =================================================================================================
    Arguments:
    
    groups  ->  Either a list, tuple, pandas DataFrame, or a numpy array that describes the groups
                in Storey FDR estimation.
    
    =================================================================================================
    Returns: A list of lists that describe the groups used in Storey FDR estimation
    
    =================================================================================================
    """
    # If the input groups are a pandas DataFrame
    if type(groups) == type(pd.DataFrame()):
        # Then convert the groups into a transposed numpy array
        groups = groups.to_numpy().transpose()
        # and use list comprehension to reformat the groups into
        # a list of pairs.
        groups = [[groups[j][i] for j in range(len(groups))] 
                  for i in range(len(groups[0])) ]
    # If the input groups are a lsit
    elif type(groups) == list:
        # Then loop over the number of lists
        for i in range(len(groups)):
            # If the type of the input groups are not
            # a list, tuple or array
            if type(groups[i]) not in [list, tuple, type(np.array([]))]:
                # Then list the element
                groups[i] = list(groups[i])
            # Otherwise
            else:
                # Just keep the list the same
                groups[i] = groups[i]
    # If the groups were given as a tuple
    elif type(groups) == tuple:
        # Then turn the groups into a lsit
        groups = list(groups)
        # and loop over the number of items in the groups
        for i in range(len(groups)):
            # and if the element is not a list, tuple, array,
            if type(groups[i]) not in [list, tuple, type(np.array([]))]:
                # then list the element and save it
                groups[i] = list(groups[i])
            # Otherwsie,
            else:
                # Keep the element the same
                groups[i] = groups[i]
    # If the input is a numpy array
    elif type(groups) == type(np.array([])):
        # Then use list comprehension to format the groups list.
        # Assumes the groups have been transposed in this instance
        groups = [[groups[j][i] for j in range(len(groups))] for i in range(len(groups[0]))]
    # At then end, return the groups list.
    return groups

def storey_check_test(test):
    
    """
    =================================================================================================
    storey_check_test(test)
    
    This function is meant to check and format the T-test type from the inputs. The function works
    almost exactly like storey_check_groups()
    
    =================================================================================================
    Arguments:
    
    test  ->  A list, array, tuple, or array describing the T-test types used for each P-value.
    
    =================================================================================================
    Returns: A properly formatted list
    
    =================================================================================================
    """
    
    # If the groups are dataframes, make the input into a list  of two-lists
    if type(test) == type(pd.DataFrame()) or type(test) == type(pd.Series([])):
        # If the input is a series or DataFrame object
        # then attempt to list it
        try:
            test = list(test)
        except:
            raise ValueError("The test dataframe is incorrectly formatted. Try:  df_name['Test']")
    # If the input type is a list
    elif type(test) == list:
        # then iterate through each element of the list
        for i in range(len(test)):
            # and if any elements are not strings, then string them
            if type(test[i]) != str:
                test[i] = str(test[i])
            else:
                test[i] = test[i]
    # If the input type is a tuple
    elif type(test) == tuple:
        # then list the test
        test = list(test)
        # and loop over the elements of test
        for i in range(len(test)):
            # If any elements are not strings, then string them
            if type(test[i]) != str:
                test[i] = str(test[i])
            else:
                test[i] = test[i]
    # If the input is a numpy araray
    elif type(test) == type(np.array([])):
        # then use list comprehension to str all elements of the array
        test = [str(test[i]) for i in range(len(test))]
    # And at the end, return the test, reformatted
    return test

def storey_check_args(pvals,
                      groups,
                      test):
    """
    =================================================================================================
    storey_check_args(pvals, groups, test)
    
    This function is meant to check the arguments passed into the storey() function, and ensures that
    FDR estimation may proceed without conflict.
    
    =================================================================================================
    Arguments:
    
    pvals   ->  A list, numpy array, dataframe, tuple of P-values
    groups  ->  A list, numpy array, dataframe, tuple of group labels, index matched to the pvalues
    test    ->  A list, numpy array, dataframe, tuple of T-tests used for calculating P-values,
                index matched to the pvals argument.
    
    =================================================================================================
    Returns: The pvalues, g_checker boolean (group checker) and the t_checker boolean (test checker)
    
    =================================================================================================
    """
    # First, type-check the inputs 
    assert type(pvals) in [list, 
                           type(np.array([] ,dtype = float)),
                           type(pd.Series()),
                           type(pd.DataFrame())], "The p-values should be given as a list or a numpy array"
    assert type(groups) in [type(None),
                            list, 
                            tuple,
                            type(np.array([])),
                            type(pd.Series()),
                            type(pd.DataFrame())], "The p-values should be given as a list, tuple, numpy array, series or dataframe."
    # Then, if the pvals were a series or DataFrame object
    if type(pvals) == type(pd.Series()) or type(pvals) == type(pd.DataFrame()):
        # Turn them into numpy arrays and transpose them
        pvals = pvals.to_numpy().transpose()
        # Then, if the length of pvals is not 1, then raise an error
        if len(pvals) != 1:
            raise ValueError("The DataFrame or Series input has more than one dimension...")
        # Otherwise, pvals are the zeroeth element of the array
        else:
            pvals = pvals[0]
    # Next, check the groups. If somethign other than NoneType
    # was provided
    if type(groups) in [list, 
                        tuple,
                        type(np.array([])),
                        type(pd.Series()),
                        type(pd.Series()),
                        type(pd.DataFrame())]:
        # Then set g_checker to True, so we will check groups
        g_checker = True
    # Otherwise, set g_checker to False, as we do not need to check groups
    else:
        g_checker = False
    # If the test is a proper typed object
    if type(test) in [list, 
                      tuple,
                      type(np.array([])),
                      type(pd.Series()),
                      type(pd.Series()),
                      type(pd.DataFrame())]:
        # Then set t_checker to True, as we need to check the test
        t_checker = True
    # Otherwise, set t_checker to False, as we do not need to check test
    else:
        t_checker = False
    # and return pvals, g_chekcer and t_checker
    return pvals, g_checker, t_checker

def storey_make_id_dict(pvals, 
                        groups,
                        test,
                        g_checker,
                        t_checker):
    
    """
    =================================================================================================
    storey_make_id_dict(pvals, groups, test, g_checker, t_checker)
    
    This function is meant to take all relevant arguments to storey() and perform checking
    operations on each of those inputs.
    
    =================================================================================================
    Arguments:
    
    For more information on pvals, groups, test, refer to storey_check_args().
    
    g_checker  ->  A boolean, determines whether a group is in need of checking
    t_checker  ->  A boolean, determines whether a test is in need of checking
    
    =================================================================================================
    Returns: A dictionary based on the pvals argument, and the groups and test arguments checked.
    
    =================================================================================================
    """
    # Initialize the idenities dictionary
    identities = {}
    # Then proceed with making the identity dictionary
    # If groups are given and tests are not given
    if g_checker != False and t_checker == False:
        # Make sure the groups are in the correct format.
        # Otherwise, terminate gracefully.
        groups = storey_check_groups(groups)
        # If there are not enough group labels given for all the pvals,
        if len(groups) != len(pvals): 
            # Just proceed without group labels
            print("Each p-value should have a corresponding label/group tuple, proceeding without labels")
            # And make a dict of lists with key = pval, value = [position]
            for i in range(len(pvals)):
                identities[f"{i}?{pvals[i]}"] = [i] 
        # Otherwise, use the labels as the keys of the dictionary
        else:
            # make a dict of lists with key = pval, value = [position, label]
            for i in range(len(pvals)):
                identities[f"{i}?{pvals[i]}"] = [i, *groups[i]]
    # If no groups were provided but tests were provieded
    elif g_checker == False and t_checker != False:
        # Make sure the tests are in the right format. Otherwise, terminate gracefully
        test = storey_check_test(test)
        # If there are not enough tests given for all the pvals,
        if len(test) != len(pvals): 
            # Just proceed without labels
            print("Each p-value should have a corresponding test, proceeding without test identifier")
            # And make a dict of lists with key = pval, value = [position]
            for i in range(len(pvals)):
                identities[f"{i}?{pvals[i]}"] = [i]
        # Otherwise, use the tests as the keys of the dictionary
        else:
            # make a dict of lists with key = pval, value = [position, label]
            for i in range(len(pvals)):
                identities[f"{i}?{pvals[i]}"] = [i, test[i]]
    # If both tests and groups are provided
    elif g_checker != False and t_checker != False:
        # Make sure they're in the right format. Otherwise, terminate gracefully
        groups = storey_check_groups(groups)
        test = storey_check_test(test)
        # If there are not enough labels given for all the pvals,
        if len(groups) != len(pvals) and len(test) != len(pvals): 
            # Just proceed without labels
            print("Each p-value should have a corresponding label/group tuple and test, proceeding without labels and test identifiers")
            # And make a dict of lists with key = pval, value = [position]
            for i in range(len(pvals)):
                identities[f"{i}?{pvals[i]}"] = [i]
        # Otherwise, use the labels as the keys of the dictionary
        elif len(groups) != len(pvals) and len(test) == len(pvals):
            print("Each p-value should have a corresponding test, proceeding without test identifiers")
            for i in range(len(pvals)):
                identities[f"{i}?{pvals[i]}"] = [i, *groups[i]]
        elif len(groups) == len(pvals) and len(test) != len(pvals):
            print("Each p-value should have a corresponding label/group tuple, proceeding without labels")
            for i in range(len(pvals)):
                identities[f"{i}?{pvals[i]}"] = [i, test[i]]
        else:
            # make a dict of lists with key = pval, value = [position, label]
            for i in range(len(pvals)):
                identities[f"{i}?{pvals[i]}"] = [i, *groups[i], test[i]]
    # If no labels are given, then just make the identities dictionary
    else:
        # by looping over the pvals
        for i in range(len(pvals)):
            # and making keys as index/pval and value index.
            identities[f"{i}?{pvals[i]}"] = [i]
    # Once checking is over, return the identities dictionary, groups and test
    return identities, groups, test

def storey_reorder_qs(qs, 
                      pvals, 
                      og_pvals):
    """
    =================================================================================================
    storey_reorder_qs(qs, pvals, og_pvals)
    
    This function is used in storey(), and is meant to take the list of qvalues, the list of pvalues,
    and the original list of pvalues, and reorder the qvalue list in correspondence with the original
    pvalue list.
                  
    =================================================================================================
    Arguments:
    
    qs        ->  A list of q values created using the pvals list
    pvals     ->  A list of p values used to estimate the q values
    og_pvals  ->  A list of p values in their original order
    
    =================================================================================================
    Returns: A list of q values which ahve been reordered to match the order of og_pvals
    
    =================================================================================================
    """
    # Initialize the list of seen pvalues and new qvalues
    seen = []
    newqs = []
    # Loop over original order of pvalues
    for i in range(len(og_pvals)):
        # If the current pvalue is already seen
        if og_pvals[i] in seen:
            # Then find the the index of this particular pvalue. It will
            # find the first instance of this pvalue in the pvals list.
            ind = pvals.index(og_pvals[i])
            # Then, see how many of these pvalues have been identified
            num = seen.count(og_pvals[i])
            # and qvalue corresponding to this pvalue is the index of
            # the pvals list up to the current pval + the number seen,
            # plus the number of elements before that list.
            ind = pvals[ind+num:].index(og_pvals[i]) + len(pvals[:ind+num])
        # If the current pvalue is not seen
        else:
            # find the index of the og_pval[i] in pvals
            ind = pvals.index(og_pvals[i])
        # move qs[ind] to newqs
        newqs.append(qs[ind])
        # Add this value to the seen list
        seen.append(og_pvals[i])
    # Once the loop is complete, the q value will be reordered based
    # on the original pvalue order.
    return newqs
    
    

def pi_0(ps, lam):
    """
    =================================================================================================
    pi_0(ps, lam)
    
    This function is used to calculate the value of pi_0 in the Storey FDR estimation algorithm
                  
    =================================================================================================
    Arguments:
    
    ps   ->  A list of pvalues
    lam  ->  A critical value of lambda
    
    =================================================================================================
    Returns: The sum of all p values greater than lambda divided by the number of p values times the
             difference and 1 and lambda.
    
    =================================================================================================
    """
    # hat(pi_0) = num(p_j >0) / len(p)/(1-lambda)
    # This just uses numpy functions to do that for us
    return np.sum(ps>lam) / (ps.shape[0]*(1-lam))

def storey(pvals, 
           pi0 = None, 
           groups = None,
           test = None):
    """
    =================================================================================================
    storey(pvals, pi0, groups, test)
    
    This function performs Storey False Discovery Rate Estimation, as described in the publication
    Statistical Significance for Genomewide Studies (Storey, Tibshirani 2003)
    https://www.pnas.org/content/pnas/100/16/9440.full.pdf
    
    =================================================================================================
    Arguments:
    
    pvals   ->  A list of pvalues. Can be unordered
    pi0     ->  A value for pi0. This does not need to be set.
    groups  ->  A list, tuple, dataframe, or numpy array describing comparison groups
    test    ->  A list, tuple, dataframe or numpy array describing the test used for Pvalue generation
    
    =================================================================================================
    Returns: A DataFrame describing the P-values, the Q-values, and whatever metadata was provided.
    
    =================================================================================================
    """
    # First, get a list of the group names if groups are provided.
    group_names = None
    if type(groups) == type(pd.DataFrame([])):
        group_names = list(groups.columns.values)
    # and check the arguments provided to the function.
    og_pvals, g_checker, t_checker  = storey_check_args(pvals, groups, test)
    #
    ##################################################
    # 1.) The pvalues should be sorted. Make a dict to keep track of the order and info
    # Use the storey_make_id_dict() function to organize the arguments
    identities, groups, test = storey_make_id_dict(og_pvals,
                                                   groups,
                                                   test,
                                                   g_checker,
                                                   t_checker)
    # Once the original sequence is preserved, sort the pvals
    pvals = sorted([float(val) for val in og_pvals])
    og_pvals = [float(val) for val in og_pvals]
    # IF a list is given, make it a numpy array
    if type(pvals) == list:
        pvals = np.array(pvals, dtype = float)
    #
    ##################################################
    # 2,3,4.) Calculate the hat(pi_0) for the lambdas and Fit the cubic spline or just set pi0s value
    # Make a range of lambda values and make an empty array
    # to hold the hat(pi_0)s
    # The lambdas are meant to be a range of values used to fit the cubic spline
    lambdas = np.arange(0.01,0.97,0.01)
    # If pi0 is not set, and the pvalues are less than 100 in length
    if len(pvals) < 100 and pi0 is None:
        # Set the hat_pi0 to 1, as it will be clsoe to 1 anyway
        hat_pi0 = 1.0
    # If pi0 is set, then set it to the input value
    elif pi0 != None:
        hat_pi0 = pi0
    # Otherwise, determine an appriopriate value.
    else:
        # Make the pi0s array
        pi0s = np.empty(lambdas.shape[0])
        # and for each lambda
        for i, cur_l in enumerate(lambdas):
            # calculate pi0 and assign it to the pi0s array
            pi0s[i] = pi_0(pvals, cur_l)
        #spline = UnivariateSpline(x=lambdas, y=pi0s, k=3)
        #hat_pi0 = spline(1)
        # Then use the scipy function splrep and splrev to
        # fit the cubic spline and calculate hat_pi0
        tck = splrep(lambdas, pi0s, k=3)
        hat_pi0 = splev(lambdas[-1], tck)
        # If the value of hat_pi0 is greater than 1
        if hat_pi0 > 1:
            # Just set it to 1, it cannot be greater than 1
            hat_pi0 = 1
    # If the value of hat_pi0 is 1
    if hat_pi0 == 1:
        # Then this algorithm is equivalent to BH
        method = "Benj-Hoch"
    # If hat_pi0 is less than 1
    else:
        # Then this is the true Storey algorithm
        method = "Storey"
    #
    ##################################################
    # 5.) Make the calculation q(p_m) = min_{t\geq p_{i}} * \dfrac{(hat_pi0 * m * t)}{num(p_j \leq t)} = hat_pi0 * p_m
    # Initialize the qs list, by multiplying hat_pi0 to
    # the highest pvalue
    qs = [hat_pi0 * pvals[-1]]
    #
    ##################################################
    # 6.) Calculate min( \dfrac{hat_pi0 * m * p_i}{i}, q(p_{i+1}))
    # Then use the calculation above to determine what the next qvalue should be
    # Loop backwards over length of pvals minus 1 to zero, since we have qm already
    for i in range(len(pvals)-1, 0, -1):
        # append the minimum of the previous q and the next prediction.
        # Note pvals[i-1], since i is mathematically counted (1,2,...,m)
        qs.append(min(hat_pi0*len(pvals)*pvals[i-1]/i, qs[-1]))
    # This list is reversed compared to the pvals list, so reverse it
    qs = [qs[-i] for i in range(1,len(qs)+1)]
    #
    ##################################################
    # 7.) The q value estimate for the ith most significant feature is q(p_i)
    # So we just need to format the output. Since the identities has the information
    # to format the qvals list,
    # align the pvalues and qvalues
    qs = storey_reorder_qs(qs, list(pvals), og_pvals)
    # Now, we need only to format the output properly.
    # if groups were provided, but test was not provided
    if g_checker != False and t_checker == False:
        # and if group_names == None,
        if group_names == None:
            # Then use the generic group names as the labels
            g_nums = [f"Group {i+1}" for i in range(len(groups[0])) ]
        # otherwise
        else:
            # Use the group names found in the beginning
            g_nums = group_names
        # and format the output as a list of lists of columns names.
        output = [["index", *g_nums, "Test", "pvalue", "method", "qvalue"]]
    # Or if the groups are not provcided and the test is provided
    elif g_checker == False and t_checker != False:
        # format the output headers as follows
        output = [["index", "Test", "pvalue", "method", "qvalue"]]
    # Or if both groups and test are provided
    elif g_checker != False and t_checker != False:
        # then determine whether group names were provided
        if group_names == None:
            # and if not, use generic group names
            g_nums = [f"Group {i+1}" for i in range(len(groups[0]))]
        # Otherwise,
        else:
            # use the group names identified eearlier
            g_nums = group_names
        # and make the output headers as follows
        output = [["index", *g_nums, "Test", "pvalue", "method", "qvalue"]]
    # IF no group names and no test was provided
    else:
        # Make the output headers minimal.
        output = [["index", "pvalue", "method", "qvalue"]]
    # Loop over the identities dictionary
    for pval, value in identities.items():
        # and add lists to the output, containing the value, the pvalue, Storey, and the qvalue
        output.append([*value, float(pval.split("?")[1]), f"{method}", qs[value[0]]])
    # Turn the outputs list into a dataframe
    returner = pd.DataFrame(np.array(output[1:]), columns = output[0])
    # remvoe the index column
    del returner["index"]
    # and return the output dataframe.
    return returner
###
#
#
#####################################################################################################################
#
#   Statistical Tests

def t_check_args(*args,
                 nan_policy = "omit"):
    
    """
    =================================================================================================
    t_check_args(*args, nan_policy)
    
    This function is meant to check the input arguments into the ttest() function.
    
    =================================================================================================
    Arguments:
    
    args        ->  Two numpy arrays, list, tuples of data.
    nan_policy  ->  How to handle "not a number" values. "omit", "propogate", "raise"
    
    =================================================================================================
    Returns: A list of the arguments, where nan values have been handled properly.
    
    =================================================================================================
    """
    # Check to make sure the number of arguments are appropriate
    assert len(args) == 2, "A T-Test only takes two arguments"
    assert nan_policy.lower() in ["omit", "propogate", "raise"], "Your nan policy is not defined properly. nans = 'omit'/'propogate'/'raise' are the only accepted answers."
    # Once the arguments have been checked, continue.
    # Initialize a list for the reformatted argument lists
    reformed = []
    for arg in args:
        # If the user elects to omit nan values
        if nan_policy.lower() == "omit":
            # Then check to see if the argument is a list
            if type(arg) == list:
                # If so, NaNs do not equal themselves, so use
                # list comprehension to filter out items that do not
                # equal themselves
                arg = [item for item in arg if item == item]
                # and add this list, as a numpy array, to the reformed list.
                reformed.append(np.array(arg))   
            # If the argument is a numpy array
            elif type(arg) == type(np.array([])):
                # Then Filter out nans using logical not operator
                arg = arg[~np.isnan(arg)]
                # And determine whether this is a one dimensional, and numeric
                check1 = len(arg) == 1 and [type(arg[0][i]) in [np.int64, np.int32, np.float64, np.float32] for i in range(len(arg[0]))] == [True for _ in range(len(arg[0]))]
                # And determine whether every argument is a number
                check2 = all([type(arg[i]) in [np.int64, np.int32, np.float64, np.float32] for i in range(len(arg))])
                # If either check1 or check2 are True,
                if check1 or check2:
                    # Add the argument tot he reformed list.
                    reformed.append(arg)
                # Otherwise, raise a value error and cease.
                else:
                    raise ValueError("Something about the given numpy array is incorrectly formatted.")
        # If the nan_policy is either propogate or raise
        elif nan_policy.lower() == "propogate" or nan_policy.lower() == "raise":
            # Then, if the argument is a list
            if type(arg) == list:
                # NaNs do not equal themselves, so determine how many nans are in the list.
                nancheck = [1 for _ in arg if item != item]
                # If nans were found during this check, and the policy is to propogate
                if sum(nancheck) > 0 and nan_policy.lower() == "propogate":
                    # Then return two nan values.
                    return np.nan, np.nan
                # If nans were found tand the policy is to raise
                elif sum(nancheck) > 0 and nan_policy.lower() == "raise":
                    # Then raise a ValueError
                    raise ValueError("A NaN was detected in one of the inputs...")
                # If no NaNs were found, then
                else:
                    # Add the argument to the reformed list.
                    reformed.append(np.array(arg))   
            # If the argument is a numpy array
            elif type(arg) == type(np.array([])):
                # Then check to see if the numpy arrays are numbers and there is only one list in the array,.
                check1 = len(arg) == 1 and [type(arg[0][i]) in [np.int64, np.int32, np.float64, np.float32] for i in range(len(arg[0]))] == [True for _ in range(len(arg[0]))]
                check2 = all([type(arg[i]) in [np.int64, np.int32, np.float64, np.float32] for i in range(len(arg))])
                # If either check1 or check2 are true 
                if check1 or check2:
                    # Filter out nans using logical not operator
                    nancheck = arg[~np.isnan(arg)]
                    # If there is a nan in the nancheck array and the policy is propogate
                    if True in nancheck and nan_policy.lower() == "propogate":
                        # Then return nan values
                        return np.nan, np.nan
                    # If there is a nan in the nancheck array and the policy is raise
                    elif True in nancheck and nan_policy.lower() == "raise":
                        # Then raise an error
                        raise ValueError("A NaN was detected in one of the inputs...")
                    # Otherwise, add the argument to the reformed list
                    else:
                        reformed.append(arg)
                # If check1 or check2 are untrue, then raise an error.
                else:
                    raise ValueError("Something about the given numpy array is incorrectly formatted.")
    # If the loop has completed and neither raised an error nor return nan values,
    # then return the reformatted lists.
    return reformed

def t_statistic(data1, 
                data2, 
                var_term, 
                t_type = "welch"):
    
    """
    =================================================================================================
    t_statistic(data1, data2, var_term, t_type)
    
    This function is meant to calculate the t statistic of an independent two sample T-test by
    the following:
        Students T-test
    T = mean(data1) - mean(data2) / (var_pooled * sqrt(1/len(data1) + 1/len(data2)))
    
            where var_pooled is defined here:
                    https://en.wikipedia.org/wiki/Student%27s_t-test#Equal_sample_sizes_and_variance
    
        Welch's T-test
    T = mean(data1) - mean(data2) / (sqrt(var_1/len(data1)+ var_2/len(data2)))
    
    This function should only be used inside of the ttest function, although if the data1 and
    data2 are numpy arrays with no nan values, there is nothing stopping you using this elsewhere.
    
    =================================================================================================
    Arguments:
    
    data1     ->  A numpy array of data
    data2     ->  A numpy array of data
    var_term  ->  A float, the variance
    t_type    ->  A string determining the type of test, "w", "s"
    
    =================================================================================================
    Returns: The T-statistic defined above, for either the Welch's or Student's T-test
    
    =================================================================================================
    """
    # Check to make sure the test type is value
    assert t_type.lower() in ["student", "s", "welch", "w"], "The type of T-Statistic is not well defined.."
    # And the numerator is always the difference of the means
    numerator = data1.mean() - data2.mean()
    # If a students T-test is being used
    if t_type.lower() in ["student", "s"]:
        # then var_term is a pooled variance, so calculate the
        # denominator accordingly.
        denominator = var_term * np.sqrt( 1/len(data1) + 1/len(data2) )
    # Or if the type is welch's T-test, then
    elif t_type.lower() in ["welch", "w"]:
        # The var term should be the entire denominator.
        denominator = var_term
    # Once the numerator and denominator are calculated
    # return the ratio of numerator/denominator
    return numerator / denominator

def t_degree_freedom(data1,
                     data2,
                     t_type = "welch"):
    """
    =================================================================================================
    t_degree_freedom(data1, data2, t_type)
    
    This function is meant to calculate the degrees of freedom for the T-test, given the data
    and the type of test being performed.
    
    This function should only be used inside of the ttest() function, although if both data1 and
    data2 are numpy arrays, this could be used outside of it.
    
    =================================================================================================
    Arguments:
    
    data1   ->  A numpy array of Data
    data2   ->  A numpy array of Data
    t_type  ->  A string determining the type of test being used, "w", "s"
    
    =================================================================================================
    Returns: The degrees of freedom of the data.
    
    =================================================================================================
    """
    # Check to make sure the test type is valid.
    assert t_type.lower() in ["student", "s", "welch", "w"], "The type of T-Test is not well defined.."
    # If the test type is student
    if t_type.lower() in ["student", "s"]:
        # The degrees of freedom are the sum of the lengths of the data minus 2
        deg_free = len(data1) + len(data2) - 2
    # If the test type is Welch's
    elif t_type.lower() in ["welch", "w"]:
        # The formula for this is complicated, the nice mathematics is here:
        # https://en.wikipedia.org/wiki/Student%27s_t-test#Equal_or_unequal_sample_sizes,_unequal_variances_(sX1_%3E_2sX2_or_sX2_%3E_2sX1)
        numerator = ( (data1.var(ddof = 1)/len(data1)) + ( data2.var(ddof = 1)/len(data2) ) )**2
        d = ( (data1.var(ddof = 1)/len(data1))**2 / (len(data1)-1) ) + ( (data2.var(ddof = 1)/len(data2))**2 / (len(data2)-1) )
        deg_free = numerator/d
    # Once the degrees of freedom have been calculated, return them.
    return deg_free
    
def t_pval(t_stat,
           deg_free,
           tails = 2):
    
    """
    =================================================================================================
    t_pval(t_stat, def_free, tails)
    
    This function is meant to calculate a P-value given a T-statistic, the degrees of freedom,
    and the tails for the test.
    
    This function is meant to be used within ttest(), although it could be used on it's own
    
    =================================================================================================
    Arguments:
    
    t_stat    ->  A float, representing the T-statistic between two data lists
    deg_free  ->  A float, the degrees of freedom corresponding to two data lists
    tails     ->  An iteger determining the tails of the distribution. (2, 1, -1)
    
    =================================================================================================
    Returns: A pvalue, the probability of the T-statistic given a T-distribution and the tails of
             said distribution.
    
    =================================================================================================
    """
    # First, check to make sure the tails argument is valid.
    assert tails in [2,1,-1], "The tails is not valid... It should be 2 (both tails), 1 (left tail), -1 (right tail)"
    # If two tails are being used
    if tails == 2:
        # Then the p-value is the probability that T > T-statistic
        # and -T-statistic < T
        pval = 2 * (1 - ss.t.cdf(abs(t_stat), df = deg_free))
    # If the left hand tail is being used
    elif tails == 1:
        # Then the P-value is the probability that -T-statistic < T
        pval = 1 - ss.t.cdf(abs(t_stat), df = deg_free)
    # If the right hand tail is being used
    elif tails == -1:
        # Then the P-value is the probability that T > T-statistic
        pval = ss.t.cdf(t_stat, df = deg_free)
    # Once the pvalue has been calculated, return it.
    return pval

def student_t(data1, 
              data2, 
              tails = 2):
    
    """
    =================================================================================================
    student_t(data1, data2, tails)
    
    This function is meant to calculate a Student's T-test given two data lists and the number of
    tails to calculate P-values
    
    =================================================================================================
    Arguments:
    
    data1   ->  A numpy array of Data
    data2   ->  A numpy array of Data
    tails   ->  An iteger determining the tails of the distribution. (2, 1, -1)
    
    =================================================================================================
    Returns: A numpy array containing the T-statistic, Degrees of Freedom, P-value, "student"
    
    =================================================================================================
    """
    # Check to make sure the input data is correct
    assert type(data1) in [list, type(np.array([]))], "The input data should be a list or a numpy array [data1]"
    assert type(data2) in [list, type(np.array([]))], "The input data should be a list or a numpy array [data2]"
    assert tails in [-1,1,2], "Valid option for 'tails' is: 1 or 2 (integer)"
    # Get the sample variance of each dataset.
    var_1 = data1.var(ddof = 1)
    var_2 = data2.var(ddof = 1)
    # Calculate the pooled variance numerator and denominator
    pooled_var_num = (len(data1) - 1) * var_1 + (len(data2) - 1) * var_2
    pooled_var_den = len(data1) + len(data2) - 2
    # And wrap them together with square root
    pooled_var = np.sqrt(pooled_var_num/pooled_var_den)
    # Then calculate the T-statistic for the data
    t_stat = t_statistic(data1, data2, pooled_var, t_type = "student")
    # And the degrees of freedom for the data
    deg_free = t_degree_freedom(data1, data2, t_type = "student")
    # and finally the P-value for this comparison
    pval = t_pval(t_stat, deg_free, tails = tails)
    # Then return the numpy array of the results.
    return np.array([t_stat, deg_free, pval, "student"], dtype =object )

def welch_t(data1, 
            data2, 
            tails = 2):
    """
    =================================================================================================
    welch_t(data1, data2, tails)
    
    This function is meant to calculate a Welch's T-test given two data lists and the numebr of tails
    to calculate P-values with. 
    
    =================================================================================================
    Arguments:
    
    data1   ->  A numpy array of Data
    data2   ->  A numpy array of Data
    tails   ->  An iteger determining the tails of the distribution. (2, 1, -1)
    
    =================================================================================================
    Returns: A numpy array containing the T-statistic, Degrees of Freedom, P-value, "welch"
    
    =================================================================================================
    """
    # Check to make sure the input data is correct
    assert type(data1) in [type(np.array([]))], "The input data should be a list or a numpy array [data1]"
    assert type(data2) in [type(np.array([]))], "The input data should be a list or a numpy array [data2]"
    assert tails in [-1,1,2], "Valid option for 'tails' is: 1 or 2 (integer)"
    # Get the sample variance of the two input data lists
    var_1 = data1.var(ddof = 1)
    var_2 = data2.var(ddof = 1)
    # Then calculate the Welch's variance
    welch_var = np.sqrt( var_1/len(data1) + var_2/len(data2) )
    # Next, calculate the T-statistic for the data
    t_stat = t_statistic(data1, data2, welch_var, t_type = "welch")
    # and the degrees of freedom for the data
    deg_free = t_degree_freedom(data1, data2, t_type = "welch")
    # and finally, the pvalue for the data
    pval = t_pval(t_stat, deg_free, tails = tails)
    # Then return the numpy array containing the output
    return np.array([t_stat, deg_free, pval, "welch"], dtype = object )

def ttest(data1, 
          data2, 
          tails = 2, 
          t_type = "d",
          nan_policy = "omit"):
    """
    =================================================================================================
    ttest(data1, data2, tails, t_type, nan_policy)
    
    This function is meant to perform a T-test given two data lists and some information about
    how to perform the analysis.
    
    =================================================================================================
    Arguments:
    
    data1       ->  A numpy array of Data
    data2       ->  A numpy array of Data
    tails       ->  An iteger determining the tails of the distribution. (2, 1, -1)
    t_type      ->  A string determining the test to be used. "welch", "student", "determine"
    nan_policy  ->  A string determining how to handle nan values.
                    "omit", "propogate", "raise"
    
    =================================================================================================
    Returns: A numpy array desrcibing the statistics performed and values calculated.
    
    =================================================================================================
    """
    # Check to make sure the test type is valid.
    assert t_type.lower() in ["student", 
                              "s", 
                              "welch", 
                              "w",
                              "determine",
                              "d"], "The type of T-Test is not well defined.."
    # Then, check the input data for validity and handle nans
    data1, data2 = t_check_args(data1,
                                data2,
                                nan_policy = nan_policy)
    # If the nan_policy is to propogate, then a nan is returned when one is found
    if nan_policy.lower() == "propogate":
        if data1 != data1 and data2 != data2:
            return np.nan
    # If the test type is a Student's T-test
    if t_type.lower() in ["student", "s"]:
        # Then return the output of student_t()
        return student_t(data1, 
                         data2, 
                         tails = tails)
    # If the test type is a Welch's T-test
    elif t_type.lower() in ["welch", "w"]:
        # Then return the output of welch_t()
        return welch_t(data1, 
                       data2, 
                       tails = tails)
    # IF the test type is to determine
    elif t_type.lower() in ["determine", "d"]:
        # Then calculate the variance of each dataset
        var1 = data1.var(ddof = 1)
        var2 = data2.var(ddof = 1)
        # Student's T-test are more appropriate when
        # the variances of the datasets are close
        if 1/2 <= var1/var2 <= 2:
            # So use a Student's T-test
            return student_t(data1, 
                             data2, 
                             tails = tails)
        # A Welch's T-test is more appropriate when
        # the variances are not close to equal
        else:
            # So use a Welch's T-test
            return welch_t(data1,
                           data2,
                           tails = tails)
    # And raise a ValueError if something happens to make this explode.
    else:
        raise ValueError("An unknown error has occured.........")

def pairwise_t(*data_list, 
               comp_labels = True, 
               omit_comb = [],
               t_type = "d",
               nan_policy = "omit"):
    
    """ 
    =================================================================================================
    pairwise_t(*data_list, comp_labels, omit_comb, t_type, nan_policy)
    
    This function is meant to perform T-tests between the the input data.
                  
    =================================================================================================
    Arguments:
    
    data_list    ->  A list of lists, where each list is a set of data.
    comp_labels  ->  A boolean describing whether or not to include comparison labels.
                        If True, then the data_list is expected to be tuples, where element
                        zero is a label, and element one is a list of data.
    omit_comb    ->  A list of lists that describe combinations to omit from comparisons.
    t_type       ->  A string describing the type of comparison to be used.
    nan_policy   ->  A string describing how to handle nan values.
    
    =================================================================================================
    Returns: A Pandas DataFrame containing the results of the T-tests.
    
    =================================================================================================
    """
    # Check to make sure the inputs are correct
    assert comp_labels in [True, False], "comp_labels argument must be a boolean"
    assert type(omit_comb) == list, "Combinations to be ommitted should be a list of 2-tuples/2-lists"
    # Next, check the omit_comb list. If it is not empty
    if omit_comb != []:
        # Initialize the seen list
        seen = []
        # and loop over each combination to omit
        for comb in omit_comb:
            # and check each combination for validity
            if len(comb) == 1:
                assert type(comb) in [list, tuple] and len(comb) == 2, "Combinations to be ommitted should be a list of 2-tuples/2-lists"
            # Once the combination is checked, add it to the seen list
            seen.append(list(comb))
    # If the omit_comb list is empty
    else:
        # Then initialize an empty seen list.
        seen = []
    # After checking the omittion combinations, initialize
    # the statistics dictioanry
    stats = []
    # And loop over the given data lists
    for data1 in data_list:
        # And loop over the data lists again
        for data2 in data_list:
            # If the current combination has not been seen
            # and the lists are not the same list
            # Then also check whether or not comparison labels are expect. If not
            if [data1,data2] not in seen and [data2,data1] not in seen and data1 != data2 and comp_labels == False:
                # And add the forward and reverse lsits to the seen list
                seen.append([data1,data2])
                seen.append([data2,data1])
                # Then run a T-Test on the data and add it
                # to the stats list
                stats.append(list(ttest(data1,
                                        data2,
                                        t_type = t_type,
                                        nan_policy = nan_policy)))
            # If comparison labels are provided
            elif [data1,data2] not in seen and [data2,data1] not in seen and data1 != data2 and comp_labels == True:
                # And add the forward and reverse lsits to the seen list
                seen.append([data1,data2])
                seen.append([data2,data1])
                # Then run a T-Test on the data and add it
                # to the stats list with the labels
                stats.append([data1[0], 
                             data2[0],
                             *list(ttest(data1[1],
                                         data2[1],
                                         t_type = t_type,
                                        nan_policy = nan_policy))])
    # If comparison labels are expected,
    if comp_labels:
        #  then add the Group columns to the output DataFrame
        return pd.DataFrame(np.array(stats), columns = ["Group 1", "Group 2", "T-statistic", "DoF", "pvalue", "Test"])
    # If comparison labels are not expected
    else:
        # Then make the simple column headers for the DataFrame
        return pd.DataFrame(np.array(stats), columns = ["T-statistic", "DoF", "pvalue", "Test"])
    
    
### Tukey HSD
    
def tukey_hsd(group_names, 
              *args ):
    """
    =================================================================================================
    tukey_hsd(group_names, *args)
    
    This function is meant to perform Tukeys Range Test (Tukey's Honest Significant Difference) on
    the input arguments. This procedyre corrects for Family Wise Error Rates.
    
    Source code:
    https://s-nako.work/2020/01/one-way-anovaanalysis-of-variance-and-multiple-comparisons-in-python/
    =================================================================================================
    Arguments:
    
    group_names  ->  The list of names identifying which comparisons are happening
    args         ->  The data as numpy arrays
    
    =================================================================================================
    Returns: The results of Tukey HSD test from SciPy.
    
    =================================================================================================
    NOTES:
    
    Eventually, I would like to code this myself, however I do not currently have the time. I find
    that using some of the built in things require strange typings/functions.
    
    =================================================================================================
    """
    # Format the args (data) for use in pairwise_tukeyhsd
    endog = np.hstack(args)
    # Initialize the groups list
    groups_list = []
    # Loop over the args list
    for i in range(len(args)):
        # Loop over the elements in the args list.
        for j in range(len(args[i])):
            # Add the group name for the i,j comparison to the groups names list.
            groups_list.append(group_names[i])
    # Turn the groups list into an array
    groups = np.array(groups_list)
    # Run pairwise_tukeyhsd
    res = pairwise_tukeyhsd(endog, groups)
    # Return the object
    return res

#
#
#####################################################################################################################
#
#     Miscellaneous Functions

def filter_nans(data,
                threshold = 3, 
                threshold_type = "data"):
    """
    =================================================================================================
    filter_nans(data, threshold, threshold_type)
    
    This function is meant to filter out the nan values from a list, based on the input arguments.
                  
    =================================================================================================
    Arguments:
    
    data            ->  A list (or iterable) of data points. The points are assumed to be numbers.
    threshold       ->  An integer describing the minimum value requirement.
    threshold_type  ->  A string describing how the threshold integer will be applied. 
                        "on_data" "on_nan"
    
    =================================================================================================
    Returns: The filtered list, or an empty list if the threshold requirements were not met.
    
    =================================================================================================
    """
    # Make sure the user gave a valid thresholding option
    assert threshold_type.lower() in ["data", 
                                      "on_data", 
                                      "on data", 
                                      "nan", 
                                      "on_nan", 
                                      "on nan"], "Threshold is either relative to NaN or data."
    assert type(data) == list, "The data should be in a list"
    # Filter NaNs, as they do not equal themselves
    filtered = [val for val in data if val == val]
    # Keep data if there are at least <threshold> data points
    if threshold_type.lower() in ["data", "on_data", "on data"]:
        if len(filtered) >= threshold:
            return filtered
        else:
            return []
    # Keep data if there are no more than <threshold> nans
    elif threshold_type.lower() in ["nan", "on_nan", "on nan"]:
        if len(data) - len(filtered) <= threshold:
            return filtered
        else:
            return []
        
    
def filter_nan_dict(data,
                    threshold = 3, 
                    threshold_type = "data"):
    """
    =================================================================================================
    filter_nan_dict(data, threshold, thershold_type)
    
    This function is meant to filter out nan values from the list-values in a dictionary. This
    function uses filter_nans() to filter lists.
                  
    =================================================================================================
    Arguments:
    
    data            ->  A dictionary of lists of data points. The points are assumed to be numbers.
    threshold       ->  An integer describing the minimum value requirement.
    threshold_type  ->  A string describing how the threshold integer will be applied. 
                        "on_data" "on_nan"
    
    =================================================================================================
    Returns: A dictionary where all values have been filtered from the list-values.
    
    =================================================================================================
    """
    # Make sure a dictionary is given as input
    assert type(data) == dict, "The data should be in a dictionary"
    # Initialize the new dictionary
    filtered_dict = {}
    # Loop over the keys/values in the dictionary
    for key, value in data.items():
        # Filter the nans
        filt_list = filter_nans(value)
        # If the list is not empty
        if filt_list != []:
            # then add it to the dictionary
            filtered_dict[key] = filt_list
        # IF the list is empty, it will not be added to the dictionary
    # Return the filtered dictionary.
    return filtered_dict

def mean(dataset):
    
    '''
    Given a dataset, return the average value.
    '''
    
    return sum(dataset) / len(dataset)

def median(dataset):
    
    '''
    Given a dataset, return the median value.
    '''
    
    # The calculation for median depends on whether the 
    # number of datapoints in the set is even or odd. 
    # If the number of datapoints is even, then we need to 
    # find the average of the two middle numbers. If the 
    # number of datapoints is odd, then we simply need
    # to find the middle one. They also need to be sorted.
    
    dataset = sorted(dataset)
    if len(dataset) % 2 == 0:                            # if the dataset is even
        index = len(dataset) // 2                        # get the middle data point
        med = (dataset[index] + dataset[index -1]) / 2   # average the middle two points
        return med                                       # return this value
    elif len(dataset) % 2 == 1:                          # if the dataset is odd
        index = len(dataset) // 2                        # get the middle point
        return dataset[index]                            # return the middle point

def count_list(dataset): 
    
    '''
    Given a dataset, count the occurence of
    each data point and return them as a dictionary. 
    '''
    
    # First we create a dictionary to hold the counts. 
    # We then loop over the elements of the dataset
    # and attempt to check the dictionary keys for them.
    # If the element has appeard, we add one to the count
    # and if the element is not in the dictionary, we 
    # add a key to the dictionary with that elements name 
    # and initialize it to 1 (since we have seen it once). 
    # At the end, return the dictionary. 
    
    dic = {}                          # Create the empty dictionary
    for element in dataset:           # Loop over the elemnts in the dataset
        try:                          # Attempt
            dic[str(element)] +=  1   # Look for the key of the element, add one to count
        except:                       # Otherwise 
            dic[str(element)] = 1     # Add a key to the dicitonary with value 1
    return dic                        # Return the dictionary

def mode(dataset):
    
    '''
    Given a dataset, returns the most frequent value
    and how many times that value appears
    '''
    
    # First, we count all of the elements and arrange them in 
    # a dictionary. Then, we create a sorted list of tuples 
    # from the key value pairs. Initialize 'pair', to hold
    # the key value pair of the highest value, and an empty
    # list to hold any of the pairs that tie. We then loop
    # over the sorted lists keys and values, and look for the 
    # highest counts. We return the highest count from the 
    # dictionary, or the highest counts if there were any ties.
    
    counted = count_list(dataset)      # Count the elements of the dataset
    sort = sorted(counted.items())     # Sort the numbers and occurences
    pair = 'hold', 0                   # Initialize the pair
    ties = []                          # Initialize the tie list
    for key, value in sort:            # Loop over key, value in sorted dictionary
        if value > pair[1]:            # If the value is greater than the pair
            pair = key, value          # Re assign the pair to the current one
            ties = []                  # Reset the tie list
        elif value == pair[1]:         # If the value is equal to the current value
            ties.append((key, value))  # Append the new key, value pair to the list
    ties.append(pair)                  # After, append the pair to the list
    svar = sorted(ties)                # Sort the list of ties
    if len(ties) > 1:                  # If there are any ties, 
        return svar                    # Return the sorted list of ties 
    elif len(ties) == 1:               # If there are no ties
        return pair                    # Return the highest value

def quantile(dataset, percentage):
    
    '''
    Given a dataset and a pecentage, the function returns
    the value under which the given percentage of the data
    lies
    '''
    
    # First, sort the dataset, then find the index at the 
    # given percentage of the list. Then, return the 
    # value of the dataset at that index. 
     
    dataset = sorted(dataset)               # Sort the dataset
    index = int(percentage * len(dataset))  # Get the index of element percentage 
    return dataset[index]                   # return the element at the index  

def interquantile_range(dataset, per_1, per_2):
    
    '''
    Given a dataset and two percentages that define 
    a range of the dataset, find the range of the 
    elements between those elements. 
    '''
    dataset = sorted(dataset)
    return quantile(dataset, per_2) - quantile(dataset, per_1)  

def data_range(dataset):
    
    '''
    Given a dataset, return the range of the elements.
    '''
    
    dataset = sorted(dataset)
    return dataset[-1] - dataset[1]

def remove_mean(dataset):
    
    '''
    Given a dataset, return the dataset subtracting
    the mean from every element.
    '''
    
    m = mean(dataset)
    return [x-m for x in dataset]

def variance(dataset): 
    
    '''
    Given a dataset, calculate the variance of the 
    parent population of the dataset. 
    '''
    
    # Calculate the data without the mean, square
    # all of the elements, then return the sum of 
    # those squares divided by the number of 
    # datapoints minus 1. 
    
    meanless = remove_mean(dataset)          # Remove the mean from the data
    squares = [x**2 for x in meanless]       # Square all of the meanless datapoints
    return sum(squares) / (len(dataset) - 1) # return the sum of the squares divided by n-1

def standard_deviation(dataset):
    
    '''
    Given a dataset, return the standard_deviation.
    '''
    
    v = variance(dataset)
    return v**0.5

def dot_product(data_1, data_2):
    
    '''
    Given two datasets of equal length, return the 
    dot product. 
    '''
    
    # First, we make sure that the lists are the same size,
    # Then we loop over the length of the lists, and sum the 
    # product of the corresponding elements of each list. 
    # Then, that sum is returned. 
    
    assert len(data_1) == len(data_2), "These lists are not the same length"
    
    sum_total = 0                            # Initialize the sum
    for i in range(len(data_1)):             # Loop over the size of the list
        sum_total += data_1[i] * data_2[i]   # Add to the sum the product of the datapoints in 1 and 2
    return sum_total                         # Return the sum

def covariance(data_1, data_2):
    
    '''
    Given two datasets, calculate the covariance between them
    '''
    
    n = len(data_1)
    return dot_product(remove_mean(data_1),remove_mean(data_2)) / (n-1)

def correlation(data_1, data_2):
    
    '''
    Given two datasets, calculate the correlation between them. 
    '''
    
    return covariance(data_1, data_2) / (standard_deviation(data_1) * standard_deviation(data_2))

def vector_sum(vectors):
    
    '''
    Given a set of vectors, return a vector which contains the 
    sum of the ith elements from each vector in index i
    '''
    
    for i in range(len(vectors)-1):
        assert len(vectors[i]) == len(vectors[i+1]), 'Vectors are not the same length'
    
    return [sum(vector[i] for vector in vectors)
           for i in range(len(vectors[0]))]

assert vector_sum([[1,2],[2,3],[3,4]]) == [6,9]

def scalar_multiply(scalar, vector):
    
    '''
    Given a scalar and a vector (list), return a vector
    where each component is multiplied by the scalar. 
    '''
    
    return [scalar * var for var in vector]

assert scalar_multiply(3, [1,2,3,4]) == [3,6,9,12]

def vector_subtract(vectors):
    
    '''
    Given a set of vectors, return the difference between
    the vectors, in index order. 
    
    This will look like:
    vectors[0] - vectors[1] - ... - vectors[n] = result
    '''
    
    for i in range(len(vectors)-1):
        assert len(vectors[i]) == len(vectors[i+1]), 'Vectors are not the same length'
    
    pass_count = 0
    result = vectors[0]
    for column in vectors:
        if column == result and pass_count == 0:
            pass_count += 1
            pass
        else:
            for i in range(len(result)):
                result[i] += -column[i]
            pass_count += 1
    return result      

assert vector_subtract([[1,2,3], [3,4,5]]) == [-2,-2,-2]

def vector_mean(vectors):
    
    '''
    Given a list of lists (which correspond to vectors, where
    each element of the vector represents a different variable)
    return the vector mean of the vectors (add each the vectors
    component-wise, divide each sum by the number of vectors)
    '''
    
    n = len(vectors)
    return scalar_multiply(1/n, vector_sum(vectors))

assert vector_mean([[1,2],[2,3],[3,4]]) == [2,3]

def scale(data_1):
    
    '''
    Given a set of datapoint sets, return the mean of each
    dataset and the standard deviation for each set. 
    
    Data points should be give as 
    x = [[x1_1, x2_1,..., xn_1],...,[x1_n, x2_n,..., xn_n]]
    if the data are not in this format, but are in the format
    x = [[x1_1, x1_2,..., x1_n],...,[xn_1, xn_2,..., xn_n]]
    apply the function reformat_starts(*args) to the data. 
    '''
    
    # First, we make sure that all of the data points
    # given are the same length, then save the size of
    # the datasets as n. We then calculate the vector 
    # mean of the data, as well as the standard deviations
    # of each data type. Then the means and SDs are returned
    
    for q in range(len(data_1) -1):
        assert len(data_1[q]) == len(data_1[q+1]), 'Data lists are different sizes'
    n = len(data_1[0])
    
    means = vector_mean(data_1)
    s_deviations = [standard_deviation([vector[i] 
                                        for vector in data_1])
                   for i in range(n)]
    return means, s_deviations

t_vectors = [[-3, -1, 1], [-1, 0, 1], [1, 1, 1]]
t_means, t_stdevs = scale(t_vectors)
assert t_means == [-1, 0, 1]
assert t_stdevs == [2, 1, 0]

def rescale(data_1):
    
    '''
    Given a set of data sets, return a list of the
    data rescaled, based on the means and the 
    standard deviations of the data. 
    '''
    
    # First, we calculate the mean and standard deviations
    # of the data, and save the size of the datasets as 
    # n. We then copy each of the vectors to the list
    # rescaled. Next, we loop over the vectors in rescaled,
    # and loop over the size of the datasets, and 
    # scale each term in v[i] based on the mean and SD
    
    means, s_deviations = scale(data_1)
    n = len(data_1[0])
    rescaled = [v[:] for v in data_1]
    for v in rescaled:
        for i in range(n):
            if s_deviations[i] > 0:
                v[i] = (v[i] - means[i]) / s_deviations[i]
    return rescaled

t2_means, t2_stdevs = scale(rescale(t_vectors))
assert t2_means == [0, 0, 1]
assert t2_stdevs == [1, 1, 0]

def unscaled(scaled_data_1, data_1, coefficients = False):
    
    '''
    Given a set of scaled datapoints, the original datapoints,
    and a truthy value for whether we are unscaling coefficients
    of regression, return the unscaled data points. 
    '''
    
    # This is basically 'rescale' in reverse, with the
    # condition of if we are unscaling coefficients. If 
    # we are unscaling coefficients, we subtract from the 
    # alpha term (v[0]) all elements in the form 
    # v[j] * mean[j] / s_deviations[j] (as described in 
    # Data Science from Scratch, 2nd Edition in Chapter 
    # 16, Logistic Regression). OTherwise, all coefficient 
    # are divided by the standard deviation term of the 
    # corresponding data. 
    
    n = len(data_1[0])
    means, s_deviations = scale(data_1)
    unscaled = [v[:] for v in scaled_data_1]
    for v in unscaled:
        for i in range(n):
            if coefficients == False:
                v[i] = v[i]*s_deviations[i] + means[i]
            elif coefficients == True:
                if i == 0:
                    for j in range(1,n):
                        if s_deviations[j] > 0:
                            v[0] = v[0] - (v[j]*means[j])/s_deviations[j]
                        else:
                            v[0] = v[0] - v[j]
                elif i != 0: 
                    if s_deviations[i] > 0: 
                        v[i] = v[i] / s_deviations[i]
                    else:
                        pass
    return unscaled

assert unscaled(rescale(t_vectors), t_vectors) == t_vectors

#
#
#####################################################################################################################