###############################################################################
# import modules:

import sys, platform, os
import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy import integrate
import copy

###############################################################################
# import CAMB:

from . import model, initialpower
from .baseconfig import CAMBError

###############################################################################
# utility function to run parameter check:


def check_params(eftcamb_params, feedback=0):
    """
    Function that initializes EFTCAMB and returns the flags that are not included in the input parameters.
    """
    # try params:
    eftcamb_params['feedback_level'] = max(feedback, 0)
    try:
        pars = model.CAMBparams()
        pars.set_cosmology(H0=67)
        pars.EFTCAMB.initialize_parameters(pars, eftcamb_params, print_header=False)
        read_par = pars.EFTCAMB.read_parameters()
        name = pars.EFTCAMB.model_name()
    except Exception as ex:
        return None, eftcamb_params, None
    # get input parameters that are actually used:
    input_dict = {key: read_par[key] for key in read_par.keys() if (key in eftcamb_params.keys()) and (read_par[key] == eftcamb_params[key])}
    # get output parameters, i.e. parameteres that the code looked for and set to default value.
    output_dict = {key: read_par[key] for key in read_par.keys() if (key not in eftcamb_params.keys()) or (key in eftcamb_params.keys() and read_par[key] != eftcamb_params[key])}
    #
    return name, input_dict, output_dict

###############################################################################
# tree data structure:


class eftcamb_tree(object):

    def __init__(self, depth=0, name=None, in_params={}, out_params={}, children=None, changed_param=None):
        """
        Helper class to hold the EFTCAMB tree
        """
        # store data:
        self.depth = depth
        self.name = name
        self.in_params = in_params.copy()
        self.out_params = out_params.copy()
        self.changed_param = copy.deepcopy(changed_param)
        # initialize empty containers:
        self.children_common_params = []
        self.numerical_params = {}
        self.model_params = {}
        # add child trees:
        self.children = []
        if children is not None:
            for child in children:
                self.add_child(child)

    def add_child(self, node):
        """
        Add a children to the tree
        """
        assert isinstance(node, eftcamb_tree)
        self.children.append(node)


def get_all_keys(all_keys, tree):
    """
    Get all unique keys in the tree
    """
    # add all keys:
    for key in tree.in_params.keys():
        all_keys.append(key)
    for key in tree.out_params.keys():
        all_keys.append(key)
    # call on childs:
    for temp_tree in tree.children:
        get_all_keys(all_keys, temp_tree)
    # get unique elements:
    all_keys = list(set(all_keys))
    #
    return all_keys


def get_max_depth(tree):
    """
    Get the tree depth
    """
    depth = -1
    for child in tree.children:
        depth = max(depth, get_max_depth(child))
    return depth + 1

###############################################################################
# tree generator:


def generate_eft_tree(eftcamb_params, tree, feedback=0, depth=0, max_depth=1, max_int=20):
    """
    Generate the EFTCAMB model tree
    """
    # check for depth:
    if depth > max_depth:
        return None
    depth += 1
    # call parameter check:
    name, input_params, diff_params = check_params(eftcamb_params, feedback=feedback-1)
    # feedback:
    if feedback > 0:
        print('*************************************')
        print('Tree depth:', depth - 1)
        print('Input parameters:', input_params)
        print('Output parameter:', diff_params)
        print('*************************************')
    # branch depending on wether the model succeeded or not:
    if diff_params is None:
        return None
    # if the model succeeded then we save the successful parameters:
    tree.name = name
    tree.depth = depth - 1
    tree.in_params = input_params.copy()
    tree.out_params = diff_params.copy()
    # remove float params from diff_params:
    real_keys = [key for key in diff_params.keys() if type(diff_params[key]) == float]
    temp_diff_params = copy.deepcopy(diff_params)
    for key in real_keys:
        del(temp_diff_params[key])
    # stop iterations if there is not other parameters (ignoring floats):
    if len(temp_diff_params) == 0:
        return {}
    # loop over options:
    for key in temp_diff_params:
        # we got an int parameter, loop over integers:
        if type(temp_diff_params[key]) == int:
            temp_params = {**input_params.copy(), **temp_diff_params.copy()}
            _counter = temp_params[key] - 1
            while _counter < max_int:
                temp_params[key] = _counter + 1
                temp_tree = eftcamb_tree(changed_param=key)
                _temp = generate_eft_tree(temp_params, temp_tree, feedback=feedback, depth=depth, max_depth=max_depth, max_int=max_int)
                # if result is none calculation failed and we assume it's going to fail for all future ones:
                if _temp is None:
                    break
                # if initialization succeeded then we add the child and continue iterating:
                tree.add_child(temp_tree)
                _counter += 1
        # we got a logical, flip it:
        elif type(temp_diff_params[key]) == bool:
            temp_params = {**input_params.copy(), **temp_diff_params.copy()}
            temp_params[key] = not temp_params[key]
            temp_tree = eftcamb_tree(changed_param=key)
            _temp = generate_eft_tree(temp_params, temp_tree, feedback=feedback, depth=depth, max_depth=max_depth, max_int=max_int)
            if _temp is not None:
                tree.add_child(temp_tree)
    #
    return diff_params


def consolidate_tree_properties(tree):
    """
    Consolidate common properties
    """
    # get all keys in the child's tree:
    def _temp(all_keys, tree):
        for child in tree.children:
            all_keys.append(copy.deepcopy(list(child.in_params.keys()) + list(child.out_params.keys())))
            _temp(all_keys, child)
    all_keys = []
    _temp(all_keys, tree)
    # look for the keys that appear everywhere:
    if len(all_keys) > 0:
        _set = set(all_keys[0])
        for _keys in all_keys:
            _set.intersection(set(_keys))
        tree.children_common_params = list(_set)
    # get numerical and model selection parameters:
    temp_params = {**tree.in_params.copy(), **tree.out_params.copy()}
    real_keys = [key for key in temp_params.keys() if type(temp_params[key]) == float]
    tree.numerical_params = {key: temp_params[key] for key in real_keys}
    tree.model_params = {key: temp_params[key] for key in temp_params if key not in real_keys}
    # if the node has childrens iterate down:
    if len(tree.children) > 0:
        for child in tree.children:
            consolidate_tree_properties(child)

###############################################################################
# helper to build the whole thing:


def eftcamb_tree_helper(**kwargs):
    """
    Compute the full EFTCAMB tree
    """
    # initialize empty tree:
    _tree = eftcamb_tree()
    # compute the tree:
    generate_eft_tree(eftcamb_params={}, tree=_tree, **kwargs)
    #
    return _tree
