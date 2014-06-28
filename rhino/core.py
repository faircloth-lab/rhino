
import os
import sys
import json
import time
import numpy
import dendropy
from collections import defaultdict

import pdb


def parse_site_rates(rate_file, correction = 1, test = False, count = 0):
    """Parse the site rate file returned from hyphy to a vector of rates"""
    # for whatever reason, when run in a virtualenv (and perhaps in other
    # cases, the file does not seem to be written quite before we try
    # to read it.  so, pause and try to re-read up to three-times.
    try:
        data = json.load(open(rate_file, 'r'))
    except IOError as e:
        if count <= 3:
            count += 1
            time.sleep(0.1)
            parse_site_rates(rate_file, correction, test, count)
        else:
            raise IOError("Cannot open {0}: {1}".format(rate_file, e))
    rates = numpy.array([line["rate"] for line in data["sites"]["rates"]])
    corrected = rates/correction
    if not test:
        data["sites"]["corrected_rates"] = [{"site":k + 1,"rate":v} \
                for k,v in enumerate(corrected)]
        json.dump(data, open(rate_file,'w'), indent = 4)
    return corrected


def correct_branch_lengths(tree_file, format, output_dir):
    """Scale branch lengths to values shorter than 100"""
    tree = dendropy.Tree.get_from_path(tree_file, format)
    depth = tree.seed_node.distance_from_tip()
    mean_branch_length = tree.length()/(2 * len(tree.leaf_nodes()) - 3)
    string_len = len(str(int(mean_branch_length + 0.5)))
    if string_len > 1:
        correction_factor = 10 ** string_len
    else:
        correction_factor = 1
    for edge in tree.preorder_edge_iter():
        if edge.length:
            edge.length /= correction_factor
    pth = os.path.join(output_dir, '{0}.corrected.newick'.format(
        os.path.basename(tree_file)
        ))
    tree.write_to_path(pth, 'newick')
    return pth, correction_factor


def get_net_pi_for_periods(pi, times):
    """Sum across the PI values for the requested times"""
    sums = numpy.nansum(pi, axis=1)[times]
    return dict(zip(times, sums))


def get_informative_sites(alignment, threshold=4):
    """Returns a list, where True indicates a site which was over the threshold
    for informativeness.
    """
    taxa = dendropy.DnaCharacterMatrix.get_from_path(alignment, 'nexus')
    results = defaultdict(int)
    for cells in taxa.vectors():
        assert len(cells) == taxa.vector_size # should all have equal lengths
        for idx, cell in enumerate(cells):
            results[idx] += 1 if str(cell).upper() in "ATGC" else 0
    return numpy.array([1 if results[x] >= threshold else numpy.nan for x in sorted(results)])


def cull_uninformative_rates(rates, inform):
    """Zeroes out rates which are uninformative"""
    return rates * inform
