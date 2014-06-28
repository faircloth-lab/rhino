#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
(c) 2014 Brant Faircloth || http://faircloth-lab.org/
All rights reserved.

This code is distributed under a 3-clause BSD license. Please see
LICENSE.txt for more information.

Created on 23 June 2014 17:37 PDT (-0700)
"""

import tempfile

import os
from rhino import io
from rhino import core
from rhino import hyphy
from rhino import parsimonator
from rhino.log import setup_logging

import pdb


def tree_worker(work):
    args, tmpdir, aln = work
    if args.tree_method == "parsimony":
        tree = parsimonator.get_parsimony_trees(tmpdir, aln)
        corrected_tree, correction = core.correct_branch_lengths(
            tree,
            "newick",
            tmpdir
        )
    return aln, corrected_tree, correction


def main(args):
    # print message
    #print welcome_message()
    # setup logging
    log, my_name = setup_logging(args)
    log.info("Getting alignments")
    # get aligment info
    alignments = io.get_alignment_files(args.alignments, args.input_format)
    # setup mp object as standin for map()
    if args.cores > 1:
        from multiprocessing import Pool, cpu_count
        pool = Pool(args.cores)
        mp = pool.map
    else:
        mp = map
    # create a base tmp dir to hold all files
    tmp_base = tempfile.mkdtemp()
    # create a temp dir to hold tree results
    tmp_tree = os.path.join(tmp_base, "trees")
    os.mkdir(tmp_tree)
    if args.tree_file:
        log.info("Correcting input tree")
        corrected_tree, correction = core.correct_branch_lengths(
            args.tree_file,
            "newick",
            tmp_tree
        )
        trees = [(aln, corrected_tree, correction) for aln in alignments]
    # infer a tree
    elif args.tree_method == 'parsimony':
        log.info("Generating parsimony tree(s)")
        # package the data
        work1 = [(args, tmp_tree, aln) for aln in alignments]
        trees = mp(tree_worker, work1)
    # create a temp dir to hold site rate results
    tmp_rates = os.path.join(tmp_base, "rates")
    os.mkdir(tmp_rates)
    # get hyphy rate template file
    hyphy_template = io.get_hyphy_conf()
    # package that data
    work2 = [(args, hyphy_template, tmp_rates, tree_data) for tree_data in trees]
    rates = mp(hyphy.rate_worker, work2)

    '''
    elif args.tree_file:
        pass
    # correct branch lengths

    # generate a vector of times given start and stops
    time_vector = tapir.get_time(0, int(tree_depth))
    params = []
    # get path to batch/template file for hyphy
    if not args.template:
        template = tapir.get_hyphy_conf()
    else:
        template = args.template
    if not args.site_rates:
        print "\nEstimating site rates and PI for files:"
        for alignment in tapir.get_files(args.alignments, '*.nex,*.nexus'):
            output = os.path.join(args.output, os.path.basename(alignment) + '.rates')
            towrite = "\n".join([alignment, tree, output])
            params.append([time_vector, args.hyphy, template, towrite, output, correction, alignment,
                args.times, args.intervals, args.threshold])
    else:
        print "Estimating PI for files (--site-rate option):"
        for rate_file in tapir.get_files(args.alignments, '*.rates'):
            params.append([time_vector, args.hyphy, template, None, rate_file,
                correction, rate_file, args.times, args.intervals,
                args.threshold])
    if not args.multiprocessing:
        pis = map(worker, params)
    else:
        from multiprocessing import Pool, cpu_count
        pool = Pool(processes = cpu_count() - 1)
        pis = pool.map(worker, params)
    # store results somewhere
    db_name = os.path.join(args.output,
        'phylogenetic-informativeness.sqlite')
    sys.stdout.write("\nStoring results in {0}...".format(db_name))
    sys.stdout.flush()
    conn, c = tapir.create_probe_db(db_name)
    tapir.insert_pi_data(conn, c, pis)
    conn.commit()
    sys.stdout.write("DONE")
    sys.stdout.flush()
    print "\n"
    c.close()
    conn.close()

    '''

    # ----------------------
    # end
    text = " Completed {} ".format(my_name)
    log.info(text.center(65, "="))
