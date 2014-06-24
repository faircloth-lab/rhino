#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
(c) 2014 Brant Faircloth || http://faircloth-lab.org/
All rights reserved.

This code is distributed under a 3-clause BSD license. Please see
LICENSE.txt for more information.

Created on 23 June 2014 17:37 PDT (-0700)
"""

from rhino import core
from rhino.log import setup_logging


def main(args):
    # print message
    print welcome_message()
    # setup logging
    log, my_name = setup_logging(args)
    log.info("Getting alignments")
    # get arguments
    args = get_args()
    # make output dir
    args.output = tapir.create_unique_dir(args.output)
    # correct branch lengths
    tree_depth, correction, tree = tapir.correct_branch_lengths(args.tree, args.tree_format, d = args.output)
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


    # ----------------------
    # end
    text = " Completed {} ".format(my_name)
    log.info(text.center(65, "="))
