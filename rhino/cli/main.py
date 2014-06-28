#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
(c) 2014 Brant Faircloth || http://faircloth-lab.org/
All rights reserved.

This code is distributed under a 3-clause BSD license. Please see
LICENSE.txt for more information.

Created on 23 June 2014 17:35 PDT (-0700)
"""

from __future__ import absolute_import
#import os
import argparse
from rhino import io


def get_args():
    parser = argparse.ArgumentParser(
        description="Partition by sites",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        '--alignments',
        required=True,
        help='A directory of PHYLIP alignments to use for tree inference.',
        action=io.FullPaths,
        type=io.is_dir
    )
    parser.add_argument(
        "--input-format",
        dest="input_format",
        choices=['fasta', 'nexus', 'phylip', 'phylip-relaxed', 'clustal', 'emboss', 'stockholm'],
        default='phylip',
        help="""The input alignment format"""
    )
    parser.add_argument(
        '--output',
        required=True,
        help='The output directory for results.',
        action=io.CreateDir
    )
    group = parser.add_mutually_exclusive_group()
    group.add_argument(
        '--tree-method',
        choices=['parsimony', 'likelihood'],
        default='parsimony',
        help='The rapid tree inference method to use.',
    )
    group.add_argument(
        '--tree-file',
        action=io.FullPaths,
        type=io.is_file,
        default=None,
        help='A phylogenetic tree to use.',
    )
    parser.add_argument(
        '--cores',
        type=int,
        default=1,
        help='The number of compute cores to use.',
    )
    parser.add_argument(
        "--verbosity",
        type=str,
        choices=["INFO", "WARN", "CRITICAL"],
        default="INFO",
        help="""The logging level to use."""
    )
    parser.add_argument(
        "--log-path",
        action=io.FullPaths,
        type=io.is_dir,
        default=None,
        help="""The path to a directory to hold logs."""
    )
    return parser.parse_args()


def main():
    from rhino.main import main
    args = get_args()
    main(args)
