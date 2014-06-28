#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
(c) 2014 Brant Faircloth || http://faircloth-lab.org/
All rights reserved.

This code is distributed under a 3-clause BSD license. Please see
LICENSE.txt for more information.

Created on 23 June 2014 17:53 PDT (-0700)
"""


import os
import sys
import glob
import shutil
import argparse
from pkg_resources import resource_filename

from Bio import AlignIO
from Bio.Alphabet import IUPAC, Gapped


import pdb


class FullPaths(argparse.Action):
    """Expand user- and relative-paths"""
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, os.path.abspath(os.path.expanduser(values)))


class CreateDir(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        # get the full path
        d = os.path.abspath(os.path.expanduser(values))
        # check to see if directory exists
        if os.path.exists(d):
            answer = raw_input("[WARNING] Output directory exists, REMOVE [Y/n]? ")
            if answer == "Y":
                shutil.rmtree(d)
            else:
                print "[QUIT]"
                sys.exit()
        # create the new directory
        os.makedirs(d)
        # return the full path
        setattr(namespace, self.dest, d)


def is_dir(dirname):
    if not os.path.isdir(dirname):
        msg = "{0} is not a directory".format(dirname)
        raise argparse.ArgumentTypeError(msg)
    else:
        return dirname


def is_file(filename):
    if not os.path.isfile:
        msg = "{0} is not a file".format(filename)
        raise argparse.ArgumentTypeError(msg)
    else:
        return filename


def mkdir(path):
    try:
        os.mkdir(path)
    except OSError:
        pass


def create_unique_dir(path, limit=100):
    """Attempts to create a directory `path`. Returns the name of the
    directory actually created, which may or may not be the same as `path`.

    e.g., if my_directory already exists, it tries to create my_directory.1,
    my_directory.2, ... until my_directory.`limit` is reached.

    Race conditions are possible.
    """
    original = path
    count = 1
    while count < limit:
        try:
            os.mkdir(path)
            return path
        except OSError as e:
            if e.errno == 17: # file exists
                path = "{0}.{1}".format(original, count)
                count += 1
            else:
                raise
    else:
        msg = "could not uniquely create directory {0}: limit `{1}` reached"
        raise Exception(msg.format(original, limit))


def get_hyphy_conf():
    return resource_filename(__name__, 'data/models_and_rates.bf')


def get_list_from_ints(string, name = 'time'):
    """Convert times input as string to a list"""
    try:
        times = [int(i) for i in string.split(',')]
    except StandardError as e:
        msg = "Cannot convert {0} to a list of integers: {1}"
        raise argparse.ArgumentTypeError(msg.format(name, e))
    return times


def get_strings_from_items(string, name = 'locus'):
    """Convert times input as string to a list"""
    try:
        times = [str(i) for i in string.split(',')]
    except StandardError as e:
        msg = "Cannot convert {0} to a list of loci: {1}"
        raise argparse.ArgumentTypeError(msg.format(name, e))
    return times


def get_list_from_ranges(string):
    """Convert ranges entered as string to nested list"""
    try:
        ranges = [[int(j) for j in i.split('-')] for i in string.split(',')]
    except StandardError as e:
        msg = "Cannot convert spans to a list of integers: {0}"
        raise argparse.ArgumentTypeError(msg.format(e))
    return ranges


def get_alignment_files(input_dir, input_format):
    extensions = get_file_extensions(input_format)
    files = []
    for ext in extensions:
        files.extend(glob.glob(os.path.join(os.path.expanduser(input_dir), '*{}*'.format(ext))))
    # ensure we collapse duplicate filenames
    return list(set(files))


def get_file_extensions(ftype):
    ext = {
        'fasta': ('.fasta', '.fsa', '.aln', '.fa'),
        'nexus': ('.nexus', '.nex'),
        'phylip': ('.phylip', '.phy'),
        'phylip-relaxed': ('.phylip', '.phy'),
        'clustal': ('.clustal', '.clw'),
        'emboss': ('.emboss',),
        'stockholm': ('.stockholm',)
    }
    return ext[ftype]


def convert_alignment(aln, format, outdir, out_format):
    outfile_name = "{}.{}".format(
        os.path.splitext(os.path.basename(aln))[0],
        out_format
    )
    outfile = os.path.join(outdir, outfile_name)
    with open(outfile, "w") as outf:
        with open(aln, "rU") as infile:
            alignment = AlignIO.parse(infile, format, alphabet=Gapped(IUPAC.ExtendedIUPACDNA()))
            AlignIO.write(alignment, outf, out_format)
    return outfile
