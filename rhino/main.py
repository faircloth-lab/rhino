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
    # setup logging
    log, my_name = setup_logging(args)
    log.info("Getting alignments")


    # ----------------------
    # end
    text = " Completed {} ".format(my_name)
    log.info(text.center(65, "="))
