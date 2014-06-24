#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
(c) 2014 Brant Faircloth || http://faircloth-lab.org/
All rights reserved.

This code is distributed under a 3-clause BSD license. Please see
LICENSE.txt for more information.

Created on 23 June 2014 17:38 PDT (-0700)
"""

from __future__ import absolute_import
import os
import subprocess

#import pdb

# get a dynamic version number, if possible.  if not running from git
# should default to static version
cwd = os.getcwd()
try:
    location = os.path.split(os.path.abspath(__file__))[0]
    os.chdir(location)
    cmd = [
        "git",
        "rev-parse",
        "--short",
        "HEAD"
    ]
    proc = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT
    )
    stdout, stderr = proc.communicate()
    if stdout and stdout.startswith("fatal:"):
        raise IOError("{}".format(stdout.strip()))
    else:
        __version__ = "git {}".format(stdout.strip())
    os.chdir(cwd)

except:
    __version__ = "1.0.0"
    if not os.getcwd == cwd:
        os.chdir(cwd)
