#!/usr/bin/env python

###############################################################################
#
#   Setuptools hook to specify config and command-line functions
#
#   AUTHOR: Krish Agarwal
#   AFFILIATION: University_of_Basel
#   CONTACT: akrish136@gmail.com
#   CREATED: 18-10-2020
#   LICENSE: Apache_2.0
#
###############################################################################

from __future__ import print_function
import logging

from multiqc.utils import config

# Initialise the main MultiQC logger
log = logging.getLogger("multiqc")


# Add default config options for the things that are used in MultiQC_NGI
def ALFA_execution_start():
    """Code to execute after the config files and
    command line flags have been parsedself.
    This setuptools hook is the earliest that will be able
    to use custom command line flags.
    """

    # Add to the search patterns used by modules
    if "ALFA" not in config.sp:
        config.update_dict(config.sp, {"ALFA": {"fn": "*ALFA_feature_counts.tsv"}})
