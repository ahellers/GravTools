"""
gravtools
=========

Code by Andreas Hellerschmied
andeas.hellerschmid@bev.gv.at

Summary
-------
Command line interface o gravtools.

This module contains functions that are accessible from the command line via console_scripts entry points of setuptool.

References
----------
.. https://python-packaging.readthedocs.io/en/latest/command-line-scripts.html
"""


from bev_legacy.schwaus import main as schwaus_main
import argparse
import os
from bev_legacy.settings import OUT_PATH, PATH_OESGN_TABLE, NAME_OESGN_TABLE


def is_file(filename):
    """Check, whether the input string is the path to an existing file."""
    if os.path.isfile(filename):
        return filename
    raise argparse.ArgumentTypeError("'{}' is not a valid file.".format(filename))


def is_dir(pathname):
    """Check, whether the input string is a valid and existing filepath."""
    if os.path.exists(pathname):
        return pathname
    raise argparse.ArgumentTypeError("'{}' is not a valid directory.".format(pathname))


def schwaus():
    """Command line interface including argument parser for the legacy BEV gravity survey processing scheme."""
    parser = argparse.ArgumentParser(prog="schwaus",
                                     description="Legacy gravity survey processing scheme at BEV.",
                                     epilog="The observatios are drift-corrected by fitting a polynomial of "
                                            "degree 1 to 3 using multiple linear regression. In the second step, the "
                                            "corrected observations are fitted onto the ÖSGN level.",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--obs-file", type=is_file, required=True, help="path and name of observation file")
    parser.add_argument("--oesgn-table", required=False, type=is_file, help="path and name of ÖSGN table",
                        default=PATH_OESGN_TABLE + NAME_OESGN_TABLE)
    parser.add_argument("--out-dir", type=is_dir, required=False,
                        help="path to output directory", default=OUT_PATH)
    args = parser.parse_args()
    # print(args)

    # Prepare paths and file names for function call:
    path_oesgn_table_str, name_oesgn_table_str = os.path.split(args.oesgn_table)
    path_obs_file_str, name_obs_file_str = os.path.split(args.obs_file)
    if path_oesgn_table_str:  # String not empty => Append '/'
        path_oesgn_table_str += '/'
    if path_obs_file_str:  # String not empty => Append '/'
        path_obs_file_str += '/'
    args.out_dir += '/'

    return schwaus_main(path_oesgn_table_str, name_oesgn_table_str, path_obs_file_str, name_obs_file_str, args.out_dir)

