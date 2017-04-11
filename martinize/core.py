#!/usr/bin/env python


import sys, logging
from . import DOC
from . import cli
from . import martinize


# EDITABLE SECTIONS ARE MARKED WITH #@#


version = "2.6"
authors = ["Djurre H. de Jong", "Jaakko J. Uusitalo", "Tsjerk A. Wassenaar"]

# This program has grown to be pretty complex.
# The routines have been organized in different files.
# For working versions, all files can be incorporated by using the catenate.py file. 
#
# Index of the program files:
#
#   1. Options and documentation                             @DOC.py
#   2. Description, options and command line parsing         @CMD.py
#   3. Helper functions and macros                           @FUNC.py
#   4. Finegrained to coarsegrained mapping                  @MAP.py
#   5. Secondary structure determination and interpretation  @SS.py
#   6. Elastic network                                       @ELN.py
#   7. Structure I/O                                         @IO.py
#   8. Topology generation                                   @TOP.py
#   9. Main                                                  @MAIN.py
#
#   Force field parameters are specified in the differentt forcefield modules,
#   e.g.: martini22_ff.py


def main(argv):

    args = sys.argv[1:]
    # Get the possible commandline arguments arguments and help text. 
    options, lists = DOC.options, DOC.lists
    # Parse commandline options.
    options = cli.option_parser(args, options, lists, version)

    martinize.main(options)

    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv))
