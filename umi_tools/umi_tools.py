'''
umi_tools.py - Tools for UMI analyses
===============================================

:Author: Tom Smith & Ian Sudbury, CGAT

***Modified by Ben Demaree, UCSF***

:Release: $Id$
:Date: |today|
:Tags: Genomics UMI

There are 6 tools:

  - whitelist
  - extract
  - group
  - dedup
  - count
  - count_tab

To get help on a specific tool, type:

    umi_tools <tool> --help

To use a specific tool, type::

    umi_tools <tool> [tool options] [tool arguments]
'''

import os
import sys
import imp


def main(argv=None):

    argv = sys.argv

    path = os.path.abspath(os.path.dirname(__file__))

    # if len(argv) == 1 or argv[1] == "--help" or argv[1] == "-h":
    #     print(globals()["__doc__"])
    #
    #     return
    #
    # elif len(argv) == 1 or argv[1] == "--version" or argv[1] == "-v":
    #     # collect umi_tools version
    #     sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))
    #     import version
    #
    #     print("UMI-tools version: %s" % version.__version__)
    #
    #     return

    command = argv[1]

    (file, pathname, description) = imp.find_module(command, [path, ])
    module = imp.load_module(command, file, pathname, description)
    # remove 'umi-tools' from sys.argv
    del sys.argv[0]
    module.main(sys.argv)

if __name__ == "__main__":
    sys.exit(main())
