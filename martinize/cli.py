# INSert membrANE
# A simple, versatile tool for building coarse-grained simulation systems
# Copyright (C) 2017  Tsjerk A. Wassenaar and contributors
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

from __future__ import absolute_import

import sys
import os
import logging
import inspect

import simopt
from simopt import MULTI, MA

from . import core
from . import DOC

#from .converters import vector, box3d, molspec


# Option list
OPTIONS = simopt.Options([
    #  level opt  attribute      type        num     default    flags    description
        """
    Input/output related options
    """,
    (0, "-v",        "verbose",        str,   1,          None,     0, "Verbosity level"),
    (0, "-f",        "input",          str,   1,          None,     0, "Input GRO or PDB file"),
    (0, "-o",        "outtop",         str,   1, "martini.top",     0, "Output topology (TOP)"),
    (0, "-x",        "outstruc",       str,   1,          None,     0, "Output coarse grained structure (PDB)"),
    (0, "-n",        "index",          str,   1,          None,     0, "Output index file with CG (and multiscale) beads."),
    (1, "-nmap",     "mapping",        str,   1,          None,     0, "Output index file containing per bead mapping."),
    (0, "-v",        "verbose",        bool,  0,         False,     0, "Verbose. Be load and noisy."),
    (1, "-ss",       "secstruc",       str,   1,          None,     0, "Secondary structure (File or string)"),
    (1, "-ssc",      "sscutoff",       float, 1,           0.5,     0, "Cutoff fraction for ss in case of ambiguity (default: 0.5)."),
    (0, "-dssp",     "dsspexe",        str,   1,          None,     0, "DSSP executable for determining structure"),
#    ("-pymol",  "pymolexe",  str,                      1,     None, "PyMOL executable for determining structure"),
    (0, "-collagen", "collagen",       bool,  0,         False,     0, "Use collagen parameters"),
    (1, "-his",      "sethischarge",   bool,  0,         False,     0, "Interactively set the charge of each His-residue."),
    (0, "-nt",       "neutraltermini", bool,  0,         False,     0, "Set neutral termini (charged is default)"),
    (1, "-cb",       "chargedbreaks",  bool,  0,         False,     0, "Set charges at chain breaks (neutral is default)"),
    (0, "-cys",      "cystines",       str,   1,          None, MULTI, "Disulphide bond (+)"),
    (1, "-link",     "links",          str,   1,          None, MULTI, "Link (+)"),
    (1, "-merge",    "merges",         str,   1,          None, MULTI, "Merge chains: e.g. -merge A,B,C (+)"),
    (0, "-name",     "name",           str,   1,          None,     0, "Moleculetype name"),
    (1, "-p",        "posres",         str,   1,        'None',     0, "Output position restraints (None/All/Backbone) (default: None)"),
    (1, "-pf",       "posrefc",        float, 1,          1000,     0, "Position restraints force constant (default: 1000 kJ/mol/nm^2)"),
    (1, "-ed",       "extdih",         bool,  0,         False,     0, "Use dihedrals for extended regions rather than elastic bonds)"),
    (1, "-sep",      "separate",       bool,  0,         False,     0, "Write separate topologies for identical chains."),
    (0, "-ff",       "forcefield",     str,   1,   'martini22',     0, "Which forcefield to use"),
# Fij = Fc exp( -a (rij - lo)**p )
    (1, "-elastic",  "elastic",        bool,  0,         False,     0, "Write elastic bonds"),
    (1, "-ef",       "elastic_fc",     float, 1,           500,     0, "Elastic bond force constant Fc"),
    (1, "-el",       "ellowerbound",   float, 1,             0,     0, "Elastic bond lower cutoff: F = Fc if rij < lo"),
    (1, "-eu",       "elupperbound",   float, 1,          0.90,     0, "Elastic bond upper cutoff: F = 0  if rij > up"),
    (1, "-ea",       "eldecay",        float, 1,             0,     0, "Elastic bond decay factor a"),
    (1, "-ep",       "elpower",        float, 1,             1,     0, "Elastic bond decay power p"),
    (1, "-em",       "elminforce",     float, 1,             0,     0, "Remove elastic bonds with force constant lower than this"),
    (1, "-eb",       "elbeads",        str,   1,          'BB',     0, "Comma separated list of bead names for elastic bonds"),
#    ("-hetatm", "hetatm",  bool,                     0,    False, "Include HETATM records from PDB file (Use with care!)"),
    (1, "-multi",    "multi",          str,   1,          None, MULTI, "Chain to be set up for multiscaling (+)"),
])




def str2atom(a):
    """ Helper function to parse atom strings given on the command line:
    resid, resname/resid, chain/resname/resid, resname/resid/atom,
    chain/resname/resid/atom, chain//resid, chain/resname/atom """
    a = a.split("/")
    if len(a) == 1:  # Only a residue number:
        return (None, None, int(a[0]), None)
    if len(a) == 2:  # Residue name and number (CYS/123):
        return (None, a[0], int(a[1]), None)
    if len(a) == 3:
        if a[2].isdigit():  # Chain, residue name, residue number
            return (None, a[1], int(a[2]), a[0])
        else:  # Residue name, residue number, atom name
            return (a[2], a[0], int(a[1]), None)
    return (a[3], a[1], int(a[2]), a[0])


def update_options(options):

    options["Version"] = ""

    # To make the program flexible, the forcefield parameters are defined
    # for multiple forcefield. We first check for a FF file in the current directory.
    # Next we check for the FF in globals (for catenated scripts). 
    # Next we check in at the location of the script and the subdiretory FF.
    try:
        options['ForceField'] = globals()[options['forcefield'].lower()]()
    except KeyError:
        try:
            _tmp = __import__(options['forcefield'].lower()+"_ff")
            options['ForceField'] = getattr(_tmp, options['forcefield'].lower())()
        except ImportError:
            try:
                # We add the directory where the script resides and a possible "ForceFields" directory to the search path
                # realpath() will make your script run, even if you symlink it :)
                cmd_folder = os.path.realpath(os.path.dirname(inspect.getfile(inspect.currentframe())))
                if cmd_folder not in sys.path:
                    sys.path.insert(0, cmd_folder)
                # use this if you want to include modules from a subfolder
                cmd_subfolder = os.path.realpath(os.path.dirname(inspect.getfile(inspect.currentframe()))) + "/ForceFields"
                if cmd_subfolder not in sys.path:
                     sys.path.insert(0, cmd_subfolder)
                _tmp = __import__(options['forcefield'].lower()+"_ff")
                options['ForceField'] = getattr(_tmp, options['forcefield'].lower())()
            except:
                logging.error("Forcefield '%s' can not be loaded." % (options['forcefield']))
                sys.exit()
        
    #    if os.path.exists(options['-ff'].value.lower()+'_ff.py'):
    #        _tmp = __import__(options['-ff'].value.lower()+"_ff")
    #        options['ForceField'] = getattr(_tmp, options['-ff'].value.lower())()
    #    elif os.path.exists('ForceFields/'+options['-ff'].value.lower()+'_ff.py'):
    #        _tmp = __import__("ForceFields."+options['-ff'].value.lower()+'_ff',fromlist="ForceFields")
    #        options['ForceField'] = getattr(_tmp, options['-ff'].value.lower())()
    #    elif options['-ff'].value.lower() in globals():
    #        options['ForceField'] = globals()[options['-ff'].value.lower()]()
    #    else:
    #        logging.error("Forcefield '%s' can not be found."%(options['-ff']))
    #        sys.exit()
    #except:
    #    logging.error("Forcefield '%s' can not be loaded."%(options['-ff']))
    #    sys.exit()

    # Process the raw options from the command line
    # Boolean options are set to more intuitive variables
    options['RetainHETATM']        = False  # options['-hetatm']
    options['MixedChains']         = False  # options['-mixed']

    options['elbeads'] = options['elbeads'].split(',')

    options['posres'] = [i.lower() for i in options['posres'].split(",")]
    if "backbone" in options['posres']: 
        options['posres'].append("BB")
    if "none" in options['posres']: 
        options['posres'] = []

    if options['ForceField'].ElasticNetwork:
        # Some forcefields, like elnedyn, always use an elatic network.
        # This is set in the forcefield file, with the parameter ElasticNetwork.
        options['elastic'] = True

    # Merges, links and cystines
    options['merges'] = "all" in options['merges'] and ["all"] or [i.split(",") for i in options['merges']]

    # Process links
    linkList   = []
    linkListCG = []
    for i in options['links']:
        ln     = i.split(",")
        a, b   = str2atom(ln[0]), str2atom(ln[1])
        if len(ln) > 3:  # Bond with given length and force constant
            bl, fc = (ln[2] and float(ln[2]) or None, float(ln[3]))
        elif len(a) == 3:  # Constraint at given distance
            bl, fc = float(a[2]), None
        else:  # Constraint at distance in structure
            bl, fc = None, None
        # Store the link, but do not list the atom name in the
        # atomistic link list. Otherwise it will not get noticed
        # as a valid link when checking for merging chains
        linkList.append(((None, a[1], a[2], a[3]), (None, b[1], b[2], b[3])))
        linkListCG.append((a, b, bl, fc))

    # Cystines
    # This should be done for all special bonds listed in the _special_ dictionary
    CystineCheckBonds = False   # By default, do not detect cystine bridges
    CystineMaxDist2   = (10*0.22)**2  # Maximum distance (A) for detection of SS bonds
    for i in options['cystines']:
        if i.lower() == "auto":
            CystineCheckBonds = True
        elif i.replace(".", "").isdigit():
            CystineCheckBonds = True
            CystineMaxDist2   = (10*float(i))**2
        else:
            # This item should be a pair of cysteines
            cysA, cysB = [str2atom(j) for j in i.split(",")]
            # Internally we handle the residue number shifted by ord(' ')<<20.
            # We have to add this to the cys-residue numbers given here as well.
            constant = 32 << 20
            linkList.append((("SG", "CYS", cysA[2]+constant, cysA[3]),
                            ("SG", "CYS", cysB[2]+constant, cysB[3])))
            linkListCG.append((("SC1", "CYS", cysA[2]+constant, cysA[3]),
                              ("SC1", "CYS", cysB[2]+constant, cysB[3]), -1, -1))

    # Now we have done everything to it, we can add Link/cystine related stuff to options
    # 'multi' is not stored anywhere else, so that we also add
    options['links'] = linkList
    options['cglinks'] = linkListCG
    options['CystineCheckBonds'] = CystineCheckBonds
    options['CystineMaxDist2']   = CystineMaxDist2

    ## LOGGING ##
    # Set the log level and communicate which options are set and what is happening
    # If 'Verbose' is set, change the logger level
    logLevel = options["verbose"] and logging.DEBUG or logging.INFO
    logging.basicConfig(format='%(levelname)-7s    %(message)s', level=logLevel)

    #logging.info('MARTINIZE, script version %s'%__version__)
    logging.info('If you use this script please cite:')
    logging.info('de Jong et al., J. Chem. Theory Comput., 2013, DOI:10.1021/ct300646g')

    logging.info("Chain termini will%s be charged"%(options['neutraltermini'] and " not" or ""))

    logging.info("Residues at chain brakes will%s be charged"%((not options['chargedbreaks']) and " not" or ""))

    if 'ForceField' in options:
        logging.info("The %s forcefield will be used."%(options['ForceField'].name))
    else:
        logging.error("Forcefield '%s' has not been implemented."%(options['forcefield']))
        sys.exit()

    if options['extdih']:
        logging.info('Dihedrals will be used for extended regions. (Elastic bonds may be more stable)')
    else:
        logging.info('Local elastic bonds will be used for extended regions.')

    if options['posres']:
        logging.info("Position restraints will be generated.")
        logging.warning("Position restraints are only enabled if -DPOSRES is set in the MDP file")

    if options['MixedChains']:
        logging.warning("So far no parameters for mixed chains are available. This might crash the program!")

    if options['RetainHETATM']:
        logging.warning("I don't know how to handle HETATMs. This will probably crash the program.")

    return options


def main(argv):
    ## TEMPORARY ---
    # Exception to be defined in martinize
    class MartinizeException(BaseException): pass
    ## <---

    ## OPTIONS
    # Parse options
    try:
        options = OPTIONS.parse(argv[1:])
        options["Arguments"] = argv[1:]
        update_options(options)
    except simopt.SimoptHelp:
        print(OPTIONS.help(argv[1:]))
        return 0
    except simopt.MissingMandatoryError as e:
        print(e)
        return 3
    except simopt.Usage as e:
        print(e)
        return 1

    ## WORK
    try:
        system = core.main(options)
    except MartinizeException as e:
        print(e)
        return 2

    ## OUTPUT
    # Build atom list
    # Build topology
    # Build index

    return 0


def cli():
    sys.exit(main(sys.argv))

