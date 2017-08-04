# MARTINIZE
# A simple, versatile tool for coarse-graining molecular systems
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


import sys, logging, random, math, os, re

from . import IO, topology, elastic, functions, mapping
from .converters import Link


def write_structure(ostream, model, title, box, chains, order):
    logging.info("Writing coarse grained structure.")
    ostream.write("MODEL %8d\n" % model)
    # TODO The title lines should not be more than 80 columns, including
    # the TITLE prefix.
    if not title.startswith('TITLE'):
        title = 'TITLE     ' + title
    ostream.write(title)
    ostream.write(IO.pdbBoxString(box))
    atid = 1
    for i in order:
        ci = chains[i]
        if ci.multiscale:
            for r in ci.residues:
                for name, resn, resi, chain, x, y, z in r:
                    ostream.write(IO.pdbOut((name, resn[:3], resi, chain, x, y, z),i=atid))
                    atid += 1
        coarseGrained = ci.cg(com=True)
        if coarseGrained:
            for name, resn, resi, chain, x, y, z, ssid in coarseGrained:
                if ci.multiscale:
                    name = "v"+name
                ostream.write(IO.pdbOut((name, resn[:3], resi, chain, x, y, z),i=atid,ssid=ssid))
                atid += 1
            ostream.write("TER\n")
        else:
            logging.warning("No mapping for coarse graining chain %s (%s); chain is skipped." % (ci.id, ci.type()))
    ostream.write("ENDMDL\n")


def average_secstruc(sslist):
    # Collect the secondary structure stuff and decide what to do with it
    # First rearrange by the residue
    ssTotal = zip(*sslist)
    ssAver  = []
    for i in ssTotal:
        si = list(set(i))
        if len(si) == 1:
            # Only one type -- consensus
            ssAver.append(si[0])
        else:
            # Transitions between secondary structure types
            i = list(i)
            si = [(1.0*i.count(j)/len(i), j) for j in si]
            si.sort()
            if si[-1][0] > options["sscutoff"]:
                ssAver.append(si[-1][1])
            else:
                ssAver.append(" ")

    ssAver = "".join(ssAver)
    logging.info('(Average) Secondary structure has been determined (see head of .itp-file).')

    return ssAver


def determine_secondary_structure(chains, options):
    ## SECONDARY STRUCTURE
    ss = ''
    if options['collagen']:
        for chain in chains:
            chain.set_ss("F")
            ss += chain.ss
    elif options["secstruc"]:
        # XXX We need error-catching here,
        # in case the file doesn't exist, or the string contains bogus.
        # If the string given for the sequence consists strictly of upper case letters
        # and does not appear to be a file, assume it is the secondary structure
        ss = options["secstruc"].replace('~', 'L').replace(' ', 'L')
        if ss.isalnum() and ss.isupper() and not os.path.exists(options["secstruc"]):
            ss = options["secstruc"]
            logging.info('Secondary structure read from command-line:\n'+ss)
        else:
            # There ought to be a file with the name specified
            ssfile = [i.strip() for i in open(options["secstruc"])]

            # Try to read the file as a Gromacs Secondary Structure Dump
            # Those have an integer as first line
            if ssfile[0].isdigit():
                logging.info('Will read secondary structure from file (assuming Gromacs ssdump).')
                ss = "".join([i for i in ssfile[1:]])
            else:
                # Get the secondary structure type from DSSP output
                logging.info('Will read secondary structure from file (assuming DSSP output).')
                pss = re.compile(r"^([ 0-9]{4}[0-9]){2}")
                ss  = "".join([i[16] for i in open(options["secstruc"]) if re.match(pss, i)])

        # Now set the secondary structure for each of the chains
        sstmp = ss
        for chain in chains:
            ln = min(len(sstmp), len(chain))
            chain.set_ss(sstmp[:ln])
            sstmp = ss[:ln]
    else:
        if options["dsspexe"]:
            method, executable = "dssp", options["dsspexe"]
        #elif options["pymol"]:
        #    method, executable = "pymol", options["pymol"]
        else:
            logging.warning("No secondary structure or determination method speficied. Protein chains will be set to 'COIL'.")
            method, executable = None, None

        for chain in chains:
            ss += chain.dss(method, executable)

        # Used to be: if method in ("dssp","pymol"): but pymol is not supported
        if method in ["dssp"]:
            logging.debug('%s determined secondary structure:\n' % method.upper()+ss)

    # Collect the secondary structure classifications for different frames
    return ss


def read_input_file(options):
    # Check whether to read from a gro/pdb file or from stdin
    # We use an iterator to wrap around the stream to allow
    # inferring the file type, without consuming lines already
    inStream = IO.streamTag(options["input"] or sys.stdin)

    # The streamTag iterator first yields the file type, which
    # is used to specify the function for reading frames
    fileType = inStream.next()
    if fileType == "GRO":
        frameIterator = IO.groFrameIterator
    else:
        frameIterator = IO.pdbFrameIterator

    # ITERATE OVER FRAMES IN STRUCTURE FILE #

    # Now iterate over the frames in the stream
    # This should become a StructureFile class with a nice .next method
    model     = 1
    cgOutPDB  = options["outstruc"] and open(options["outstruc"], "w")
    ssTotal   = []
    cysteines = []
    for title, atoms, box in frameIterator(inStream):
        if fileType == "PDB":
            # The PDB file can have chains, in which case we list and process them specifically
            # TER statements are also interpreted as chain separators
            # A chain may have breaks in which case the breaking residues are flagged
            chains = [IO.Chain(options, [i for i in IO.residues(chain)]) for chain in IO.pdbChains(atoms)]
        else:
            # The GRO file does not define chains. Here breaks in the backbone are
            # interpreted as chain separators.
            residuelist = [residue for residue in IO.residues(atoms)]
            # The breaks are indices to residues
            broken = IO.breaks(residuelist)
            # Reorder, such that each chain is specified with (i,j,k)
            # where i and j are the start and end of the chain, and
            # k is a chain identifier
            chains = zip([0]+broken, broken+[len(residuelist)], range(len(broken)+1))
            chains = [IO.Chain(options, residuelist[i:j], name=chr(65+k)) for i, j, k in chains]

        for chain in chains:
            chain.multiscale = "all" in options['multi'] or chain.id in options['multi']

        # Check the chain identifiers
        if model == 1 and len(chains) != len(set([i.id for i in chains])):
            # Ending down here means that non-consecutive blocks of atoms in the
            # PDB file have the same chain ID. The warning pertains to PDB files only,
            # since chains from GRO files get a unique chain identifier assigned.
            logging.warning("Several chains have identical chain identifiers in the PDB file.")

        # Check if chains are of mixed type. If so, split them.
        # Note that in some cases HETATM residues are part of a
        # chain. This will get problematic. But we cannot cover
        # all, probably.
        #***
        if True or not options['MixedChains']:
            demixedChains = []
            for chain in chains:
                demixedChains.extend(chain.split())
            chains = demixedChains

        n = 1
        logging.info("Found %d chains:" % len(chains))
        for chain in chains:
            logging.info("  %2d:   %s (%s), %d atoms in %d residues." % (n, chain.id, chain._type, chain.natoms, len(chain)))
            n += 1

        # Check all chains
        keep = []
        for chain in chains:
            if chain.type() == "Water":
                logging.info("Removing %d water molecules (chain %s)." % (len(chain), chain.id))
            elif chain.type() in ("Protein", "Nucleic"):
                keep.append(chain)
            # This is currently not active:
            elif False and options['RetainHETATM']: #***
                keep.append(chain)
            else:
                logging.info("Removing HETATM chain %s consisting of %d residues." % (chain.id, len(chain)))
        chains = keep

        # Here we interactively check the charge state of resides
        # Can be easily expanded to residues other than HIS
        for chain in chains:
            for i, resname in enumerate(chain.sequence):
                if resname == 'HIS' and options['sethischarge']:
                    choices = {0: 'HIH', 1: 'HIS'}
                    choice = IO.getChargeType(resname, i, choices)
                    chain.sequence[i] = choice

        # Check which chains need merging
        if model == 1:
            order, merge = IO.check_merge(chains, options['merges'], options['links'], 
                                          options['CystineCheckBonds'] and options['CystineMaxDist2'])

        # Get the total length of the sequence
        seqlength = sum([len(chain) for chain in chains])
        logging.info('Total size of the system: %s residues.' % seqlength)

        ssTotal.append(determine_secondary_structure(chains, options))

        # Write the coarse grained structure if requested
        if cgOutPDB:
            write_structure(cgOutPDB, model, title, box, chains, order)

        # Gather cysteine sulphur coordinates
        cyslist = [cys["SG"] for chain in chains for cys in chain["CYS"]]
        cysteines.append([cys for cys in cyslist if cys])

        model += 1

    ## DONE READING STRUCTURE
    return chains, atoms, ssTotal, cysteines, merge, order
    

def write_index(indexfile, chains, order):
    # Write the index file if requested.
    # Mainly of interest for multiscaling.
    # Could be improved by adding separate groups for BB, SC, etc.

    logging.info("Writing index file.")
    # Lists for All-atom, Virtual sites and Coarse Grain.
    NAA, NVZ, NCG = [], [], []
    atid = 1
    print(order)
    for i in order:
        ci = chains[i]
        coarseGrained = ci.cg(force=True)
        if ci.multiscale:
            NAA.extend([" %5d" % (a+atid) for a in range(ci.natoms)])
            atid += ci.natoms
        if coarseGrained:
            if ci.multiscale:
                NVZ.extend([" %5d" % (a+atid) for a in range(len(coarseGrained))])
            else:
                NCG.extend([" %5d" % (a+atid) for a in range(len(coarseGrained))])
            atid += len(coarseGrained)
    outNDX = open(indexfile, "w")
    outNDX.write("\n[ AA ]\n"+"\n".join([" ".join(NAA[i:i+15]) for i in range(0, len(NAA), 15)]))
    outNDX.write("\n[ VZ ]\n"+"\n".join([" ".join(NVZ[i:i+15]) for i in range(0, len(NVZ), 15)]))
    outNDX.write("\n[ CG ]\n"+"\n".join([" ".join(NCG[i:i+15]) for i in range(0, len(NCG), 15)]))
    outNDX.close()

    return


def write_mapping_index(filename, atoms):
    logging.info("Writing trajectory index file.")
    atid = 1
    outNDX = open(filename, "w")
    # Get all AA atoms as lists of atoms in residues
    # First we skip hetatoms and unknowns then iterate over beads
    # In DNA the O3' atom is mapped together with atoms from the next residue
    # This stores it until we get to the next residue
    o3_shift = ''
    for i_count, i in enumerate(IO.residues(atoms)): ## 'atoms' contains last frame read
        if i[0][1] in ("SOL", "HOH", "TIP"):
            continue
        if not i[0][1] in mapping.CoarseGrained.mapping.keys():
            continue
        nra = 0
        names = [j[0] for j in i]
        # This gives out a list of atoms in residue, each tuple has other
        # stuff in it that's needed elsewhere so we just take the last
        # element which is the atom index (in that residue)
        for j_count, j in enumerate(mapping.mapIndex(i)):
            outNDX.write('[ Bead %i of residue %i ]\n' % (j_count+1, i_count+1))
            line = ''
            for k in j:
                if names[k[2]] == "O3'":
                    line += '%s ' % (str(o3_shift))
                    o3_shift = k[2]+atid
                else:
                    line += '%i ' % (k[2]+atid)
            line += '\n'
            nra += len(j)
            outNDX.write(line)
        atid += nra

    return


def write_topology(filename, name, molecules, moleculeTypes):
    # WRITING THE MASTER TOPOLOGY
    # Output stream
    top  = filename and open(filename, 'w') or sys.stdout

    # ITP file listing
    itps = '\n'.join(['#include "%s.itp"' % molecule for molecule in set(moleculeTypes.values())])

    # Molecule listing
    logging.info("Output contains %d molecules:" % len(molecules))
    n = 1
    for molecule in molecules:
        chainInfo = (n, moleculeTypes[molecule], len(molecule) > 1 and "s" or " ", " ".join([i.id for i in molecule]))
        logging.info("  %2d->  %s (chain%s %s)" % chainInfo)
        n += 1
    molecules   = '\n'.join(['%s \t 1' % moleculeTypes[molecule] for molecule in molecules])

    # Set a define if we are to use rubber bands
    # useRubber   = options['elastic'] and "#define RUBBER_BANDS" or ""

    # XXX Specify a better, version specific base-itp name.
    # Do not set a define for position restrains here, as people are more used to do it in mdp file?
    top.write(
'''#include "martini.itp"

%s

[ system ]
; name
Martini system from %s

[ molecules ]
; name        number
%s''' % (itps, name, molecules))

    logging.info('Written topology files')

    return


def cystine_bridges(cys, options):
    logging.info("Checking for cystine bridges, based on sulphur (SG) atoms lying closer than %.4f nm" % math.sqrt(options['CystineMaxDist2']/100))

    cyscoord  = zip(*[[j[4:7] for j in i] for i in cys])
    cysteines = [i[:4] for i in cys[0]]

    bl, kb    = options['ForceField'].special[(("SC1", "CYS"), ("SC1", "CYS"))]

    bridges = []

    # Check the distances and add the cysteines to the link list if the
    # SG atoms have a distance smaller than the cutoff.
    rlc = range(len(cysteines))
    for i in rlc[:-1]:
        for j in rlc[i+1:]:
            # Checking the minimum distance over all frames
            # But we could also take the maximum, or the mean
            d2 = min([functions.distance2(a, b) for a, b in zip(cyscoord[i], cyscoord[j])])
            if d2 <= options['CystineMaxDist2']:
                a, b = cysteines[i], cysteines[j]
                bridges.append(Link(a=("SC1", "CYS", a[2], a[3]), 
                                    b=("SC1", "CYS", b[2], b[3]), 
                                    length=bl, fc=kb))
                a, b = (a[0], a[1], a[2]-(32 << 20), a[3]), (b[0], b[1], b[2]-(32 << 20), b[3])
                logging.info("Detected SS bridge between %s and %s (%f nm)" % (a, b, math.sqrt(d2)/10))
    
    return bridges


def elastic_network(atoms, coords, options):
    rubberType = options['ForceField'].EBondType
    rubberList = elastic.rubberBands(
        [(i[0], j) for i, j in zip(atoms, coords) if i[4] in options['elbeads']],
        options['ellowerbound'], options['elupperbound'],
        options['eldecay'], options['elpower'],
        options['elastic_fc'], options['elminforce'])
    return [ topology.Bond(i, options=options, type=rubberType, category="Rubber band") for i in rubberList]



def links(linklist, atoms, forcefield, options):

    # Run through the link list and add connections (links = cys bridges or hand specified links)
    out = []
    for lnk in linklist:
        atomA, atomB, bondlength, forceconst = lnk.a, lnk.b, lnk.length, lnk.fc
        if bondlength == -1 and forceconst == -1:
            bondlength, forceconst = forcefield.special[(atomA[:2], atomB[:2])]
        # Check whether this link applies to this group
        atomA = atomA in atoms and atoms.index(atomA)+1
        atomB = atomB in atoms and atoms.index(atomB)+1
        if atomA and atomB:
            cat = (forceconst is None) and "Constraint" or "Link"
            out.append(topology.Bond(
                (atomA, atomB),
                options    = options,
                type       = 1,
                parameters = (bondlength, forceconst),
                category   = cat,
                comments   = "Cys-bonds/special link"))

    return out


def do_topology(options, chains, ssTotal, cysteines, merge):

    ssAver = average_secstruc(ssTotal)

    # Divide the secondary structure according to the division in chains
    # This will set the secondary structure types to be used for the
    # topology.
    for chain in chains:
        chain.set_ss(ssAver[:len(chain)])
        ssAver = ssAver[len(chain):]

    # Now the chains are complete, each consisting of a residuelist,
    # and a secondary structure designation if the chain is of type 'Protein'.
    # There may be mixed chains, there may be HETATM things.
    # Water has been discarded. Maybe this has to be changed at some point.
    # The order in the coarse grained files matches the order in the set of chains.
    #
    # If there are no merges to be done, i.e. no global Elnedyn network, no
    # disulphide bridges, no links, no distance restraints and no explicit merges,
    # then we can write out the topology, which will match the coarse grained file.
    #
    # If there are merges to be done, the order of things may be changed, in which
    # case the coarse grained structure will not match with the topology...

    # CYSTINE BRIDGES #
    # Extract the cysteine coordinates (for all frames) and the cysteine identifiers
    if options['CystineCheckBonds']:
        options["links"].extend(cystine_bridges(cysteines, options))

    # REAL ITP STUFF #
    # Check whether we have identical chains, in which case we
    # only write the ITP for one...
    # This means making a distinction between chains and
    # moleculetypes.

    molecules = [tuple([chains[i] for i in j]) for j in merge]

    # At this point we should have a list or dictionary of chains
    # Each chain should be given a unique name, based on the value
    # of options["-o"] combined with the chain identifier and possibly
    # a number if there are chains with identical identifiers.
    # For each chain we then write an ITP file using the name for
    # moleculetype and name + ".itp" for the topology include file.
    # In addition we write a master topology file, using the value of
    # options["-o"], with an added extension ".top" if not given.

    # XXX *NOTE*: This should probably be gathered in a 'Universe' class
    itp = 0
    moleculeTypes = {}
    for mi in range(len(molecules)):
        mol = molecules[mi]
        # Check if the moleculetype is already listed
        # If not, generate the topology from the chain definition
        if mol not in moleculeTypes or options['separate']:
            # Name of the moleculetype
            # XXX: The naming should be changed; now it becomes Protein_X+Protein_Y+...
            name = "+".join([chain.getname(options['name']) for chain in mol])
            moleculeTypes[mol] = name

            # Write the molecule type topology
            top = topology.Topology(mol[0], options=options, name=name)
            for m in mol[1:]:
                top += topology.Topology(m, options=options)

            # Have to add the connections, like the connecting network
            # Gather coordinates
            mcg, coords = zip(*[(j[:4], j[4:7]) for m in mol for j in m.cg(force=True)])
            mcg         = list(mcg)

            top.bonds.extend(links(options['links'], mcg, options["ForceField"], options))

            # Elastic Network
            # The elastic network is added after the topology is constructed, since that
            # is where the correct atom list with numbering and the full set of
            # coordinates for the merged chains are available.
            if options['elastic']:
                top.bonds.extend(elastic_network(top.atoms, coords, options))

            # Write out the MoleculeType topology
            destination = options["outtop"] and open(moleculeTypes[mol]+".itp", 'w') or sys.stdout
            destination.write(str(top))

            itp += 1

        # Check whether other chains are equal to this one
        # Skip this step if we are to write all chains to separate moleculetypes
        if not options['separate']:
            for j in range(mi+1, len(molecules)):
                if not molecules[j] in moleculeTypes and mol == molecules[j]:
                    # Molecule j is equal to a molecule mi
                    # Set the name of the moleculetype to the one of that molecule
                    moleculeTypes[molecules[j]] = moleculeTypes[mol]

    logging.info('Written %d ITP file%s' % (itp, itp > 1 and "s" or ""))

    name = options["input"] and options["input"] or "stdin"
    write_topology(options["outtop"], name, molecules, moleculeTypes)

    return


def main(options):
    chains, atoms, ssTotal, cysteines, merge, order = read_input_file(options)

    if options["index"]:
        write_index(options["index"], chains, order)

    # Write the index file for mapping AA trajectory if requested
    if options["mapping"]:
        write_mapping_index(options["mapping"], atoms)

    if options['outtop']:
        do_topology(options, chains, ssTotal, cysteines, merge)

    # Maybe there are forcefield specific log messages?
    options['ForceField'].messages()

    # The following lines are always printed (if no errors occur).
    print "\n\tThere you are. One MARTINI. Shaken, not stirred.\n"

    Q = random.choice(open(os.path.join(os.path.dirname(__file__), "quotes.txt")).readlines()).split('::')
    print "\n", Q[1].strip(), "\n%80s" % ("--"+Q[0].strip()), "\n"

    return 0
