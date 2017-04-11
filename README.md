# Martinize


MARTINIZE.py is a script to create Coarse Grain Martini input files of
proteins, ready for use in the molecular dynamics simulations package
Gromacs. For more information on the Martini forcefield, see:
www.cgmartini.nl
and read our papers:
Monticelli et al., J. Chem. Theory Comput., 2008, 4(5), 819-834
de Jong et al., J. Chem. Theory Comput., 2013, DOI:10.1021/ct300646g

Primary input/output
--------------------
The input file (-f) should be a coordinate file in PDB or GROMOS
format. The format is inferred from the structure of the file. The
input can also be provided through stdin, allowing piping of
structures. The input structure can have multiple frames/models. If an output
structure file (-x) is given, each frame will be coarse grained,
resulting in a multimodel output structure. Having multiple frames may
also affect the topology. If secondary structure is determined
internally, the structure will be averaged over the frames. Likewise,
interatomic distances, as used for backbone bond lengths in Elnedyn
and in elastic networks, are also averaged over the frames available.

If an output file (-o) is indicated for the topology, that file will
be used for the master topology, using #include statements to link the
moleculetype definitions, which are written to separate files. If no
output filename is given, the topology and the moleculetype
definitions are written to stdout.

Secondary structure
-------------------
The secondary structure plays a central role in the assignment of atom
types and bonded interactions in MARTINI. Martinize allows
specification of the secondary structure as a string (-ss), or as a
file containing a specification in GROMACS' ssdump format
(-ss). Alternatively, DSSP can be used for an on-the-fly assignment of
the secondary structure. For this, the option -dssp has to be used
giving the location of the executable as the argument.
The option -collagen will set the whole structure to collagen. If this
is not what you want (eg only part of the structure is collagen, you
can give a secondary structure file/string (-ss) and specifiy collagen
as "F". Parameters for collagen are taken from: Gautieri et al.,
J. Chem. Theory Comput., 2010, 6, 1210-1218.
With multimodel input files, the secondary structure as determined with
DSSP will be averaged over the frames. In this case, a cutoff
can be specified (-ssc) indicating the fraction of frames to match a
certain secondary structure type for designation.

Topology
--------
Several options are available to tune the resulting topology. By
default, termini are charged, and chain breaks are kept neutral. This
behaviour can be changed using -nt and -cb, respectively.

Disulphide bridges can be specified using -cys. This option can be
given multiple times on the command line. The argument is a pair of
cysteine residues, using the format
chain/resn/resi,chain/resn/resi.
It is also possible to let martinize detect cysteine pairs based on a
cut-off distance of 0.22nm, by giving the keyword 'auto' as argument to -cys.
Alternatively, a different cut-off distance can be specified, which
will also trigger a search of pairs satisfying the distance
criterion (eg: -cys 0.32).

In addition to cystine bridges, links between other atoms can be
specified using -link. This requires specification of the atoms, using
the format
chain/resi/resn/atom,chain/resi/resn/atom,bondlength,forceconstant.
If only two atoms are given, a constraint will be added with length
equal to the (average) distance in the coordinate file. If a bond
length is added, but no force constant, then the bondlength will be
used to set a constraint.

Linking atoms requires that the atoms are part of the same
moleculetype. Therefore any link between chains will cause the chains
to be merged. Merges can also be specified explicitly, using the
option -merge with a comma-separated list of chain identifiers to be
joined into one moleculetype. The option -merge can be used several
times. Note that specifying a chain in several merge groups will cause
all chains involved to be merged into a single moleculetype.

The moleculetype definitions are written to topology include (.itp)
files, using a name consisting of the molecule class (e.g. Protein)
and the chain identifier. With -name a name can be specified instead.
By default, martinize only writes a moleculetype for each unique
molecule, inferred from the sequence and the secondary structure
definition. It is possible to force writing a moleculetype definition
for every single molecule, using -sep.

The option -p can be used to write position restraints, using the
force constant specified with -pf, which is set to 1000 kJ/mol
by default.

For stability, elastic bonds are used to retain the structure of
extended strands. The option -ed causes dihedrals to be used
instead.

Different forcefields can be specified with -ff. All the parameters and
options belonging to that forcefield  will be set (eg. bonded interactions,
BB-bead positions, Elastic Network, etc.). By default martini 2.1 is
used.

Elastic network
---------------
Martinize can write an elastic network for atom pairs within a cutoff
distance. The force constant (-ef) and the upper distance bound (-eu)
can be speficied. If a force field with an intrinsic Elastic
network is specified (eg. Elnedyn) with -ff, -elastic in implied and
the default values for the force constant and upper cutoff are used.
However, these can be overwritten.

Multiscaling
------------
Martinize can process a structure to yield a multiscale system,
consisting of a coordinate file with atomistic parts and
corresponding, overlaid coarsegrained parts. For chains that are
multiscaled, rather than writing a full moleculetype definition,
additional [atoms] and [virtual_sitesn] sections are written, to
be appended to the atomistic moleculetype definitions.
The option -multi can be specified multiple times, and takes a chain
identifier as argument. Alternatively, the keyword 'all' can be given
as argument, causing all chains to be multiscaled.

========================================================================


         -f  Input file (PDB|GRO)
         -o  Output topology (TOP)
         -x  Output coarse grained structure (PDB)
         -n  Output index file with CG (and multiscale) beads.
      -nmap  Output index file containing per bead mapping.
         -v  Verbose. Be load and noisy.
         -h  Display this help.
        -ss  Secondary structure (File or string)
       -ssc  Cutoff fraction for ss in case of ambiguity (default: 0.5).
      -dssp  DSSP executable for determining structure
  -collagen  Use collagen parameters
       -his  Interactively set the charge of each His-residue.
        -nt  Set neutral termini (charged is default)
        -cb  Set charges at chain breaks (neutral is default)
       -cys  Disulphide bond (+)
      -link  Link (+)
     -merge  Merge chains: e.g. -merge A,B,C (+)
      -name  Moleculetype name
         -p  Output position restraints (None/All/Backbone) (default: None)
        -pf  Position restraints force constant (default: 1000 kJ/mol/nm^2)
        -ed  Use dihedrals for extended regions rather than elastic bonds)
       -sep  Write separate topologies for identical chains.
        -ff  Which forcefield to use: 
   -elastic  Write elastic bonds
        -ef  Elastic bond force constant Fc
        -el  Elastic bond lower cutoff: F = Fc if rij < lo
        -eu  Elastic bond upper cutoff: F = 0  if rij > up
        -ea  Elastic bond decay factor a
        -ep  Elastic bond decay power p
        -em  Remove elastic bonds with force constant lower than this
        -eb  Comma separated list of bead names for elastic bonds
     -multi  Chain to be set up for multiscaling (+)


