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

##########################
## 4 # FG -> CG MAPPING ##  -> @MAP <-
##########################
import functions


dnares3 = " DA DC DG DT"
dnares1 = " dA dC dG dT"
rnares3 = "  A  C  G  U"
rnares1 = " rA rC rG rU"

# Amino acid nucleic acid codes:
# The naming (AA and '3') is not strictly correct when adding DNA/RNA, but we keep it like this for consistincy.
AA3     = functions.spl("TRP TYR PHE HIS HIH ARG LYS CYS ASP GLU ILE LEU MET ASN PRO HYP GLN SER THR VAL ALA GLY"+dnares3+rnares3) #@#
AA1     = functions.spl("  W   Y   F   H   H   R   K   C   D   E   I   L   M   N   P   O   Q   S   T   V   A   G"+dnares1+rnares1) #@#

# Dictionaries for conversion from one letter code to three letter code v.v.
AA123, AA321 = functions.hash(AA1, AA3), functions.hash(AA3, AA1)

# Residue classes:
protein = AA3[:-8]   # remove eight to get rid of DNA/RNA here.
water   = functions.spl("HOH SOL TIP")
lipids  = functions.spl("DPP DHP DLP DMP DSP POP DOP DAP DUP DPP DHP DLP DMP DSP PPC DSM DSD DSS")
nucleic = functions.spl("DAD DCY DGU DTH ADE CYT GUA THY URA DA DC DG DT")

residueTypes = dict(
    [(i, "Protein") for i in protein ] +
    [(i, "Water")   for i in water   ] +
    [(i, "Lipid")   for i in lipids  ] +
    [(i, "Nucleic") for i in nucleic ]
    )


# Determine average position for a set of weights and coordinates
# This is a rather specific function that requires a list of items
# [(m,(x,y,z),id),..] and returns the weighted average of the
# coordinates and the list of ids mapped to this bead
def aver(b):
    mwx, ids = zip(*[((m*x, m*y, m*z), i) for m, (x, y, z), i in b])      # Weighted coordinates
    tm  = sum(zip(*b)[0])                                                 # Sum of weights
    return [sum(i)/tm for i in zip(*mwx)], ids                            # Centre of mass


# Return the CG beads for an atomistic residue, using the mapping specified above
# The residue 'r' is simply a list of atoms, and each atom is a list:
# [ name, resname, resid, chain, x, y, z ]
def map(r, ff, ca2bb = False):
    p = ff.mapping[r[0][1]]     # Mapping for this residue
    if ca2bb: p[0] = ["CA"]                # Elnedyn maps BB to CA, ca2bb is False or True
    # Get the name, mass and coordinates for all atoms in the residue
    a = [(i[0], ff.mass.get(i[0][0], 0), i[4:]) for i in r]
    # Store weight, coordinate and index for atoms that match a bead
    q = [[(m, coord, a.index((atom, m, coord))) for atom, m, coord in a if atom in i] for i in p]

    # Bead positions
    return zip(*[aver(i) for i in q])


# Mapping for index file
def mapIndex(r, ff, ca2bb = False):
    p = ff.mapping[r[0][1]]        # Mapping for this residue
    if ca2bb: p[0] = ["CA"]                   # Elnedyn maps BB to CA, ca2bb is False or True
    # Get the name, mass and coordinates for all atoms in the residue
    a = [(i[0], ff.mass.get(i[0][0], 0), i[4:]) for i in r]
    # Store weight, coordinate and index for atoms that match a bead
    return [[(m, coord, a.index((atom, m, coord))) for atom, m, coord in a if atom in i] for i in p]
