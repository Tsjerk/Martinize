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


def atom(a):
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


def atoms(s):
    """Return a list of atoms from the command line"""
    return [atom(i) for i in s.split(",")]


class Link:
    def __init__(self, s=None, a=None, b=None, length=None, fc=None):
        """Parse a link string: two atoms and optional parameters"""
        self.a = a
        self.b = b
        self.length = length
        self.fc = fc

        if s:
            self.from_string(s)

    def from_string(s):
        ln     = s.split(",")
        self.a = ln[0]
        self.b = ln[1]

        if len(ln) > 3:  # Bond with given length and force constant
            self.length = float(ln[2]) if ln[2] else None
            self.force = float(ln[3])
        elif len(a) == 3:  # Constraint at given distance
            self.length = float(a[2])


