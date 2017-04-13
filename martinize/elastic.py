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

#########################
## 7 # ELASTIC NETWORK ##  -> @ELN <-
#########################
import math
import functions

## ELASTIC NETWORK ##

# Only the decay function is defined here, the network
# itself is set up through the Topology class


# The function to determine the decay scaling factor for the elastic network
# force constant, based on the distance and the parameters provided.
# This function is very versatile and can be fitted to most commonly used
# profiles, including a straight line (rate=0)
def decayFunction(distance, shift, rate, power):
    return math.exp(-rate*math.pow(distance-shift, power))


def rubberBands(atomList, lowerBound, upperBound, decayFactor, decayPower, forceConstant, minimumForce):
    out = []
    u2  = upperBound**2
    while len(atomList) > 3:
        bi, xi = atomList.pop(0)
        for bj, xj in atomList[2:]:
            # Mind the nm/A conversion -- This has to be standardized! Global use of nm?
            d2 = functions.distance2(xi, xj)/100

            if d2 < u2:
                dij  = math.sqrt(d2)
                fscl = decayFunction(dij, lowerBound, decayFactor, decayPower)*forceConstant
                if fscl > minimumForce:
                    out.append({"atoms": (bi, bj), "parameters": (dij, fscl)})
    return out
