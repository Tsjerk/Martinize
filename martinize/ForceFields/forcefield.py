
from .. import functions, secstruc
from .. import mapping

FORCE_FIELD_COLLECTION = {}


class _MetaForceField(type):
    def __init__(cls, name, bases, classdict):
        super(_MetaForceField, cls).__init__(name, bases, classdict)
        try:
            ff_name = classdict['name']
        except KeyError:
            raise KeyError('Class "{}" does not define a name.'.format(name))
        if ff_name is not None:
            ff_name = ff_name.lower()
            if ff_name in FORCE_FIELD_COLLECTION:
                raise RuntimeError('Force field "{}" already exists'.format(ff_name))
            FORCE_FIELD_COLLECTION[ff_name] = cls


class Forcefield(object):
    __metaclass__ = _MetaForceField
    name = None

    def __init__(self):
        self.set_mapping()
        self.setup()
        self.finish()

    def setup(self):
        raise NotImplementedError()

    def set_mapping(self):
        # Standard mapping groups
        # Protein backbone
        bb        = "N CA C O H H1 H2 H3 O1 O2"                                                                    #@#
        # Lipid tails
        palmitoyl1    = functions.nsplit("C1B C1C C1D C1E", "C1F C1G C1H C1I", "C1J C1K C1L C1M", "C1N C1O C1P")              #@#
        palmitoyl2    = functions.nsplit("C2B C2C C2D C2E", "C2F C2G C2H C2I", "C2J C2K C2L C2M", "C2N C2O C2P")              #@#
        oleyl1        = functions.nsplit("C1B C1C C1D C1E", "C1F C1G C1H", "C1I C1J", "C1K C1L C1M C1N", "C1O C1P C1Q C1R")   #@#
        oleyl2        = functions.nsplit("C2B C2C C2D C2E", "C2F C2G C2H", "C2I C2J", "C2K C2L C2M C2N", "C2O C2P C2Q C2R")   #@#
        #lauroyl1      = []
        #stearoyl1     = []
        #arachidonoyl1 = []
        #linoleyl1     = []
        #hexanoyl1     = []
        # Lipid head groups
        #phoshpatidylcholine      =
        phosphatydilethanolamine = functions.nsplit("N H1 H2 H3 CA", "CB P OA OB OC OD", "CC CD OG C2A OH", "CE OE C1A OF")      #@#
        phosphatidylglycerol     = functions.nsplit("H1 O1 CA H2 O2 CB", "CC P OA OB OC OD", "CD CE OG C2A OH", "CF OE C1A OF")  #@#
        #phosphatidylserine       =

        self.dna_bb = "P OP1 OP2 O5' O3'", "C5' O4' C4'", "C3' O3' C2' C1'"

        # This is the mapping dictionary
        # For each residue it returns a list, each element of which
        # lists the atom names to be mapped to the corresponding bead.
        # The order should be the standard order of the coarse grained
        # beads for the residue. Only atom names matching with those
        # present in the list of atoms for the residue will be used
        # to determine the bead position. This adds flexibility to the
        # approach, as a single definition can be used for different
        # states of a residue (e.g., GLU/GLUH).
        # For convenience, the list can be specified as a set of strings,
        # converted into a list of lists by 'functions.nsplit' defined above.
        self.mapping = {
            "ALA":  functions.nsplit(bb + " CB"),
            "CYS":  functions.nsplit(bb, "CB SG"),
            "ASP":  functions.nsplit(bb, "CB CG OD1 OD2"),
            "GLU":  functions.nsplit(bb, "CB CG CD OE1 OE2"),
            "PHE":  functions.nsplit(bb, "CB CG CD1 HD1", "CD2 HD2 CE2 HE2", "CE1 HE1 CZ HZ"),
            "GLY":  functions.nsplit(bb),
            "HIS":  functions.nsplit(bb, "CB CG", "CD2 HD2 NE2 HE2", "ND1 HD1 CE1 HE1"),
            "HIH":  functions.nsplit(bb, "CB CG", "CD2 HD2 NE2 HE2", "ND1 HD1 CE1 HE1"),     # Charged Histidine.
            "ILE":  functions.nsplit(bb, "CB CG1 CG2 CD CD1"),
            "LYS":  functions.nsplit(bb, "CB CG CD", "CE NZ HZ1 HZ2 HZ3"),
            "LEU":  functions.nsplit(bb, "CB CG CD1 CD2"),
            "MET":  functions.nsplit(bb, "CB CG SD CE"),
            "ASN":  functions.nsplit(bb, "CB CG ND1 ND2 OD1 OD2 HD11 HD12 HD21 HD22"),
            "PRO":  functions.nsplit(bb, "CB CG CD"),
            "HYP":  functions.nsplit(bb, "CB CG CD OD"),
            "GLN":  functions.nsplit(bb, "CB CG CD OE1 OE2 NE1 NE2 HE11 HE12 HE21 HE22"),
            "ARG":  functions.nsplit(bb, "CB CG CD", "NE HE CZ NH1 NH2 HH11 HH12 HH21 HH22"),
            "SER":  functions.nsplit(bb, "CB OG HG"),
            "THR":  functions.nsplit(bb, "CB OG1 HG1 CG2"),
            "VAL":  functions.nsplit(bb, "CB CG1 CG2"),
            "TRP":  functions.nsplit(bb, "CB CG CD2", "CD1 HD1 NE1 HE1 CE2", "CE3 HE3 CZ3 HZ3", "CZ2 HZ2 CH2 HH2"),
            "TYR":  functions.nsplit(bb, "CB CG CD1 HD1", "CD2 HD2 CE2 HE2", "CE1 HE1 CZ OH HH"),
            "POPE": phosphatydilethanolamine + palmitoyl1 + oleyl2,
            "DOPE": phosphatydilethanolamine + oleyl1     + oleyl2,
            "DPPE": phosphatydilethanolamine + palmitoyl1 + palmitoyl2,
            "POPG": phosphatidylglycerol     + palmitoyl1 + oleyl2,
            "DOPG": phosphatidylglycerol     + oleyl1     + oleyl2,
            "DPPG": phosphatidylglycerol     + palmitoyl1 + palmitoyl2,
            "DA": functions.nsplit("P OP1 OP2 O5' O3' O1P O2P", "C5' O4' C4'", "C3' C2' C1'", "N9 C4", "C8 N7 C5", "C6 N6 N1", "C2 N3"),
            "DG": functions.nsplit("P OP1 OP2 O5' O3' O1P O2P", "C5' O4' C4'", "C3' C2' C1'", "N9 C4", "C8 N7 C5", "C6 O6 N1", "C2 N2 N3"),
            "DC": functions.nsplit("P OP1 OP2 O5' O3' O1P O2P", "C5' O4' C4'", "C3' C2' C1'", "N1 C6", "C5 C4 N4", "N3 C2 O2"),
            "DT": functions.nsplit("P OP1 OP2 O5' O3' O1P O2P", "C5' O4' C4'", "C3' C2' C1'", "N1 C6", "C5 C4 O4 C7 C5M", "N3 C2 O2"),
            }

        # Generic names for side chain beads
        self.residue_bead_names = functions.spl("BB SC1 SC2 SC3 SC4")
        # Generic names for DNA beads
        self.residue_bead_names_dna = functions.spl("BB1 BB2 BB3 SC1 SC2 SC3 SC4")

        # This dictionary contains the bead names for all residues,
        # following the order in 'mapping'
        self.names  = {
            "POPE": "NH3 PO4 GL1 GL2 C1A C2A C3A C4A C1B C2B D3B C4B C5B".split(),
            "POPG": "GLC PO4 GL1 GL2 C1A C2A C3A C4A C1B C2B D3B C4B C5B".split()
            }
        # Add default bead names for all amino acids
        self.names.update([(i, ("BB", "SC1", "SC2", "SC3", "SC4")) for i in mapping.AA3])

        # Add the default bead names for all DNA nucleic acids
        self.names.update([(i, ("BB1", "BB2", "BB3", "SC1", "SC2", "SC3", "SC4")) for i in mapping.nucleic])

        # This dictionary allows determining four letter residue names
        # for ones specified with three letters, e.g., resulting from
        # truncation to adhere to the PDB format.
        # Each entry returns a prototypical test, given as a string,
        # and the residue name to be applied if eval(test) is True.
        # This is particularly handy to determine lipid types.
        # The test assumes there is a local or global array 'atoms'
        # containing the atom names of the residue in correct order.
        self.restest = {
            "POP": [('atoms[0] == "CA"', "POPG"),
                    ('atoms[0] == "N"',  "POPE")]
            }

        # Crude mass for weighted average. No consideration of united atoms.
        # This will probably give only minor deviations, while also giving less headache
        self.mass = {'H': 1, 'C': 12, 'N': 14, 'O': 16, 'S': 32, 'P': 31, 'M': 0}


    def finish(self):
        print "A"
        ## BACKBONE BEAD TYPE ##                                                                    
        # Dictionary of default bead types (*D)                                                     
        self.bbBeadDictD  = functions.hash(secstruc.bbss,self.bbdef)                                                             
        # Dictionary of dictionaries of types for specific residues (*S)                            
        self.bbBeadDictS  = dict([(i,functions.hash(secstruc.bbss,self.bbtyp[i])) for i in self.bbtyp.keys()])                        
        print "B"
        ## BB BOND TYPE ##                                                                          
        # Dictionary of default abond types (*D)                                                    
        self.bbBondDictD = functions.hash(secstruc.bbss,zip(self.bbldef,self.bbkb))                                                   
        # Dictionary of dictionaries for specific types (*S)                                        
        self.bbBondDictS = dict([(i,functions.hash(secstruc.bbss,zip(self.bbltyp[i],self.bbkbtyp[i]))) for i in self.bbltyp.keys()])       
        # This is tricky to read, but it gives the right bondlength/force constant
        print "C"
        ## BBB ANGLE TYPE ##                                                                        
        # Dictionary of default angle types (*D)                                                    
        self.bbAngleDictD = functions.hash(secstruc.bbss,zip(self.bbadef,self.bbka))                                                  
        # Dictionary of dictionaries for specific types (*S)                                        
        self.bbAngleDictS = dict([(i,functions.hash(secstruc.bbss,zip(self.bbatyp[i],self.bbkatyp[i]))) for i in self.bbatyp.keys()])      
        print "D"    
        ## BBBB DIHEDRAL TYPE ##                                                                    
        # Dictionary of default dihedral types (*D)                                                 
        self.bbDihedDictD = functions.hash(secstruc.bbss,zip(self.bbddef,self.bbkd,self.bbdmul))                                           
        # Dictionary of dictionaries for specific types (*S)                                        
        self.bbDihedDictS = dict([(i,functions.hash(secstruc.bbss,zip(self.bbdtyp[i],self.bbkdtyp[i]))) for i in self.bbdtyp.keys()])      
        print "E"
        
    # The following function returns the backbone bead for a given residue and                   
    # secondary structure type.                                                                 
    # 1. Look up the proper dictionary for the residue                                          
    # 2. Get the proper type from it for the secondary structure                                
    # If the residue is not in the dictionary of specials, use the default                      
    # If the secondary structure is not listed (in the residue specific                         
    # dictionary) revert to the default.                                                        
    def bbGetBead(self,r1,ss="C"):                                                                   
        return self.bbBeadDictS.get(r1,self.bbBeadDictD).get(ss,self.bbBeadDictD.get(ss))                      
    
    
    def bbGetBond(self,r,a,ss):
        # Retrieve parameters for each residue from table defined above
        b1 = self.bbBondDictS.get(r[0],self.bbBondDictD).get(ss[0],self.bbBondDictD.get(ss[0]))
        b2 = self.bbBondDictS.get(r[1],self.bbBondDictD).get(ss[1],self.bbBondDictD.get(ss[1]))
        # Determine which parameters to use for the bond
        return ( (b1[0]+b2[0])/2, min(b1[1],b2[1]) )

    def bbGetAngle(self,r,ca,ss):
        # PRO in helices is dominant
        if r[1] == "PRO" and ss[1] in "H123":
            return self.bbAngleDictS["PRO"].get(ss[1])
        else:
            # Retrieve parameters for each residue from table defined above
            a = [ self.bbAngleDictS.get(r[0],self.bbAngleDictD).get(ss[0],self.bbAngleDictD.get(ss[0])),
                  self.bbAngleDictS.get(r[1],self.bbAngleDictD).get(ss[1],self.bbAngleDictD.get(ss[1])),
                  self.bbAngleDictS.get(r[2],self.bbAngleDictD).get(ss[2],self.bbAngleDictD.get(ss[2])) ]
            # Sort according to force constant
            a.sort(key=lambda i: (i[1],i[0]))
            # This selects the set with the smallest force constant and the smallest angle
            return a[0]

    def messages(self):
        '''Prints any force-field specific logging messages.'''
