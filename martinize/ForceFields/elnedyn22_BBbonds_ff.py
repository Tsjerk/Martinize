from .forcefield import Forcefield
from .. import secstruc,functions,IO 
from ..functions import distance2, cos_angle
from ..IO import d2r

################################
## 6 # FORCE FIELD PARAMETERS ##  -> @FF <-
################################

class elnedyn22(Forcefield):
    '''The forcefield has been implemented with some changes compared to the published parameters:
    - Backbone-Backbone bonds are constraints in stead of strong bonds.
    - Trp has an extra constrain added to the sidechain
    - The Backbone-Sidechain bonds with high force constants are replaced by constraints except for Trp and His.
    '''
    ff = True
    name = 'elnedyn22_BBbonds'

    def setup(self):
        
        # Charged types:
        self.charges = {"Qd":1, "Qa":-1, "SQd":1, "SQa":-1, "RQd":1, "AQa":-1}                                                           #@#
        
        
        #----+---------------------+
        ## A | BACKBONE PARAMETERS |
        #----+---------------------+
        #
        # bbss  lists the one letter secondary structure code
        # bbdef lists the corresponding default backbone beads
        # bbtyp lists the corresponding residue specific backbone beads
        #
        # bbd   lists the structure specific backbone bond lengths
        # bbkb  lists the corresponding bond force constants
        #
        # bba   lists the structure specific angles
        # bbka  lists the corresponding angle force constants
        #
        # bbd   lists the structure specific dihedral angles
        # bbkd  lists the corresponding force constants
        #
        # -=NOTE=- 
        #  if the secondary structure types differ between bonded atoms
        #  the bond is assigned the lowest corresponding force constant 
        #
        # -=NOTE=-
        # if proline is anywhere in the helix, the BBB angle changes for 
        # all residues
        #
        
        ###############################################################################################
        ## BEADS ##                                                                          #                 
        #                               F     E     H     1     2     3     T     S     C    # secstruc one letter   
        self.bbdef    =     functions.spl(" N0   Nda    N0    Nd    Na   Nda   Nda    P5    P5")  # Default beads   #@#
        self.bbtyp    = {                                                                    #                 #@#
                     "ALA": functions.spl(" C5    N0    C5    N0    N0    N0    N0    P4    P4"), # ALA specific    #@#
                     "PRO": functions.spl(" C5    N0    C5    N0    Na    N0    N0    P4    P4"), # PRO specific    #@#
                     "HYP": functions.spl(" C5    N0    C5    N0    N0    N0    N0    P4    P4")  # HYP specific    #@#
        }                                                                                    #                 #@#
        ## BONDS ##                                                                          #                 
        self.bbldef   =             (.365, .350, .350, .350, .350, .350, .350, .350, .350)   # BB bond lengths #@#
        self.bbkb     =             (1250, 1250, 1250, 1250, 1250, 1250,  500,  400,  400)   # BB bond kB      #@#
        self.bbltyp   = {}                                                                   #                 #@#
        self.bbkbtyp  = {}                                                                   #                 #@#
        ## ANGLES ##                                                                         #                 
        self.bbadef   =             (119.2, 134,   96,   96,   96,   96,  100,  130,  127)   # BBB angles      #@#
        self.bbka     =             ( 150,   25,  700,  700,  700,  700,   25,   25,   25)   # BBB angle kB    #@#
        self.bbatyp   = {                                                                    #                 #@#
                    "PRO":          ( 119.2,134,   98,   98,   98,   98,  100,  130,  127),  # PRO specific    #@#
                    "HYP":          ( 119.2,134,   98,   98,   98,   98,  100,  130,  127)   # PRO specific    #@#
        }                                                                                    #                 #@#
        self.bbkatyp  = {                                                                    #                 #@#
                    "PRO":          ( 150,   25,  100,  100,  100,  100,   25,   25,   25),  # PRO specific    #@#
                    "HYP":          ( 150,   25,  100,  100,  100,  100,   25,   25,   25)   # PRO specific    #@#
        }                                                                                    #                 #@#
        ## DIHEDRALS ##                                                                      #                 
        self.bbddef   =             (90.7,    0, -120, -120, -120, -120)                     # BBBB dihedrals  #@#
        self.bbkd     =             ( 100,   10,  400,  400,  400,  400)                     # BBBB kB         #@#
        self.bbdmul   =             (   1,    1,    1,    1,    1,    1)                     # BBBB mltplcty   #@#
        self.bbdtyp   = {}                                                                   #                 #@#
        self.bbkdtyp  = {}                                                                   #                 #@#
                                                                                             #                 
        ###############################################################################################               
        
        # Some Forcefields use the Ca position to position the BB-bead (me like!)
        self.ca2bb = True 
        
        # BBS angle, equal for all ss types                                                         
        # Connects BB(i-1),BB(i),SC(i), except for first residue: BB(i+1),BB(i),SC(i)               
        #                      ANGLE   Ka                                                                
        self.bbsangle =      [   100,  25]                                                          #@#
        
        # Bonds for extended structures (more stable than using dihedrals)                          
        #               LENGTH FORCE                                                                
        self.ebonds   = {                                                                           #@#
               'short': [ .640, 2500],                                                              #@#
               'long' : [ .970, 2500]                                                               #@#
        }                                                                                           #@#
        
        
        #----+-----------------------+
        ## B | SIDE CHAIN PARAMETERS |
        #----+-----------------------+
        
        # Sidechain parameters for Elnedyn. (read from cg-2.1.dat). 
        # For HIS the order of bonds is changed and a bond with fc=0 is added.
        # In the elnedyn2, TRP has an extra, cross-ring constraint
        self.sidechains = {
        #RES#   BEADS                      BONDS                                                                    ANGLES                          DIHEDRALS
        'TRP': [functions.spl("SC4 SNd SC5 SC5"), [(0.255,73000), (0.220,None), (0.250,None), (0.280,None), (0.255,None), (0.35454,None)], [(142,30), (143,20), (104,50)], [(180,200)]],
        'TYR': [functions.spl("SC4 SC4 SP1"),     [(0.335, 6000), (0.335,6000), (0.240,None), (0.310,None), (0.310,None)], [(70,100), (130, 50)]],
        'PHE': [functions.spl("SC5 SC5 SC5"),     [(0.340, 7500), (0.340,7500), (0.240,None), (0.240,None), (0.240,None)], [(70,100), (125,100)]],
        'HIS': [functions.spl("SC4 SP1 SP1"),     [(0.195, None), (0.193,None), (0.295,None), (0.216,None)],               [(135,100),(115, 50)]],
        'HIH': [functions.spl("SC4 SP1 SP1"),     [(0.195, None), (0.193,None), (0.295,None), (0.216,None)],               [(135,100),(115, 50)]],
        'ARG': [functions.spl("N0 Qd"),           [(0.250,12500), (0.350,6200)],                                           [(150,15)]],
        'LYS': [functions.spl("C3 Qd"),           [(0.250,12500), (0.300,9700)],                                           [(150,20)]],
        'CYS': [functions.spl("C5"),              [(0.240, None)]],
        'ASP': [functions.spl("Qa"),              [(0.255, None)]],
        'GLU': [functions.spl("Qa"),              [(0.310, 2500)]],
        'ILE': [functions.spl("AC1"),              [(0.225,13250)]],
        'LEU': [functions.spl("AC1"),              [(0.265, None)]],
        'MET': [functions.spl("C5"),              [(0.310, 2800)]],
        'ASN': [functions.spl("P5"),              [(0.250, None)]],
        'PRO': [functions.spl("C3"),              [(0.190, None)]],
        'GLN': [functions.spl("P4"),              [(0.300, 2400)]],
        'SER': [functions.spl("P1"),              [(0.195, None)]],
        'THR': [functions.spl("P1"),              [(0.195, None)]],
        'VAL': [functions.spl("AC2"),              [(0.200, None)]],
        'GLY': [],
        'ALA': [],
        }
        
        # Not all (eg Elnedyn) forcefields use backbone-backbone-sidechain angles and BBBB-dihedrals.
        self.UseBBSAngles        = False 
        self.UseBBBBDihedrals    = False

        # Martini 2.2p has polar and charged residues with seperate charges.
        self.polar   = []
        self.charged = []

        # If masses or charged diverge from standard (45/72 and -/+1) they are defined here.
        self.mass_charge = {
        #RES   MAsecstruc               CHARGE
        }

        # Defines the connectivity between between beads
        # Connectivity records for Elnedyn (read from cg-2.1.dat). 
        # For HIS the order of bonds is changed and a bond with fc=0 is added.
        self.connectivity = {
        #RES       BONDS                                             ANGLES                            DIHEDRALS       V-SITE
        "TRP":     [[(0, 1), (1, 2), (2, 4), (4, 3), (3, 1), (1, 4)],[(0, 1, 2), (0, 1, 4), (0, 1, 3)],[(1, 2, 3, 4)]],
        "TYR":     [[(0, 1), (0, 2), (1, 2), (1, 3), (2, 3)],        [(0, 1, 2), (0, 1, 3)]],
        "PHE":     [[(0, 1), (0, 2), (1, 2), (1, 3), (2, 3)],        [(0, 1, 2), (0, 1, 3)]],
        "HIS":     [[(0, 1), (1, 2), (1, 3), (2, 3)],        [(0, 1, 2), (0, 1, 3)]],
        "HIH":     [[(0, 1), (1, 2), (1, 3), (2, 3)],        [(0, 1, 2), (0, 1, 3)]],
        "GLN":     [[(0,1)]],
        "ASN":     [[(0,1)]],
        "SER":     [[(0,1)]],
        "THR":     [[(0,1)]],
        "ARG":     [[(0,1),(1,2)],                         [(0,1,2)]],
        "LYS":     [[(0,1),(1,2)],                         [(0,1,2)]],
        "ASP":     [[(0,1)]],
        "GLU":     [[(0,1)]],
        "CYS":     [[(0,1)]],
        "ILE":     [[(0,1)]],
        "LEU":     [[(0,1)]],
        "MET":     [[(0,1)]],
        "PRO":     [[(0,1)]],
        "HYP":     [[(0,1)]],
        "VAL":     [[(0,1)]],
        "ALA":     [],
        "GLY":     [],
        }
       
        #----+----------------+
        ## C | SPECIAL BONDS  |
        #----+----------------+
        
        self.special = {
            # Used for sulfur bridges
            # ATOM 1         ATOM 2          BOND LENGTH   FORCE CONSTANT
            (("SC1","CYS"), ("SC1","CYS")):     (0.24,         None),
            }
       
        # By default use an elastic network
        self.ElasticNetwork = True 

        # Elastic networks bond shouldn't lead to exclusions (type 6) 
        # But Elnedyn has been parametrized with type 1.
        self.EBondType = 1
        
        self.ss = secstruc

    def bbGetBond(self,r,ca,ss):
        import math
        # The 150000 forceconstant gave an error message, turning to constraints would be better.
        return ( math.sqrt(distance2(ca[0],ca[1]))/10., None   )
    
    def bbGetAngle(self,r,ca,ss):
        import math
        # Elnedyn takes angles from structure, with fc=40
        return (math.acos(cos_angle([i-j for i,j in zip(ca[0],ca[1])],[i-j for i,j in zip(ca[2],ca[1])]))/d2r, 40)

    def messages(self):
        '''Prints any force-field specific logging messages.'''
        import logging
        logging.info('The elnedyn forcefield has been implemented with some changes compared to the published parameters:')
        #logging.info('- Backbone-Backbone bonds are constraints in stead of high force constant bonds.')
        logging.info('- Backbone-Backbone bonds use high force constant bonds instead of constraints.')
        logging.info('- Trp has an extra constrain added to the sidechain.')
        logging.info('- The Backbone sidechain bonds with high force constants are replaced by constraints except for Trp and His.')
        logging.info('- Cysteine bonds are 0.24 nm constraints, instead of the published 0.39nm/5000kJ/mol.')
        logging.warning('Elnedyn topologies might not give numerical stable simulations with a 20fs timestep.')
        logging.warning('This can be solved by setting all S-type bead masses to 72amu.')
        pass
