
from .. import functions, secstruc

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
        self.setup()
        self.finish()

    def setup(self):
        raise NotImplementedError()

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
