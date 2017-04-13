
import types, os
    
    
# Parameters can be defined for multiple forcefields
# We look for them within the script...
forcefields = [str(ff).split('.')[-1] for ff in globals().values() if (type(ff) == types.ClassType and hasattr(ff,"ff"))]
# ... in the local directory, ....
forcefields += [ff[:-3] for ff in os.listdir(".") if ff[-6:] == "_ff.py"]
# ... and in the GMXDATA dir.
if os.environ.has_key("GMXDATA"):
    forcefields += [ff[:-3] for ff in os.listdir(os.environ["GMXDATA"]+"/top/") if ff[-6:] == "_ff.py"]

# Martini Quotes
martiniq = [
    ("Robert Benchley",
     "Why don't you get out of that wet coat and into a dry martini?"),
    ("James Thurber",
     "One martini is all right, two is two many, three is not enough"),
    ("Philip Larkin",
     "The chromatic scale is what you use to give the effect of drinking a quinine martini and having an enema simultaneously."),
    ("William Emerson, Jr.",
     "And when that first martini hits the liver like a silver bullet, there is a sigh of contentment that can be heard in Dubuque."),
    ("Alec Waugh",
     "I am prepared to believe that a dry martini slightly impairs the palate, but think what it does for the soul."),
    ("Gerald R. Ford",
     "The three-martini lunch is the epitome of American efficiency. Where else can you get an earful, a bellyful and a snootful at the same time?"),
    ("P. G. Wodehouse",
     "He was white and shaken, like a dry martini."),
    ]



