
import importlib

_QB = importlib.import_module('.QuESTBase', package = 'QuESTPy')

if _QB.QuESTLib is not None:
    __all__ = ['QuESTTypes','QuESTFunc']
    from .QuESTTypes import *
    from .QuESTFunc import *
else:
    __all__ = ['QuESTBase']
del _QB
