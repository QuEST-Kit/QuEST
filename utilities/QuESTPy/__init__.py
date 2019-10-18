"""
Module to import and set up QuEST's Python interface environment.

If QuESTLibdir is set correctly, simply importing the module will perform all the setup.

If not, only QuESTBase will be imported to allow manual import of libQuEST
"""
import importlib

_QB = importlib.import_module('.QuESTBase', package='QuESTPy')

if _QB.QuESTLib is not None:
    __all__ = ['QuESTTypes', 'QuESTFunc']
else:
    __all__ = ['QuESTBase']
del _QB
