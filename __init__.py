from importlib import import_module
from os.path import dirname, basename, isfile, join, isdir
from os import listdir
import glob

# add all the python files

gpath = dirname(__file__)
modules = glob.glob(join(gpath, "*.py"))
__all__ = [ basename(f)[:-3] for f in modules if isfile(f) and not f.endswith('__init__.py')]

# add all the folders that has python files

modules = listdir(gpath)
folders = [ basename(f) for f in modules if isdir(join(gpath, f))  and not ("__" in f)]
for folder in folders:
    try:
        pm = import_module(folder)
        __all__ += [folder]
    except:
        pass
