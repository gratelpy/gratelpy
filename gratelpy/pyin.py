""" Stuff needed to compile GraTeLPy with PyInstaller """

import os
import sys

def resource_path(base_path, relative_path):
    """ Get absolute path to resource, works for dev and for PyInstaller 
    (adapted from http://stackoverflow.com/a/13790741/2042696)
    """
    try:
        # PyInstaller creates a temp folder and stores path in _MEIPASS
        base_path = sys._MEIPASS
    except Exception:
        pass

    return os.path.join(base_path, relative_path)
