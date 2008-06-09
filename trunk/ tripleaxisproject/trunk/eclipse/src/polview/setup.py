#!/usr/bin/env python

#
# The setup to create a single exe file.
#

import os, sys

from distutils.core import setup

from distutils.filelist import findall

import matplotlib

import py2exe

    
class Target:
    def __init__(self, **kw):
        self.__dict__.update(kw)
        # for the versioninfo resources
        self.version = "0.1"
        self.company_name = "NIST"
        self.copyright = "Public Domain"
        self.name = "PolApp"
        
#
# Adapted from http://www.py2exe.org/index.cgi/MatPlotLib
# to use the MatPlotLib.
#
matplotlibdatadir = matplotlib.get_data_path()
matplotlibdata = findall(matplotlibdatadir)
matplotlibdata_files = []

for f in matplotlibdata:
    dirname = os.path.join('matplotlibdata', f[len(matplotlibdatadir)+1:])
    matplotlibdata_files.append((os.path.split(dirname)[0], [f]))

#
# packages
#
#packages = [
#    'matplotlib', 'pytz'
#   ]

packages = [
    'matplotlib'
    ]



includes = []
excludes = ['OpenGL','enthought','BaseHTTPServer'] 

dll_excludes = [
    'libgdk_pixbuf-2.0-0.dll', 
    'libgobject-2.0-0.dll',
    'libgdk-win32-2.0-0.dll',
    ]

## This is the client for PARK: run on wx
target_wx_client = Target(
    description = 'Polarization Correcter',
    script = 'PolApp.py',
    #other_resources = [(RT_MANIFEST, 1, manifest_template % dict(prog="AppJob"))],
    #icon_resources = [(1, "images/ball.ico")],
    dest_base = "PolApp"
    )



setup(
    windows=[target_wx_client],
    console=[],
    
    options={
        'py2exe': {
            'dll_excludes': dll_excludes,
            'packages' : packages,
            'includes':includes,
            'excludes':excludes,
            "compressed": 1,
            "optimize": 0,
            "bundle_files":3,
            'typelibs' : [('{EAB22AC0-30C1-11CF-A7EB-0000C05BAE0B}', 0, 1, 1)],
            },
    },
    data_files=matplotlibdata_files
    
    # Do something like this to add images to the distribution
    #data_files=[ ("prog",["kategorien.xml",])]
)
