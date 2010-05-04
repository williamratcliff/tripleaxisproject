from distutils.core import setup
import py2exe
import glob


opts = {
    'py2exe': {'excludes': ['BaseHTTPServer','enthought'],
               "includes" : ["matplotlib.backends",  "matplotlib.backends.backend_qt4agg",
                               "matplotlib.figure","pylab", "numpy", "matplotlib.numerix.fft",
                               "matplotlib.numerix.linear_algebra", "matplotlib.numerix.random_array",
                               "matplotlib.backends.backend_tkagg"],
                'dll_excludes': ['libgdk-win32-2.0-0.dll','_gtkagg', '_tkagg', '_agg2', '_cairo', 
                                 '_cocoaagg','_fltkagg', '_gtk', '_gtkcairo', 
                                 'libgobject-2.0-0.dll']
              }
       }


data_files = [(r'mpl-data', glob.glob(r'C:\Python25\Lib\site-packages\matplotlib\mpl-data\*.*')),
                    # Because matplotlibrc does not have an extension, glob does not find it (at least I think that's why)
                    # So add it manually here:
                  (r'mpl-data', [r'C:\Python25\Lib\site-packages\matplotlib\mpl-data\matplotlibrc']),
                  (r'mpl-data\images',glob.glob(r'C:\Python25\Lib\site-packages\matplotlib\mpl-data\images\*.*')),
                  (r'mpl-data\fonts',glob.glob(r'C:\Python25\Lib\site-packages\matplotlib\mpl-data\fonts\*.*'))]

setup(windows=[{"script" : "PolApp.py"}], options=opts,   data_files=data_files)