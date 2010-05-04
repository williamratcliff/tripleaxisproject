import ctypes, numpy

Embed = numpy.dtype([('n',numpy.intc),('k',numpy.intc),
                     ('A',float,(4,4)),('B',float,(4,4))])
a = Embed()
lib = numpy.ctypeslib.load_library('struct')
#lib = ctypes.cdll['struct.dyld']
lib.call(a.ctypes.data)
print a['A'],a['n']
