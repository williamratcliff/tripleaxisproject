#!/usr/bin/env python

"""
Example usage of BindArtist, the MPL extension to allow canvas objects.

Note: the BindArtist api is not yet settled, so this application is likely
to change.
"""

from pylab import *
from binder import BindArtist
from bspline3 import BSpline3
        
# TODO: should profile editors subclass Artist, or maybe Widget?
class BSplinePlot:
    """Edit a parametric B-Spline profile"""
    def __init__(self, ax, spline):
        self.ax = ax
        self.fig = ax.figure
        self.spline = spline

        # Initialize callbacks
        self.callback = BindArtist(self.fig)
        self.callback.clearall()
        self.ax.cla()
        self.callback('dclick',self.ax,self.onAppend)

        # Add artists to the canvas
        [self.hspline] = self.ax.plot([],[],marker='',
                                      linestyle='-',color='green')        
        self.hpoints = []
        for x,y in self.spline:
            self._appendknot(x,y)

        # Draw the canvas
        self.draw()
    
    def _appendknot(self,x,y):
        """
        Append a knot to the list of knots on the canvas.  This had better
        maintain order within the spline if callbacks are to work properly.
        """
        [h] = self.ax.plot([x],[y],marker='s',markersize=8,
                      linestyle='',color='yellow',pickradius=5)
        self.callback('enter',h,self.onHilite)
        self.callback('leave',h,self.onHilite)
        self.callback('drag',h,self.onDrag)
        self.callback('dclick',h,self.onRemove)
        self.hpoints.append(h)
        return True
    
    def _removeknot(self, artist):
        """
        Remove the knot associated with the artist.
        """
        i = self.hpoints.index(artist)
        del self.spline[i]
        del self.hpoints[i]
        artist.remove()
 
    def draw(self):
        """
        Recompute the spline curve and show the canvas.
        """
        x,y = self.spline.sample()
        self.hspline.set_data(x,y)
        self.fig.canvas.draw_idle()

    def onHilite(self,ev):
        """
        Let the user know which point is being edited.
        """
        if ev.action == 'enter':
            ev.artist.set_color('lightblue')
        else:
            ev.artist.set_color('yellow')
        self.fig.canvas.draw_idle()
        return True

    def onDrag(self,ev):
        """
        Move the selected control point.
        """
        if ev.inaxes == self.ax:
            i = self.hpoints.index(ev.artist)
            self.spline[i] = ev.xdata, ev.ydata
            ev.artist.set_data([ev.xdata],[ev.ydata])
            self.draw()
        return True

    def onRemove(self,ev):
        """
        Remove the selected point.
        """
        # Don't delete the last one
        if len(self.spline) > 1:
            self._removeknot(ev.artist)
            self.draw()
        return True
        
    def onAppend(self,ev):
        """
        Append a new control point to the end of the spline.
        """
        self.spline.append(ev.xdata,ev.ydata)
        self._appendknot(ev.xdata,ev.ydata)
        self.draw()
        return True
        
def demo():
    import numpy
    Ix = numpy.array([1,2,3,4,5,6],'f')
    Iy = numpy.array([1,2,3,2,1,2],'f')
    model = BSplinePlot(subplot(111),BSpline3(Ix,Iy))
    show()

if __name__ == "__main__": demo()
