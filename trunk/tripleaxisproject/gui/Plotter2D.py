"""
    Class that manages what is being displayed on the 2D plot
    of a SANS modeling application
"""

import wx
from PlotPanel import PlotPanel
import numpy, pylab, math
import ModelParameters
import state
from binder import BindArtist
from SlicerParameters import SlicerEvent
import detector_dialog
from copy import deepcopy

# Debug printout
import config
from config import printEVT
(DetectorParsEvent, EVT_DETECTOR)   = wx.lib.newevent.NewEvent()
(InternalEvent, EVT_INTERNAL)   = wx.lib.newevent.NewEvent()

def createDetectorParsEvent(qmax=0.1, npts=10, zmin=None, zmax=None, beam=None, info_only=False):
    return DetectorParsEvent(qmax=qmax,
                             npts=npts,
                             zmin=zmin,
                             zmax=zmax,
                             beam=beam,
                             info_only=info_only)

DEFAULT_QMAX = 0.5
DEFAULT_QSTEP = 0.01
DEFAULT_BEAM = 0.05

class ModelPanel2D(PlotPanel):
    def __init__(self, parent, id = -1, color = None,\
        dpi = None, style = wx.NO_FULL_REPAINT_ON_RESIZE, **kwargs):
        PlotPanel.__init__(self, parent, id = id, style = style, **kwargs)

        self.parent = parent
        self.qmax = DEFAULT_QMAX
        self.qstep = DEFAULT_QSTEP
        self.data = None
        self.scale = 'log'
        
        # Color map
        self.cmap_min = None
        self.cmap_max = None
        
        # Beam stop
        self.beamstop_radius = DEFAULT_BEAM
        
        # Slicer
        self.slicer = None
        self.connect = BindArtist(self.subplot.figure)
        self.axes_frozen = False
        self.slicer_z = 5

        self.x = pylab.arange(-self.qmax, self.qmax+self.qstep*0.01, self.qstep)
        self.y = pylab.arange(-self.qmax, self.qmax+self.qstep*0.01, self.qstep)

        # Label axes
        self.subplot.set_xlabel(r'$\rm{q_x}\ (A^{-1})$')
        self.subplot.set_ylabel(r'$\rm{q_y}\ (A^{-1})$')

        # Bindings
        self.parent.Bind(state.EVT_2DREFRESH, self._onEVT_2DREPLOT)
        self.parent.Bind(ModelParameters.EVT_MODEL, self.onEVT_MODEL)
        self.parent.Bind(EVT_DETECTOR, self._onEVT_DETECTOR)
        self.parent.Bind(EVT_INTERNAL, self._onEVT_INTERNAL)
        self.subplot.figure.canvas.mpl_connect('button_press_event',self._onPick)

    def _onPick(self, event):
        if event.button == 1 and not event.xdata == None and not event.ydata == None:
            # Get bin number
            i_x = 0
            i_y = 0
            for i in range(len(self.x)):
                if math.fabs(event.xdata-self.x[i])<self.qstep/2.0:
                    i_x = i
            for i in range(len(self.y)):
                if math.fabs(event.ydata-self.y[i])<self.qstep/2.0:
                    i_y = i
            # Get intensity
            if not self.data == None:
                wx.PostEvent(self.parent, 
                         config.StatusBarEvent(message="Qx = %5f; Qy = %5f; I(Qx,Qy) = %5f" % (event.xdata, event.ydata, self.data[i_x][i_y])))
            else:
                wx.PostEvent(self.parent, 
                         config.StatusBarEvent(message="Qx = %5f; Qy = %5f" % (event.xdata, event.ydata)))

    def onEVT_MODEL(self, event):
        """
            Process EVT_MODEL events
            When the model is updated, update the plot
            #TODO: do we need to process this event. Why can't 
            this object only receive the replot event?
            
            @param event: EVT_MODEL event
        """
        event.Skip()
        printEVT("Mng2D.onEVT_MODEL")
        # Fire an event instead of calling
        event.data.get_data_2d(self.qmax, self.x, self.y)
        
    def _onEVT_2DREPLOT(self, event):
        """
            Data is ready to be displayed
            @param event: data event
        """
        event.Skip()
        printEVT("ModelPanel2D.replot")
        output = event.output
        
        # Store data
        self.data = output
        if not self.data == None:
            self._plot_image()
        
    def get_corrected_data(self):
        output = deepcopy(self.data)
        self._remove_beamspot(output)
        return output
    
    def _plot_image(self, event=None):
        tmp = deepcopy(self.data)
        
        # Remove beam spot
        self._remove_beamspot(tmp)
        
        zmin = self.cmap_min
        zmax = self.cmap_max
        if self.scale == 'log':
            if zmin:
                zmin = math.log(zmin)
            if zmax:
                zmax = math.log(zmax)

        if self.scale=='log': 
            # Find minimum
            if self.cmap_min:
                tmp[tmp<=0] = min
            else:  
                min = 1
                for row in tmp:
                    for value in row:
                        if value<min and value>0:
                            min = value
                    
                tmp[tmp<=0] = min
            tmp[tmp>0] = numpy.log(tmp[tmp>0])
            
        # Extent should include half a step on each side
        # because points are calculated in the middle of the bin
        im = self.subplot.imshow(tmp, interpolation='nearest', origin='lower',
                    vmin=zmin, vmax=zmax,
                    cmap=pylab.cm.jet, extent=(-self.qmax-self.qstep/2.0, self.qmax+self.qstep/2.0,
                                               -self.qmax-self.qstep/2.0, self.qmax+self.qstep/2.0))

        self.subplot.figure.canvas.draw_idle()
        
        # Update the slicer
        if not self.slicer == None:
            self.slicer.update_and_post()

    def _remove_beamspot(self, image):
        # Get bin number
        i_x = 0
        i_y = 0
        for i_x in range(len(self.x)):
            for i_y in range(len(self.y)):
                if (self.x[i_x]*self.x[i_x]+self.y[i_y]*self.y[i_y]) \
                    < self.beamstop_radius * self.beamstop_radius:
                    image[i_x][i_y] = 0.0
      

    def freeze_axes(self):
        self.axes_frozen = True
        
    def thaw_axes(self):
        self.axes_frozen = False

    def update(self, draw=True):
        """
            Respond to changes in the model by recalculating the 
            profiles and resetting the widgets.
        """
        self.slicer.update()
        self.draw()
            

            

    def onContextMenu(self, event):
        # Slicer plot popup menu
        popupmenu = wx.Menu()
        
        popupmenu.Append(301, "&Line slicer", 
                              'Switch to a line slicer')
        popupmenu.Append(302, "&Annulus slicer", 
                              'Switch to an annulus slicer')
        popupmenu.Append(303, "&Clear slicer",
                              'Clear the current slicer')
        popupmenu.AppendSeparator()
        popupmenu.Append(304, "&Edit Parameters", "Edit detector and plot parameters")
        popupmenu.Append(313,'&Save image', 'Save image as PNG')
        popupmenu.AppendSeparator()
        popupmenu.Append(315, '&Toggle Linear/Log scale')

        wx.EVT_MENU(self, 301, self.onLineSlicer)
        wx.EVT_MENU(self, 302, self.onAnnulusSlicer)
        wx.EVT_MENU(self, 303, self.onClearSlicer)
        wx.EVT_MENU(self, 304, self._onEditDetector)
        wx.EVT_MENU(self, 313, self.onSaveImage)
        wx.EVT_MENU(self, 315, self._onToggleScale)

        pos = event.GetPosition()
        pos = self.ScreenToClient(pos)
        self.PopupMenu(popupmenu, pos)

    def _getEmptySlicerEvent(self):
        return SlicerEvent(type=None,
                           params=None,
                           obj_class=None)


    def _setSlicer(self, slicer):
        # Clear current slicer
        printEVT("Plotter2D._setSlicer %s" % slicer)
        
        if not self.slicer == None:  
            self.slicer.clear()            
            
        self.slicer_z += 1
        self.slicer = slicer(self, self.subplot, zorder=self.slicer_z)
        self.subplot.set_ylim(-self.qmax, self.qmax)
        self.subplot.set_xlim(-self.qmax, self.qmax)
        self.update()
        self.slicer.update()
        
        # Post slicer event
        event = self._getEmptySlicerEvent()
        event.type = self.slicer.__class__.__name__
        event.obj_class = self.slicer.__class__
        event.params = self.slicer.get_params()
        wx.PostEvent(self.parent, event)
        
    def _onEVT_INTERNAL(self, event):
        """
            I don't understand why Unbind followed by a Bind
            using a modified self.slicer doesn't work.
            For now, I post a clear event followed by
            a new slicer event...
        """
        self._setSlicer(event.slicer)

    def onLineSlicer(self, event):
        from LineSlicer import LineInteractor
        self.onClearSlicer(event)
        wx.PostEvent(self.parent, InternalEvent(slicer=LineInteractor))
        #self._setSlicer(LineInteractor)
        
    def onAnnulusSlicer(self, event):
        from AnnulusSlicer import AnnulusInteractor
        self.onClearSlicer(event)
        wx.PostEvent(self.parent, InternalEvent(slicer=AnnulusInteractor))
        #self._setSlicer(AnnulusInteractor)
        
    def onClearSlicer(self, event):
        if not self.slicer==None:
            self.slicer.clear()
            self.subplot.figure.canvas.draw()
            self.slicer = None
        
            # Post slicer None event
            event = self._getEmptySlicerEvent()
            wx.PostEvent(self.parent, event)
            

        
    def _onEditDetector(self, event):
        dialog = detector_dialog.DetectorDialog(None, -1, "")
        dialog.setContent(len(self.x), self.qmax, self.beamstop_radius)
        if dialog.ShowModal() == wx.ID_OK:
            evt = dialog.getContent()
            wx.PostEvent(self.parent, createDetectorParsEvent(qmax=evt.qmax,
                                                              npts=evt.npts,
                                                              zmin=evt.zmin,
                                                              zmax=evt.zmax,
                                                              beam=evt.beam))

        dialog.Destroy()

    def _onEVT_DETECTOR(self, event):
        printEVT("modelPanel2D._onEVTDetector")
        event.Skip()
        if not event.qmax == None:
            self.qmax  = event.qmax
        else:
            self.qmax= DEFAULT_QMAX
            
        if not event.npts == None:
            self.qstep = (2.0*self.qmax)/(event.npts-1)
        else:
            self.qstep = DEFAULT_QSTEP
            
        self.cmap_min = event.zmin
        self.cmap_max = event.zmax
        
        if not event.beam == None:
            self.beamstop_radius = event.beam
        else:
            self.beamstop_radius = DEFAULT_BEAM

        self.x = pylab.arange(-self.qmax, self.qmax+self.qstep*0.01, self.qstep)
        self.y = pylab.arange(-self.qmax, self.qmax+self.qstep*0.01, self.qstep)

    def _onToggleScale(self, event):
        if self.data == None: 
            return
        
        if self.scale == 'log':
            self.scale = 'linear'
        else:
            self.scale = 'log'

        self._plot_image()
        