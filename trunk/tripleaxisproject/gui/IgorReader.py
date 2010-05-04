import wx
from PlotPanel import PlotPanel
from copy import deepcopy
import os
import pylab, math, numpy
from binder import BindArtist
import config
import SlicerParameters
import Plotter2D
           
class DataFrame(wx.Frame):
        
    def __init__(self, parent, id, title):
        wx.Frame.__init__(self, parent, id, title, wx.DefaultPosition, wx.Size(370,370))
        
        self.parent = parent
        self.SetIcon(wx.Icon('images/ball.ico', wx.BITMAP_TYPE_ICO))
        
        self.two_d_plot = DataPanel2D(self, -1, style=wx.RAISED_BORDER)
        self._setup_layout()
        
        # Bindings
        self.parent.Bind(SlicerParameters.EVT_SLICER, self.two_d_plot._onEVT_SLICER)
        wx.EVT_CLOSE(self, self._onClose)

    def _onClose(self, event):
        if not self.two_d_plot.slicer == None:
            self.two_d_plot.slicer.clear()
        self.parent.Unbind(SlicerParameters.EVT_SLICER)
        self.two_d_plot.parent.Unbind(Plotter2D.EVT_DETECTOR)
        self.Destroy()
        
    def _setup_layout(self):
        
        
        # Two columns of panels
        sizer = wx.GridBagSizer(0,0)
        
        
        # Add panels
        sizer.Add(self.two_d_plot, (0,0), (1,1), flag=wx.EXPAND | wx.ALL, border=0)
        sizer.SetItemMinSize(self.two_d_plot, 250, 250)
        sizer.AddGrowableRow(0)
        sizer.AddGrowableCol(0)
        
        
        self.SetSizer(sizer)
        self.Centre()
        
    def read(self, filename):
        reader = DataReader(filename)
        Z, xmin, xmax, ymin, ymax = reader.read()
        self.two_d_plot.x = reader.x
        self.two_d_plot.y = reader.y
        self.two_d_plot.plot_data(Z, xmin, xmax, ymin, ymax)
        
        

           
class DataPanel2D(PlotPanel):
    def __init__(self, parent, id = -1, color = None,\
        dpi = None, style = wx.NO_FULL_REPAINT_ON_RESIZE, **kwargs):
        PlotPanel.__init__(self, parent, id = id, style = style, **kwargs)

        self.parent = parent.parent
        self.data = None
        self.scale = 'log'
        self.qmax = 0.5
        
        # Extent of the plit
        self.xmin = -0.1
        self.xmax = 0.1
        self.ymin = -0.1
        self.ymax = 0.1

        # Color map
        self.cmap_min = None
        self.cmap_max = None
        
        # Slicer
        self.slicer = None
        self.connect = BindArtist(self.subplot.figure)
        self.axes_frozen = False
        self.slicer_z = 5
        
        self.x = []
        self.y = []

        # Label axes
        self.subplot.set_xlabel(r'$\rm{q_x}\ (A^{-1})$')
        self.subplot.set_ylabel(r'$\rm{q_y}\ (A^{-1})$')
        
        self.parent.Bind(Plotter2D.EVT_DETECTOR, self._onEVT_DETECTOR)
        self.subplot.figure.canvas.mpl_connect('button_press_event',self._onPick)
        

    def _onPick(self, event):
        if event.button == 1 and not event.xdata == None and not event.ydata == None:
            # Get bin number
            i_x = 0
            i_y = 0
            
            for i in range(len(self.x)):
                if math.fabs(event.xdata-self.x[i]) < math.fabs(event.xdata-self.x[i_x]):
                    i_x = i
            for i in range(len(self.y)):
                if math.fabs(event.ydata-self.y[i]) < math.fabs(event.ydata-self.y[i_y]):
                    i_y = i
            # Get intensity
            if not self.data == None:
                wx.PostEvent(self.parent, 
                         config.StatusBarEvent(message="Qx = %5f; Qy = %5f; I(Qx,Qy) = %5f" % (event.xdata, event.ydata, self.data[i_x][i_y])))
            else:
                wx.PostEvent(self.parent, 
                         config.StatusBarEvent(message="Qx = %5f; Qy = %5f" % (event.xdata, event.ydata)))


    def _onEVT_DETECTOR(self, event):
        config.printEVT("DataPanel2D._onEVTDetector")
        event.Skip()
            
        self.cmap_min = event.zmin
        self.cmap_max = event.zmax

        self.plot_data(self.data, self.xmin, self.xmax, self.ymin, self.ymax)

    def _onEVT_SLICER(self, event):
        event.Skip()
        config.printEVT("DataPanel2D._onEVT_SLICER")
        # Clear current slicer
        if not self.slicer == None:                
            self.slicer.clear()
            
        self.subplot.figure.canvas.draw()
        self.slicer = None
        if event.obj_class == None:
            return
            
        self.slicer_z += 1
        self.slicer = event.obj_class(self, self.subplot, zorder=self.slicer_z)
        self.subplot.set_ylim(self.ymin, self.ymax)
        self.subplot.set_xlim(self.xmin, self.xmax)
        self.slicer.set_params(event.params)
        self.slicer.update_and_post()
        self.update()
        
        
    def plot_data(self, matrix, qx_min, qx_max, qy_min, qy_max):
        
        self.xmin = qx_min
        self.xmax = qx_max
        self.ymin = qy_min
        self.ymax = qy_max
        self.data = matrix
        
        self._plot_image()
        
    def _plot_image(self):        
        
        output = deepcopy(self.data)

        zmin = self.cmap_min
        zmax = self.cmap_max
        if self.scale == 'log':
            if zmin:
                zmin = math.log(zmin)
            if zmax:
                zmax = math.log(zmax)

            output[output>0] = numpy.log(output[output>0])
        
        im = self.subplot.imshow(output, interpolation='nearest', origin='lower',
                    vmin=zmin, vmax=zmax,
                    cmap=pylab.cm.jet, extent=(self.xmin, self.xmax,
                                               self.ymin, self.ymax))
        #self.Update()
        self.Refresh()
        #self.subplot.figure.canvas.draw_idle()
        
        # Update the slicer
        if not self.slicer == None:
            self.slicer.update_and_post()
        
    def get_corrected_data(self):
        return self.data
    
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
            

    def _onToggleScale(self, event):
        if self.data == None: 
            return
        
        if self.scale == 'log':
            self.scale = 'linear'
        else:
            self.scale = 'log'

        self._plot_image()

    def onContextMenu(self, event):
        # Slicer plot popup menu
        popupmenu = wx.Menu()
        
        popupmenu.Append(313,'&Save image', 'Save image as PNG')
        popupmenu.AppendSeparator()
        popupmenu.Append(315, '&Toggle Linear/Log scale')

        wx.EVT_MENU(self, 313, self.onSaveImage)
        wx.EVT_MENU(self, 315, self._onToggleScale)

        pos = event.GetPosition()
        pos = self.ScreenToClient(pos)
        self.PopupMenu(popupmenu, pos)
  
                
class DataReader:
    """ Simple data reader for Igor data files """
    
    def __init__(self, filename):
        """ Init
            @param filename: Name of Igor data file to read
        """
        self.file = filename
        self.x = []
        self.y = []
        self.image = []
        
    def read(self):
        """ Read file """
        # Check if the file is there
        if not os.path.isfile(self.file):
            raise ValueError, \
            "Specified file %s is not a regular file" % self.file
        
        # Read file
        f = open(self.file,'r')
        buf = f.read()
        
        # Get content
        dataStarted = False
        
        
        lines = buf.split('\n')
        itot = 0
        self.x = []
        self.y = []
        self.image = []
        
        ncounts = 0
        
        #x = pylab.arange(0, 128, 1)
        #y = pylab.arange(0, 128, 1)
        x = pylab.arange(-.5, .5, 1.0/128)
        y = pylab.arange(-.5, .5, 1.0/128)
        X, Y = pylab.meshgrid(x, y)
        Z = deepcopy(X)
        
        xmin = None
        xmax = None
        ymin = None
        ymax = None
        
        i_x = 0
        i_y = -1
        
        isInfo = False
        isCenter = False
        for line in lines:
            
            # Find setup info line
            if isInfo:
                isInfo = False
                line_toks = line.split()
                # Wavelength in Angstrom
                wavelength = float(line_toks[1])
                # Distance in meters
                distance = float(line_toks[3])
                
            if line.count("LAMBDA")>0:
                isInfo = True
                
            # Find center info line
            if isCenter:
                isCenter = False                
                line_toks = line.split()
                # Center in bin number
                center_x = float(line_toks[0])
                center_y = float(line_toks[1])

            if line.count("BCENT")>0:
                isCenter = True
                
        
            # Find data start
            if line.count("***")>0:
                dataStarted = True
                
                # Check that we have all the info
                if wavelength == None \
                    or distance == None \
                    or center_x == None \
                    or center_y == None:
                    raise ValueError, "Missing information in data file"
                
            if dataStarted == True:

                
                try:
                    value = float(line)
                except:
                    continue
                
                # Get bin number
                if math.fmod(itot, 128)==0:
                    i_x = 0
                    i_y += 1
                else:
                    i_x += 1
                    
                Z[i_y][i_x] = value
                ncounts += 1 
                
                # Det 640 x 640 mm
                # Q = 4pi/lambda sin(theta/2)
                # Bin size is 0.5 cm
                theta = (i_x-center_x+1)*0.5 / distance / 100.0
                qx = 4.0*math.pi/wavelength * math.sin(theta/2.0)
                if xmin==None or qx<xmin:
                    xmin = qx
                if xmax==None or qx>xmax:
                    xmax = qx
                
                theta = (i_y-center_y+1)*0.5 / distance / 100.0
                qy = 4.0*math.pi/wavelength * math.sin(theta/2.0)
                if ymin==None or qy<ymin:
                    ymin = qy
                if ymax==None or qy>ymax:
                    ymax = qy
                
                
                if not qx in self.x:
                    self.x.append(qx)
                if not qy in self.y:
                    self.y.append(qy)
                
                itot += 1
                  
                  
        theta = 0.25 / distance / 100.0
        xstep = 4.0*math.pi/wavelength * math.sin(theta/2.0)
        
        theta = 0.25 / distance / 100.0
        ystep = 4.0*math.pi/wavelength * math.sin(theta/2.0)
        
        # Store q max 
        if xmax>ymax:
            self.qmax = xmax
        else:
            self.qmax = ymax
        
        print xmin, xmax, ymin, ymax, xstep, ystep
        print center_x, center_y
        print self.x
        print self.y
        print len(self.x)
        print len(self.y)
  
        print "Read %g points from file %s" % (ncounts, self.file)
        return Z, xmin-xstep/2.0, xmax+xstep/2.0, ymin-ystep/2.0, ymax+ystep/2.0
     
class ViewApp(wx.App):
    def OnInit(self):
        frame = DataFrame(None, -1, 'SANS Viewer')  
        frame.read("MAR07259_strained.ASC")  
        frame.Show(True)
        self.SetTopWindow(frame)
        return True
        

if __name__ == "__main__": 
    app = ViewApp(0)
    app.MainLoop()        
            
   