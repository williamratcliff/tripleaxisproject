"""
This demo demonstrates how to embed a matplotlib (mpl) plot 
into a PyQt4 GUI application, including:

* Using the navigation toolbar
* Adding data to the plot
* Dynamically modifying the plot's properties
* Processing mpl events
* Saving the plot to a file from a menu

The main goal is to serve as a basis for developing rich PyQt GUI
applications featuring mpl plots (using the mpl OO API).

Eli Bendersky (eliben@gmail.com)
License: this code is in the public domain
Last modified: 19.01.2009
"""
global iter
iter = 0
import sys, os, random, copy
from PyQt4.QtCore import *
from PyQt4.QtGui import *

from math import log10,pow 
import matplotlib, numpy,math
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar
from matplotlib.figure import Figure
import d10reader2 as d10reader
import matplotlib.colors
from matplotlib.font_manager import FontProperties
from matplotlib.colors import Normalize,LogNorm
import ticker
import dragable_rectangle
from matplotlib.patches import Rectangle
from filetable import DataModel, DataDelegate,MainForm

#(FILE_NUMBER, H, K, L, TTHETA, OMEGA) = range(6)

#from cmapmenu import CMapMenu


def _rescale(lo,hi,step,pt=None,bal=None,scale='linear'):
    """
    Rescale (lo,hi) by step, returning the new (lo,hi)
    The scaling is centered on pt, with positive values of step
    driving lo/hi away from pt and negative values pulling them in.
    If bal is given instead of point, it is already in [0,1] coordinates.

    This is a helper function for step-based zooming.
    """
    # Convert values into the correct scale for a linear transformation
    # TODO: use proper scale transformers
    if scale=='log':
        lo,hi = log10(lo),log10(hi)
        if pt is not None: pt = log10(pt)

    # Compute delta from axis range * %, or 1-% if persent is negative
    if step > 0:
        delta = float(hi-lo)*step/100
    else:
        delta = float(hi-lo)*step/(100-step)

    # Add scale factor proportionally to the lo and hi values, preserving the
    # point under the mouse
    if bal is None:
        bal = float(pt-lo)/(hi-lo)
    lo = lo - bal*delta
    hi = hi + (1-bal)*delta

    # Convert transformed values back to the original scale
    if scale=='log':
        lo,hi = pow(10.,lo),pow(10.,hi)

    return (lo,hi)

def bbox_union(bboxes):
    """
    Return a Bbox that contains all of the given bboxes.
    """
    from matplotlib.transforms import Bbox
    if len(bboxes) == 1:
        return bboxes[0]

    x0 = numpy.inf
    y0 = numpy.inf
    x1 = -numpy.inf
    y1 = -numpy.inf

    for bbox in bboxes:
        xs = bbox.intervalx
        ys = bbox.intervaly
        x0 = min(x0, numpy.min(xs))
        y0 = min(y0, numpy.min(ys))
        x1 = max(x1, numpy.max(xs))
        y1 = max(y1, numpy.max(ys))

    return Bbox([[x0, y0], [x1, y1]])



class AppForm(QMainWindow):
    def __init__(self, parent=None):
        QMainWindow.__init__(self, parent)
        self.setWindowTitle('Demo: PyQt with matplotlib')

        self.create_menu()
        self.create_main_frame()
        self.create_status_bar()

        self.textbox.setText('12 20 12 20')
        self.on_draw()
        self.on_draw_tth()

    def save_plot(self):
        file_choices = "PNG (*.png)|*.png"
        
        path = unicode(QFileDialog.getSaveFileName(self, 
                        'Save file', '', 
                        file_choices))
        if path:
            self.canvas.print_figure(path, dpi=self.dpi)
            self.statusBar().showMessage('Saved to %s' % path, 2000)
    
    def on_about(self):
        msg = """ A demo of using PyQt with matplotlib:
        
         * Use the matplotlib navigation bar
         * Add values to the text box and press Enter (or click "Draw")
         * Show or hide the grid
         * Drag the slider to modify the width of the bars
         * Save the plot to a file using the File menu
         * Click on a bar to receive an informative message
        """
        QMessageBox.about(self, "About the demo", msg.strip())
    
    def on_pick(self, event):
        # The event received here is of the type
        # matplotlib.backend_bases.PickEvent
        #
        # It carries lots of information, of which we're using
        # only a small amount here.
        # 
        box_points = event.artist.get_bbox().get_points()
        msg = "You've clicked on a bar with coords:\n %s" % box_points
        
        QMessageBox.information(self, "Click!", msg)
        
    def save_rect(self):
        self.rx=self.dr.rect.get_x()
        self.ry=self.dr.rect.get_y()
        self.rwidth=self.dr.rect.get_width()
        self.rheight=self.dr.rect.get_height()
    
    def on_draw_tth(self):
        """ Redraws the figure
        """
        #str = unicode(self.textbox.text())
        #self.data = map(int, str.split())
        # clear the axes and redraw the plot anew
        #
        
        ax=self.tth_axes
        #self.save_rect()
        ax.clear()
        
        for curr_ax in self.axes:
            curr_ax.grid(alpha=0.8, visible=self.grid_cb.isChecked())
        currentframe=self.spinbox.value()
        #snap=self.data.frames[currentframe]+1 #To deal with annoyance of log scale!!!
        
        dataset=numpy.array(self.data.frames)
        snap=numpy.zeros((dataset.shape[1],dataset.shape[0]),'Float64')
        for i in range(dataset.shape[2]):
            snap=snap+dataset[:,:,i].T
            
        snap=snap+1  #to deal with annoyance of log scale!
        
        #ax.imshow(im)
        #self.canvas.draw()
        if self.vscale == 'linear':
            norm = Normalize(numpy.min(snap),numpy.max(snap))
        else:
            norm = LogNorm(numpy.min(snap),numpy.max(snap))
        try:
            im = ax.pcolorfast(snap,norm=norm)
        except AttributeError:
            im = ax.pcolormesh(snap,shading='flat',norm=norm)
        self.im_tth = im
        #self.colormapper.callbacksSM.connect('changed',self.on_colormap_changed)
        #self.vmin = min(self.vmin, numpy.min(snap))
        #self.vmax = max(self.vmax, numpy.max(snap))
        #self.colormapper.set_clim(vmin=self.vmin,vmax=self.vmax)
        #self.colorbar.set_clim(vmin=self.vmin,vmax=self.vmax)
        #self.autoaxes()
        #self.attach_rectangle()
        self.tth_canvas.draw_idle()
        #self.onIntegrate()
        #self.on_int_draw()
        
    def on_draw(self):
        """ Redraws the figure
        """
        #str = unicode(self.textbox.text())
        #self.data = map(int, str.split())
        # clear the axes and redraw the plot anew
        #
        
        ax=self.axes[0]
        self.save_rect()
        ax.clear()
        
        for curr_ax in self.axes:
            curr_ax.grid(alpha=0.8, visible=self.grid_cb.isChecked())
        currentframe=self.spinbox.value()
        snap=self.data.frames[currentframe]+1 #To deal with annoyance of log scale!!!
        
        #if self.vscale=='log':
        #    snap=copy.deepcopy(numpy.log10(snap))
        
        #ax.imshow(im)
        #self.canvas.draw()
        norm = self.colormapper.norm
        try:
            im = ax.pcolorfast(snap,norm=norm)
        except AttributeError:
            im = ax.pcolormesh(snap,shading='flat',norm=norm)
        self.im = im
        def on_changed(m):
            print "changed",m
            self.colorbar.set_cmap(m.get_cmap())
            self.colorbar.set_clim(m.get_clim())
            self.colorbar.update_bruteforce(m)
            im.set_cmap(m.get_cmap())
            im.set_clim(m.get_clim())
            #self.canvas.draw_idle()
        #self.colormapper.callbacksSM.connect('changed',on_changed)
        #self.vmin = min(self.vmin, numpy.min(snap))
        #self.vmax = max(self.vmax, numpy.max(snap))
        #self.colormapper.set_clim(vmin=self.vmin,vmax=self.vmax)
        #self.colorbar.set_clim(vmin=self.vmin,vmax=self.vmax)
        #self.autoaxes()
        self.attach_rectangle()
        self.canvas.draw_idle()
        self.onIntegrate()
        self.on_int_draw()
        
    def on_int_draw(self):
        self.onIntegrate()
        ax=self.integration_axes
        ax.clear()
        omega=self.omega
        intensity=self.intensity
        intensity_err=self.intensity_err
        points=ax.errorbar(omega,intensity,yerr=intensity_err,marker='s',mfc='red',mec='red')
        ax.set_ylabel('Intensity (arb. units)')
        ax.set_xlabel(r'$\omega$')
        ax.semilogy()
        self.integration_canvas.draw_idle()
        
    def on_colormap_changed(self,m):
        print "changed",m
        #self.colorbar.set_cmap(m.get_cmap())
        #self.colorbar.set_clim(m.get_clim())
        #self.colorbar.update_bruteforce(m)
        #self.axes[0].set_cmap(m.get_cmap())
        #self.axes[0].set_clim(m.get_clim())
        #self.canvas.draw_idle()
        
    def init_rect(self,rwidth=8, rheight=8,rx=8,ry=8):
        self.rwidth=32/4
        self.rheight=32/4
        self.rx=32/2-rwidth/2
        self.ry=32/2-rheight/2
        self.rvisible=True
        rect=Rectangle((self.rx,self.ry),self.rwidth,self.rheight,fill=True, fc='none',ec='black', visible=self.rvisible)
        ax=self.axes[0]
        ax.add_patch(rect)
        #ax.axis([0,10,0,10])
        self.dr=dragable_rectangle.DraggableRectangle(rect)
        self.dr.myconnect()
        self.dr.connect(self.dr,SIGNAL("rectangleMoved"),self.on_int_draw)
  
    def attach_rectangle(self):
        #rect=Rectangle((self.rx,self.ry),self.rwidth,self.rheight,fill=True, fc='none',ec='black', visible=self.rvisible)
        ax=self.axes[0]
        self.dr.mydisconnect()
        ax.add_patch(self.dr.rect)
        #ax.axis([0,10,0,10])
        #self.dr=dragable_rectangle.DraggableRectangle(rect)
        self.dr.myconnect()
        
        
    def set_rectangle(self):
        self.dr.rect.set_x(self.rx)
        self.dr.rect.set_y(self.ry)
        self.dr.rect.set_width(self.rwidth)
        self.dr.rect.set_height(self.rheight)
        
        
    def on_change_box(self):
        text=str(self.textbox.text())
        xmin,xmax,ymin,ymax=numpy.array(text.split(),'float64')
        self.rx=xmin
        self.ry=ymin
        self.rwidth=math.fabs(xmax-xmin)
        self.rheight=math.fabs(ymax-ymin)
        self.set_rectangle()
        self.on_draw()
        
    def create_frame_tth_window(self):
        self.tth_frame=QWidget()
        self.tth_dock=QDockWidget("Two Theta",self)
        self.tth_dock.setObjectName("TwoTheta")
        self.tth_dock.setAllowedAreas(Qt.LeftDockWidgetArea |Qt.RightDockWidgetArea)      
        self.tth_dock.setWidget(self.tth_frame)
        self.addDockWidget(Qt.LeftDockWidgetArea,self.tth_dock)
        #create a window to show the values integrated over chi
        self.tth_fig=Figure((5.0, 4.0), dpi=self.dpi)
        self.tth_canvas = FigureCanvas(self.tth_fig)
        self.tth_canvas.setParent(self.tth_frame)
        
        #self.tth_canvas.mpl_connect('key_press_event',self.onKeyPressTth)
        self.tth_canvas.mpl_connect('button_press_event',self.onButtonPressTth)
        #self.tth_canvas.mpl_connect('scroll_event',self.onWheelTth)
        
        self.tth_axes=self.tth_fig.add_subplot(111)
        #self.tth_canvas.mpl_connect('pick_event', self.on_pick)
        self.tth_canvas.setFocusPolicy( Qt.ClickFocus )
        self.tth_canvas.setFocus()
        self.tth_toolbar= NavigationToolbar(self.tth_canvas, self.tth_frame)
        
        vbox_tth = QVBoxLayout()
        vbox_tth.addWidget(self.tth_canvas)
        vbox_tth.addWidget(self.tth_toolbar)
        self.tth_frame.setLayout(vbox_tth)


    def create_main_frame(self):
        self.main_frame = QWidget()
        self.integration_frame = QWidget()
        
        self.integration_dock=QDockWidget("integration",self)
        self.integration_dock.setObjectName("integration")
        self.integration_dock.setAllowedAreas(Qt.LeftDockWidgetArea |Qt.RightDockWidgetArea)
        
        self.integration_dock.setWidget(self.integration_frame)
        self.addDockWidget(Qt.RightDockWidgetArea,self.integration_dock)
        #self.integration_dock.setFloating(True)
        #self.integration_dock.setFeatures(QDockWidget.DockWidgetClosable|
        #                                  QDockWidget.DockWidgetMovable|
        #                                  ^QDockWidget.DockWidgetClosable)
        
        
        # Create the mpl Figure and FigCanvas objects. 
        # 5x4 inches, 100 dots-per-inch
        #
        self.dpi = 100
        self.fig = Figure((5.0, 4.0), dpi=self.dpi)
        self.canvas = FigureCanvas(self.fig)
        self.canvas.setParent(self.main_frame)
        
        #create a window to show the values integrated over the displayed rectangle
        self.integration_fig=Figure((5.0, 4.0), dpi=self.dpi)
        self.integration_canvas = FigureCanvas(self.integration_fig)
        self.integration_canvas.setParent(self.integration_frame)
        
        #self.create_frame_tth_window()
        
        
        self.canvas.mpl_connect('key_press_event',self.onKeyPress)
        self.canvas.mpl_connect('button_press_event',self.onButtonPress)
        self.canvas.mpl_connect('scroll_event',self.onWheel)
        
        
        
        # Since we have only one plot, we can use add_axes 
        # instead of add_subplot, but then the subplot
        # configuration tool in the navigation toolbar wouldn't
        # work.
        #
        self.axes = [self.fig.add_subplot(111)]
        
        
        # Create the colorbar
        # Provide an empty handle to attach colormap properties
        self.coloraxes = self.fig.add_axes([0.91, 0.2, 0.04, 0.6])
        self.colormapper = matplotlib.image.FigureImage(self.fig)
        self.colormapper.set_array(numpy.ones(1))
        self.colorbar = self.fig.colorbar(self.colormapper,self.coloraxes)
        self.vmin, self.vmax =numpy.inf,-numpy.inf
        self.colormapper.set_clim(vmin=1,vmax=10)
        self.init_rect()
        self.attach_rectangle()
        
        self.integration_axes=self.integration_fig.add_subplot(111)
        #TODO link the axes together
        
        
        myfilestr=r'192613.dat' 
        myfilestr=r'192601.dat' 
        myfilestr=r'192612.dat' 
        myfilestr=r'192632.dat' 
        myfilestr=r'192623.dat' 
        myfilestr=r'192609.dat' 
       # myfilestr=r'192661.dat' 
        self.data=d10reader.reader2(myfilestr) 
        self.data.frames=numpy.array(self.data.frames)
        self.original_data=copy.deepcopy(self.data)
        self.xscale = 'linear'
        self.yscale = 'linear'
        self.zscale = 'linear'
        self.vmin, self.vmax =numpy.min(self.data.frames)+0.5,numpy.max(self.data.frames)                 #numpy.inf,-numpy.inf
        self.colormapper.set_clim(vmin=self.vmin,vmax=self.vmax)
        self.set_vscale('log')
        self.grid=True       
        
        
        
        # Bind the 'pick' event for clicking on one of the bars
        #
        self.canvas.mpl_connect('pick_event', self.on_pick)
        self.canvas.setFocusPolicy( Qt.ClickFocus )
        self.canvas.setFocus()
        
        self.integration_canvas.setFocusPolicy( Qt.ClickFocus )
        self.integration_canvas.setFocus()
        
        
        
        # Create the navigation toolbar, tied to the canvas
        #
        self.mpl_toolbar = NavigationToolbar(self.canvas, self.main_frame)
        
        self.integration_toolbar= NavigationToolbar(self.integration_canvas, self.integration_frame)
        
        # Other GUI controls
        # 
        self.textbox = QLineEdit()
        self.textbox.setMinimumWidth(200)
        self.connect(self.textbox, SIGNAL('editingFinished ()'), self.on_change_box)
        
        self.draw_button = QPushButton("&Box")
        self.connect(self.draw_button, SIGNAL('clicked()'), self.on_change_box)
        
        self.grid_cb = QCheckBox("Show &Grid")
        self.grid_cb.setChecked(False)
        self.connect(self.grid_cb, SIGNAL('stateChanged(int)'), self.on_draw)
        
        #slider_label = QLabel('Bar width (%):')
        #self.slider = QSlider(Qt.Horizontal)
        #self.slider.setRange(0, len(self.data.frames)-1)
        #self.slider.setValue(len(self.data.frames)/2)
        #self.slider.setTracking(True)
        #self.slider.setTickPosition(QSlider.TicksBothSides)
        #self.connect(self.slider, SIGNAL('valueChanged(int)'), self.on_draw)
        
        #Swap out slider for spinbox
        
        spinbox_label=QLabel('Frame (#):')
        self.spinbox=QSpinBox()
        self.spinbox.setRange(0, len(self.data.frames)-1)
        self.spinbox.setValue(len(self.data.frames)/2)
        #self.spinbox.setTracking(True)
        #self.spinbox.setTickPosition(QSlider.TicksBothSides)
        self.connect(self.spinbox, SIGNAL('valueChanged(int)'), self.on_draw)
        #
        # Layout with box sizers
        # 
        hbox = QHBoxLayout()
        
        for w in [  self.textbox, self.draw_button, self.grid_cb,
                    spinbox_label, self.spinbox]:
            hbox.addWidget(w)
            hbox.setAlignment(w, Qt.AlignVCenter)
        
        vbox = QVBoxLayout()
        vbox.addWidget(self.canvas)
        vbox.addWidget(self.mpl_toolbar)
        vbox.addLayout(hbox)
        
        vbox_int = QVBoxLayout()
        vbox_int.addWidget(self.integration_canvas)
        vbox_int.addWidget(self.integration_toolbar)
        
        
        
        self.main_frame.setLayout(vbox)
        self.integration_frame.setLayout(vbox_int)
        self.create_frame_tth_window()
        self.setCentralWidget(self.main_frame)
        #self.integration_frame.activateWindow()
        #self.integration_dock.DockWidgetClosable=False
        self.integration_dock.setFeatures(QDockWidget.DockWidgetMovable|QDockWidget.DockWidgetFloatable)
        self.tth_dock.setFeatures(QDockWidget.DockWidgetMovable|QDockWidget.DockWidgetFloatable)

        #self.setSizeGripEnabled(True)
        
        self.form = MainForm(self)
        self.form.connect(self.form,SIGNAL("fileClicked"),self.on_reload_frame)
        self.form.resize(850, 620)
        self.form.show()
        self.form.raise_()
        self.form.activateWindow()

    def on_reload_frame(self,myfilestr):
        print 'myfilestr', myfilestr
        self.data=d10reader.reader2(myfilestr) 
        self.data.frames=numpy.array(self.data.frames)
        self.original_data=copy.deepcopy(self.data)
        self.xscale = 'linear'
        self.yscale = 'linear'
        self.zscale = 'linear'
        self.vmin, self.vmax =numpy.min(self.data.frames)+0.5,numpy.max(self.data.frames)                 #numpy.inf,-numpy.inf
        self.colormapper.set_clim(vmin=self.vmin,vmax=self.vmax)
        self.set_vscale('log')
        self.on_draw()
        self.on_draw_tth()
        self.onIntegrate()
        self.on_int_draw()
        
        
    def create_status_bar(self):
        self.status_text = QLabel("This is a demo")
        self.statusBar().addWidget(self.status_text, 1)
        
    def create_menu(self):        
        self.file_menu = self.menuBar().addMenu("&File")
        
        load_file_action = self.create_action("&Save plot",
            shortcut="Ctrl+S", slot=self.save_plot, 
            tip="Save the plot")
        quit_action = self.create_action("&Quit", slot=self.close, 
            shortcut="Ctrl+Q", tip="Close the application")
        
        self.add_actions(self.file_menu, 
            (load_file_action, None, quit_action))
        
        self.help_menu = self.menuBar().addMenu("&Help")
        about_action = self.create_action("&About", 
            shortcut='F1', slot=self.on_about, 
            tip='About the demo')
        
        self.add_actions(self.help_menu, (about_action,))

    def add_actions(self, target, actions):
        for action in actions:
            if action is None:
                target.addSeparator()
            else:
                target.addAction(action)

    def create_action(  self, text, slot=None, shortcut=None, 
                        icon=None, tip=None, checkable=False, 
                        signal="triggered()"):
        action = QAction(text, self)
        if icon is not None:
            action.setIcon(QIcon(":/%s.png" % icon))
        if shortcut is not None:
            action.setShortcut(shortcut)
        if tip is not None:
            action.setToolTip(tip)
            action.setStatusTip(tip)
        if slot is not None:
            self.connect(action, SIGNAL(signal), slot)
        if checkable:
            action.setCheckable(True)
        return action
    
    def onKeyPressTth(self,event):
        #if not event.inaxes: return
        # Let me zoom and unzoom even without the scroll wheel
        if event.key == 'a':
            setattr(event,'step',5.)
            self.onWheel(event)
        elif event.key == 'z':
            setattr(event,'step',-5.)
            self.onWheel(event)
        elif event.key == 'l':
            if self.vscale=='linear':
                self.set_vscale('log')
                self.on_draw()
                self.on_draw_tth()
            else:
                self.set_vscale('linear')
                self.on_draw()
                self.on_draw_tth()
    
    def onKeyPress(self,event):
        #if not event.inaxes: return
        # Let me zoom and unzoom even without the scroll wheel
        if event.key == 'a':
            setattr(event,'step',5.)
            self.onWheel(event)
        elif event.key == 'z':
            setattr(event,'step',-5.)
            self.onWheel(event)
        elif event.key == 'l':
            if self.vscale=='linear':
                self.set_vscale('log')
                self.on_draw()
                self.on_draw_tth()
            else:
                self.set_vscale('linear')
                self.on_draw()
                self.on_draw_tth()
        elif event.key == 'r':
            print "reseting vmin,vmax"
            #self.vmin=self.data.frames[self.spinbox.value()].min()
            #self.vmax=self.data.frames[self.spinbox.value()].max()
            self.colorbar.set_clim(vmin=self.vmin,vmax=self.vmax)
            #self.set_vscale(self.vscale)
            #self.colorbar.set_clim(vmin=self.vmin,vmax=self.vmax)
            self.on_draw()
            
        elif event.key == 'b':
            self.onIntegrate()
            self.on_int_draw()
    
    def onIntegrate(self):
        print 'integrating'
        pvertices=self.dr.rect.get_verts()
        vertices=self.axes[0].transData.inverted().transform(pvertices)
        print vertices
        xmin=vertices[:,0].min()
        xmax=vertices[:,0].max()
        ymin=vertices[:,1].min()
        ymax=vertices[:,1].max()
        print xmin,xmax,ymin,ymax
        intensity=[]
        intensity_err=[]
        i=0
        for frame in self.data.frames:
            val=0
            
            for x in range(xmin,xmax):
                for y in range(ymin,ymax):
                    val=val+frame[x,y]
                    
            #if self.data.angle2==[]:
                #tth=self.data.analyzer_tth
            #else:
                #tth=self.data.angle2[i]
            #lfactor=math.sin(math.radians(tth))
            #lfactor=1.0
            intensity.append(val)
            #intensity_err.append(math.sqrt(val))
            i=i+1
        tth=self.data.tth
        lfactor=math.sin(math.radians(tth))
        print 'lfactor',lfactor,'tth',tth
        lfactor=1.0
        intensity_err=numpy.sqrt(numpy.array(intensity))*lfactor
        intensity=numpy.array(intensity)*lfactor
        omega=numpy.array(self.data.angle1)
        outdata=numpy.array([omega,intensity,intensity_err])
        numpy.savetxt('myfile.txt', outdata.T, fmt="%12.6G %12.6G %12.6G")
        self.omega=omega
        self.intensity=intensity
        self.intensity_err=intensity_err
        
        

    def onButtonPress(self,event):
        # TODO: finish binder so that it allows tagging of fully qualified
        # events on an artist by artist basis.
        if event.inaxes == self.coloraxes:
            #self.colormapper.set_clim(vmin=self.vmin,vmax=self.vmax)
            self.canvas.draw_idle()
            #self.tth_canvas.draw_idle()
        elif event.inaxes != None:
            #self.autoaxes()
            self.canvas.draw_idle()
            #self.tth_canvas.draw_idle()
            
    def onButtonPressTth(self,event):
        # TODO: finish binder so that it allows tagging of fully qualified
        # events on an artist by artist basis.
        if event.inaxes == self.coloraxes:
            #self.colormapper.set_clim(vmin=self.vmin,vmax=self.vmax)
            #self.canvas.draw_idle()
            #self.tth_canvas.draw_idle()
            pass
        elif event.inaxes != None:
            frame=event.xdata
            self.spinbox.setValue(int(frame))
            self.on_draw()
            #self.autoaxes()
            #self.canvas.draw_idle()
            #self.tth_canvas.draw_idle()

    def onWheel(self, event):
        """
        Process mouse wheel as zoom events
        """
        ax = event.inaxes
        step = event.step

        # Icky hardcoding of colorbar zoom.
        if ax == self.coloraxes:
            # rescale colormap: the axes are already scaled to 0..1,
            # so use bal instead of pt for centering
            lo,hi = self.colormapper.get_clim()
            lo,hi = _rescale(lo,hi,step,bal=event.ydata,scale=self.vscale)
            print 'in axes',lo,hi,self.vmin,self.vmax
            self.colormapper.set_clim(lo,hi)
            self.colorbar.set_cmap(self.colormapper.get_cmap())
            self.colorbar.set_clim(self.colormapper.get_clim())
            self.colorbar.update_bruteforce(self.colormapper)
            self.im.set_cmap(self.colormapper.get_cmap())
            self.im.set_clim(self.colormapper.get_clim())
            print 'out axes'
            #self.im_tth.set_cmap(self.colormapper.get_cmap())
            #self.im_tth.set_clim(self.colormapper.get_clim())
        elif ax != None:
            # Event occurred inside a plotting area
            lo,hi = ax.get_xlim()
            lo,hi = _rescale(lo,hi,step,pt=event.xdata)
            ax.set_xlim((lo,hi))

            lo,hi = ax.get_ylim()
            lo,hi = _rescale(lo,hi,step,pt=event.ydata)
            ax.set_ylim((lo,hi))
        else:
            # Check if zoom happens in the axes
            xdata,ydata = None,None
            x,y = event.x,event.y
            for ax in self.axes:
                insidex,_ = ax.xaxis.contains(event)
                if insidex:
                    xdata,_ = ax.transAxes.inverse_xy_tup((x,y))
                    #print "xaxis",x,"->",xdata
                insidey,_ = ax.yaxis.contains(event)
                if insidey:
                    _,ydata = ax.transAxes.inverse_xy_tup((x,y))
                    #print "yaxis",y,"->",ydata
            if xdata is not None:
                lo,hi = ax.get_xlim()
                lo,hi = _rescale(lo,hi,step,bal=xdata)
                ax.set_xlim((lo,hi))
            if ydata is not None:
                lo,hi = ax.get_ylim()
                lo,hi = _rescale(lo,hi,step,bal=ydata)
                ax.set_ylim((lo,hi))
            
        self.canvas.draw_idle()
        #self.tth_canvas.draw_idle()

        
    def autoaxes(self):
        return
        #bbox = bbox_union([ax.dataLim for ax in self.axes])
        #xlims = bbox.intervalx
        #ylims = bbox.intervaly
        #self.pp.axis(list(xlims)+list(ylims))

    # These are properties which the user should control but for which
    # a particular plottable might want to set a reasonable default.
    # For now leave control with the plottable.
    def set_xscale(self, scale='linear'):
        for axes in self.axes: axes.set_xscale(scale)
        self.xscale = scale

    def get_xscale(self):
        return self.xscale

    def set_yscale(self, scale='linear'):
        for axes in self.axes: axes.set_yscale(scale)
        self.yscale = scale

    def get_yscale(self):
        return self.yscale

    def set_vscale(self, scale='linear'):
        """Alternate between log and linear colormap"""
        if scale == 'linear':
            vmapper = Normalize(*self.colormapper.get_clim())
            vformat = matplotlib.ticker.ScalarFormatter()
            vlocate = matplotlib.ticker.AutoLocator()
        else:
            vmapper = LogNorm(*self.colormapper.get_clim())
            vformat = matplotlib.ticker.LogFormatterMathtext(base=10,labelOnlyBase=False)
            vlocate = matplotlib.ticker.LogLocator(base=10)
        self.colormapper.set_norm(vmapper)
        self.colorbar.formatter = vformat
        self.colorbar.locator = vlocate

        self.vscale = scale

    def get_vscale(self):
        return self.vscale

    ## The following is plottable functionality
    def properties(self,prop):
        """Set some properties of the graph.

        The set of properties is not yet determined.
        """
        # The particulars of how they are stored and manipulated (e.g., do
        # we want an inventory internally) is not settled.  I've used a
        # property dictionary for now.
        #
        # How these properties interact with a user defined style file is
        # even less clear.

        # Properties defined by plot
        self.xbox.set_text(r"$%s$" % prop["xlabel"])
        self.ybox.set_text(r"$%s$" % prop["ylabel"])
        self.vbox.set_text(r"$%s$" % prop["vlabel"])
        self.tbox.set_text(r"$%s$" % prop["title"])

        # Properties defined by user
        #self.axes.grid(True)

    def xaxis(self,label,units):
        """xaxis label and units.

        Axis labels know about units.

        We need to do this so that we can detect when axes are not
        commesurate.  Currently this is ignored other than for formatting
        purposes.
        """
        if units != "": label = label + " (" + units + ")"
        self.xbox.set_text(r"$%s$" % (label))
        pass

    def yaxis(self,label,units):
        """yaxis label and units."""
        if units != "": label = label + " (" + units + ")"
        self.ybox.set_text(r"$%s$" % (label))
        pass

    def vaxis(self,label,units):
        """vaxis label and units."""
        if units != "": label = label + " (" + units + ")"
        self.vbox.set_text(r"$%s$" % (label))
        pass



def main():
    app = QApplication(sys.argv)
    form = AppForm()
    form.show()
    app.exec_()


if __name__ == "__main__":
    main()
