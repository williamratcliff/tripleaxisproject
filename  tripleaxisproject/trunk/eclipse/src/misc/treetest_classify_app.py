import treetest_classify
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
import matplotlib.colors
from matplotlib.font_manager import FontProperties
from matplotlib.colors import Normalize,LogNorm
import ticker



class AppForm(QMainWindow):
    def __init__(self, parent=None):
        QMainWindow.__init__(self, parent)
        self.setWindowTitle('Demo: PyQt with matplotlib')

        self.create_menu()
        self.create_main_frame()
        self.create_status_bar()
        #self.on_draw()

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
        

    def onClear(self):
        ax=self.axes[0]
        ax.clear()
        return
        

    def on_draw(self):
        """ Redraws the figure
        """
        #str = unicode(self.textbox.text())
        #self.data = map(int, str.split())
        # clear the axes and redraw the plot anew
        #
        
        
        
        if self.x==None or self.y==None or self.yerr==None:
            return
        ax=self.axes[0]
        #ax.clear()
        x=self.x
        y=self.y
        yerr=self.yerr
        xlabel=self.xlabel
        ylabel=self.ylabel
        
        #for curr_ax in self.axes:
        #    curr_ax.grid(alpha=0.8, visible=self.grid_cb.isChecked())
        self.im=ax.errorbar(x,y,yerr=yerr,fmt='rs')
        #ax.set_xlabel(r'$\omega$')
        #ax.set_ylabel('Intensity (arb. units)')
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        #ax.semilogy()
        self.canvas.draw_idle()
        #self.on_int_draw()  
        



    def create_main_frame(self):
        self.main_frame = QWidget()
        #self.integration_frame = QWidget()
        
        self.tree=treetest_classify.myTreeView()
        self.tree.connect(self.tree,SIGNAL("plot"),self.onPlot)
        self.tree.connect(self.tree,SIGNAL("clearplot"),self.onClear)
        self.tree_dock=QDockWidget("integration",self)
        self.tree_dock.setObjectName("Files")
        self.tree_dock.setAllowedAreas(Qt.LeftDockWidgetArea |Qt.RightDockWidgetArea)
        
        self.tree_dock.setWidget(self.tree)
        self.addDockWidget(Qt.RightDockWidgetArea,self.tree_dock)
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
        self.xscale = 'linear'
        self.yscale = 'linear'
        self.zscale = 'linear'
        self.grid=True       
        
        
        
        # Bind the 'pick' event for clicking on one of the bars
        #
        self.canvas.mpl_connect('pick_event', self.on_pick)
        self.canvas.setFocusPolicy( Qt.ClickFocus )
        self.canvas.setFocus()
      
        # Create the navigation toolbar, tied to the canvas
        #
        self.mpl_toolbar = NavigationToolbar(self.canvas, self.main_frame)
     
        vbox = QVBoxLayout()
        vbox.addWidget(self.canvas)
        vbox.addWidget(self.mpl_toolbar)
      
        self.main_frame.setLayout(vbox)
        self.setCentralWidget(self.main_frame)
        self.tree_dock.setFeatures(QDockWidget.DockWidgetMovable|QDockWidget.DockWidgetFloatable)

        #self.setSizeGripEnabled(True)
        
        #self.form = MainForm(self)
        #self.form.connect(self.form,SIGNAL("fileClicked"),self.on_reload_frame)
        #self.form.resize(850, 620)
        #self.form.show()
        #self.form.raise_()
        #self.form.activateWindow()
        
    def onButtonPress(self,event):
        # TODO: finish binder so that it allows tagging of fully qualified
        # events on an artist by artist basis.
 
        if event.inaxes != None:
            #self.autoaxes()
            self.canvas.draw_idle()
            #self.tth_canvas.draw_idle()


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
    
    def onPlot(self,plotdict):
        scantype=plotdict['scantype']
        #if scantype=='th':
        print 'onPlot'
        self.x=plotdict['data']['x']
        self.xlabel=plotdict['xlabel']
        self.y=plotdict['data']['y']
        self.yerr=plotdict['data']['yerr']
        self.ylabel=plotdict['ylabel']
        self.on_draw()
    
    def onKeyPress(self,event):
        #if not event.inaxes: return
        # Let me zoom and unzoom even without the scroll wheel
        if event.key == 'a':
            setattr(event,'step',5.)
            self.onWheel(event)
        elif event.key == 'z':
            setattr(event,'step',-5.)
            self.onWheel(event)

    

        #intensity_err=numpy.sqrt(numpy.array(intensity))*lfactor
        #intensity=numpy.array(intensity)*lfactor
        #omega=numpy.array(self.data.angle1)
        #outdata=numpy.array([omega,intensity,intensity_err])
        #numpy.savetxt('myfile.txt', outdata.T, fmt="%12.6G %12.6G %12.6G")

        


    def onWheel(self, event):
        """
        Process mouse wheel as zoom events
        """
        ax = event.inaxes
        step = event.step


        if ax != None:
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







def main():
    app = QApplication(sys.argv)
    form = AppForm()
    form.show()
    app.exec_()


if __name__ == "__main__":
    main()
