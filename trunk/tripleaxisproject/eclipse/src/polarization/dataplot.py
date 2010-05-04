import wx.lib.newevent

#import PropertyDialog
(FunctionParamEvent, EVT_FUNC_PARS) = wx.lib.newevent.NewEvent()
(FunctionFitEvent, EVT_FUNC_FIT) = wx.lib.newevent.NewEvent()
(FitDoneEvent, EVT_FIT_DONE) = wx.lib.newevent.NewEvent()
#(FunctionScaleEvent,EVT_PROPERTY) = wx.lib.newevent.NewEvent()
#import PropertyDialog
from PlotPanel import PlotPanel
import numpy as N
import pylab
from plottables import Plottable, Graph, Data1D
import os
import readncnr2 as readncnr
#import math

class MyApp(wx.App):
    def __init__(self, redirect=False, filename=None, useBestVisual=False, clearSigInt=True):
        wx.App.__init__(self,redirect,filename,clearSigInt)


    def OnInit(self):
        return True



class PolarizedPanel1D(PlotPanel):
    def __init__(self, parent, id = -1, color = None,\
        dpi = None, style = wx.NO_FULL_REPAINT_ON_RESIZE, **kwargs):
        PlotPanel.__init__(self, parent, id = id, style = style, **kwargs)

        self.parent = parent

        self.figure.subplots_adjust(bottom=.25)

        self.set_yscale('log')
        self.set_xscale('linear')
        xtrans="x"
        ytrans="Log(y)"
        self.setTrans(xtrans,ytrans)
##        self.x = pylab.arange(self.qmin, self.qmax+self.qstep*0.01, self.qstep)
##        # Error on x
##        self.dx = numpy.zeros(len(self.x))
##        # Intensity values
##        y  = numpy.ones(len(self.x))
##        print y
##        # Error on y
##        self.dy = numpy.zeros(len(self.x))

        # Plottables
##        self.file_data = Data1D(x=[], y=[], dx=None, dy=None)
##        self.file_data.name = "Loaded 1D data"
##
##        self.file_data1 = Data1D(x=[], y=[], dx=None, dy=None)
##        self.file_data1.name = "y= exp(A + bx**2)"

        # Graph
        self.graph = Graph()
        self.graph.xaxis('\\rm{q} ', 'A^{-1}')
        self.graph.yaxis("\\rm{Intensity} ","cm^{-1}")
##        self.graph.add(self.file_data)
##        #self.graph.add(self.file_data1)
##        self.graph.render(self)



    def onPlot(self):
        self.file_data.name = "Loaded 1D data"
        self.graph.xaxis('\\rm{q} ', 'A^{-1}')
        self.graph.yaxis("\\rm{Intensity} ","cm^{-1}")

        self.graph.add(self.file_data)
        self.graph.render(self)
        self.subplot.figure.canvas.draw_idle()


    def returnPlottable(self):
        self.file_data1 = Data1D(x=[], y=[], dx=[], dy=[])
        self.file_data1.name = "y= exp(A + bx**2)"
        return self.file_data1
    def onContextMenu(self, event):
        """
            Pop up a context menu
        """
        # Slicer plot popup menu
        #slicerpop = wx.Menu()
        #slicerpop.Append(314, "&Save 1D model points (%s)" % self.file_data.name,
        #                      'Save randomly oriented data currently displayed')

        #slicerpop.Append(316, '&Load 1D data file')

        #slicerpop.Append(315, '&Toggle Linear/Log intensity scale Y-axis')
        PlotPanel.onContextMenu(self,event)
        #pos = event.GetPosition()
        #pos = self.ScreenToClient(pos)
        #self.PopupMenu(slicerpop, pos)








    def _onEVT_FUNC_PARS(self, event):

        """
            Plot exp(cstA + x^(2) * cstB)
        """
        temp=[]
        fittings.Parameter(self.model, 'A', event.cstA)
        fittings.Parameter(self.model, 'B', event.cstB)
        if self.file_data.x:
            for x_i in self.file_data.x:
                temp.append(self.model.run(x_i))
            self.file_data1.y =temp
            self.file_data1.x= self.file_data.x
        else:
            # xtemp has a default value in case the user doesn't load data
            xtemp = [1, 2, 3, 4, 5, 6]
            for x_i in xtemp:
                temp.append(self.model.run(x_i))
            self.file_data1.x =xtemp
            self.file_data1.y =temp
        self.file_data1.reset_view()
        self.graph.add(self.file_data1)
        self.graph.render(self)
        self.subplot.figure.canvas.draw_idle()

    def _onLoad1DData(self, event):
        """
            Load a data file
        """
        path = None
        dlg = wx.FileDialog(self, "Choose a file", os.getcwd(), "", "*.txt", wx.OPEN)
        if dlg.ShowModal() == wx.ID_OK:
            path = dlg.GetPath()
            mypath = os.path.basename(path)
            print mypath
        dlg.Destroy()

        file_x = []
        file_y = []
        file_dy = []
        file_dx = []
        if not path == None:
            self.path =path
            input_f =  open(path,'r')
            buff = input_f.read()
            lines = buff.split('\n')
            for line in lines:
                try:
                    toks = line.split()
                    x = float(toks[0])
                    y = float(toks[1])
                    #dx = math.sqrt(x)
                    dx=1/x
                    if dx >= x:
                        dx = 0.9*x
                    #dy = math.sqrt(y)
                    dy=1/y
                    if dy >= y:
                        dy = 0.9*y
                    file_x.append(x)
                    file_y.append(y)
                    file_dy.append(dy)
                    file_dx.append(dx)

                except:
                    print "READ ERROR", line

        # Sanity check
        if not len(file_x) == len(file_dx):
            raise ValueError, "X and dX have different length"
        if not len(file_y) == len(file_dy):
            raise ValueError, "y and dy have different length"
        # reset the graph before loading
        self.graph.reset()
        self.file_data.x = file_x
        self.file_data.y = file_y
        self.file_data.dy = file_dy
        #self.file_data.dy = None

        #self.file_data.dx = file_dx
        self.file_data.dx = None

        self.file_data.reset_view()

        self.file_data.name = "Loaded 1D data"
        self.graph.xaxis('\\rm{q} ', 'A^{-1}')
        self.graph.yaxis("\\rm{Intensity} ","cm^{-1}")

        # Set the scale
        self.set_yscale('log')
        self.set_xscale('linear')
        #Add the default transformation of x and y into Property Dialog
        if self.get_xscale()=='log':
            xtrans="Log(x)"
        if self.get_xscale()=='linear':
            xtrans="x"
        if self.get_yscale()=='log':
            ytrans="Log(y)"
        if self.get_yscale()=='linear':
            ytrans="y"
        self.setTrans(xtrans,ytrans)

        #Plot the data
        self.graph.add(self.file_data)
        self. _onEVT_FUNC_PROPERTY()

        #self.graph.render(self)
        #self.subplot.figure.canvas.draw_idle()


    def _onSquaredQ(self, event):
        """
            This method plots Q**2
        """
        self.graph.xaxis('\\rm{q}^2 ', 'A^{-2}')
        self.set_xscale('squared')

        self.graph.render(self)
        self.subplot.figure.canvas.draw_idle()


    def _onLinearQ(self, event):
        """
            This method plots Q
        """
        self.graph.xaxis('\\rm{q} ', 'A^{-1}')
        self.set_xscale('linear')
        self.graph.render(self)
        self.subplot.figure.canvas.draw_idle()

    def _onToggleScale(self, event):
        """
            Toggle the scale of the y-axis
        """
        if self.get_yscale() == 'log':
            self.set_yscale('linear')
        else:
            self.set_yscale('log')
        self.subplot.figure.canvas.draw_idle()

    def _onLogQ(self, event):
        """
            Plot log(q)
        """
        self.set_xscale('log')
        self.graph.xaxis('\\rm{log(q)} ', 'A^{-1}')

        self.graph.render(self)
        self.subplot.figure.canvas.draw_idle()


class TestFrame(wx.Frame):
    def __init__(self,parent,id):
        wx.Frame.__init__(self,parent,id,'Plot Panel',size=(640,200),style=wx.DEFAULT_FRAME_STYLE^wx.CLOSE_BOX)
        self.Bind(wx.EVT_CLOSE,self.OnCloseWindow)
        p = PolarizedPanel1D(self, -1, style=0)
        bs = wx.BoxSizer(wx.VERTICAL)
        bs.Add(p, 1, wx.GROW|wx.ALL|wx.EXPAND, 5)

        self.p=p
        #self.load_data()
        self.SetSizer(bs)

    def OnCloseWindow(self,event):
        self.Destroy()

    def load_data(self):
        myfilestr=r'c:\bifeo3xtal\jan8_2008\9175\fieldscanminusplusreset53631.bt7'
        mydatareader=readncnr.datareader()
        mydata=mydatareader.readbuffer(myfilestr)
        qx1=N.array(mydata.data['qx'])
        I1=N.array(mydata.data['detector'])
        dI1=N.sqrt(I1)
        data={}
        data['x']=qx1
        data['y']=I1
        data['dy']=dI1
        self.data1=data

        myfilestr=r'c:\bifeo3xtal\jan8_2008\9175\fieldscansplusminusreset53630.bt7'
        mydatareader2=readncnr.datareader()
        mydata2=mydatareader2.readbuffer(myfilestr)
        qx2=N.array(mydata2.data['qx'])
        I2=N.array(mydata2.data['detector'])
        dI2=N.sqrt(I2)
        data={}
        data['x']=qx2
        data['y']=I2
        data['dy']=dI2
        self.data2=data
        self.file_data = Data1D(x=data['x'], y=data['y'], dx=None, dy=data['dy'])
        self.file_data.name = "MyData2"
        print 'data added'
        self.p.graph.add(self.file_data)
        print 'graph added'
        self.p.graph.render(self.p)
        print 'rendered'
        self.p.subplot.figure.canvas.draw_idle()
        print 'drawing'



if __name__=='__main__':
    app=MyApp()
    frame=TestFrame(parent=None,id=-1)
    frame.Show()
    app.MainLoop()

