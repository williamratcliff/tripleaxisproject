import wx
import wxaddons.sized_controls as sc
import  wx.lib.intctrl
import  wx.grid as  gridlib
import numpy as N
import sys,os

class MyApp(wx.App):
    def __init__(self, redirect=False, filename=None, useBestVisual=False, clearSigInt=True):
        wx.App.__init__(self,redirect,filename,clearSigInt)


    def OnInit(self):
        return True
    
    
class mFormValidator(wx.PyValidator):
    def __init__(self,data,key):
        wx.PyValidator.__init__(self)
        self.data=data
        self.key=key
        #print 'FormValidator init', key
        #self.TransferToWindow()

    def Clone(self):
        return mFormValidator(self.data,self.key)

    def Validate(self,win):
        return True

    def TransferToWindow(self):
        #print 'Form TransferToWindow',self.key
        
        checkctrl=self.GetWindow()
        ##print 'checkctrl',checkctrl
        ##print self.__dict__
        #print 'key',self.key
        #print 'dict', self.data._array[self.key]
        #print 'data',self.data.__dict__[self.key]
        checkctrl.SetValue(self.data._array[self.key])
        return True

    def TransferFromWindow(self):
        #print 'TransferFromWindow'
        checkctrl=self.GetWindow()
        self.data._array[self.key]=int(checkctrl.GetValue())
        return True    
    
class FormDialog(sc.SizedDialog):
    def __init__(self, parent, id):
        sc.SizedDialog.__init__(self, None, -1, "Calculate Dispersion",
                        style=wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER)#| wx.WS_EX_VALIDATE_RECURSIVELY)

        pane = self.GetContentsPane()
        pane.SetSizerType("vertical")
        FilePane = sc.SizedPanel(pane, -1)
        FilePane.SetSizerType("vertical")
        FilePane.SetSizerProps(expand=True)
        self.spinfile=''
        self.spinfilectrl=wx.StaticText(FilePane, -1, "Spin File:%s"%(self.spinfile,))
        b = wx.Button(FilePane, -1, "Browse", (50,50))
        self.Bind(wx.EVT_BUTTON, self.OnOpen, b)


        InteractionFilePane = sc.SizedPanel(pane, -1)
        InteractionFilePane.SetSizerType("vertical")
        InteractionFilePane.SetSizerProps(expand=True)
        self.interactionfile=''
        self.interactionfilectrl=wx.StaticText(FilePane, -1, "Interaction File:%s"%(self.interactionfile,))
        b = wx.Button(FilePane, -1, "Browse", (50,50))
        self.Bind(wx.EVT_BUTTON, self.OnOpenInt, b)




        self.kx=1.0
        self.ky=0.0
        self.kz=0.0
        DirectionPane = sc.SizedPanel(pane, -1)
        DirectionPane.SetSizerType("vertical")
        DirectionPane.SetSizerProps(expand=True)
        wx.StaticText(DirectionPane, -1, "Scan Direction")
        
        DirectionsubPane = sc.SizedPanel(pane, -1)
        DirectionsubPane.SetSizerType("horizontal")
        DirectionsubPane.SetSizerProps(expand=True)
        wx.StaticText(DirectionsubPane, -1, "qx")
        qx=wx.TextCtrl(DirectionsubPane, -1, str(self.kx),validator=mFormValidator(self.kx,self.kx))
        self.Bind(wx.EVT_TEXT, self.Evtqx, qx)
        wx.StaticText(DirectionsubPane, -1, "qy")
        qy=wx.TextCtrl(DirectionsubPane, -1, str(self.ky),validator=mFormValidator(self.ky,self.ky))
        self.Bind(wx.EVT_TEXT, self.Evtqx, qy)
        wx.StaticText(DirectionsubPane, -1, "qz")
        qz=wx.TextCtrl(DirectionsubPane, -1, str(self.kz),validator=mFormValidator(self.ky,self.ky))
        self.Bind(wx.EVT_TEXT, self.Evtqx, qz)
        print 'Directed'



        # add dialog buttons
        self.SetButtonSizer(self.CreateStdDialogButtonSizer(wx.OK | wx.CANCEL))
        self.TransferDataToWindow()
        # a little trick to make sure that you can't resize the dialog to
        # less screen space than the controls need
        self.Fit()
        self.SetMinSize(self.GetSize())


    def Evtqx(self,evt):
        self.kx=float(wx.wxTextCtrl.GetValue(self))
 


    def OnOpen(self,event):
        # Create the dialog. In this case the current directory is forced as the starting
        # directory for the dialog, and no default file name is forced. This can easilly
        # be changed in your program. This is an 'open' dialog, and allows multitple
        # file selections as well.
        #
        # Finally, if the directory is changed in the process of getting files, this
        # dialog is set up to change the current working directory to the path chosen.

        #defaultDir=os.getcwd()
        #defaultDir=r'C:\polcorrecter\data'
        defaultDir=wx.ConfigBase.Get().GetPath()
        wildcard="files (*.txt)|*.txt|All files (*.*)|*.*"
        dlg = wx.FileDialog(
            self, message="Choose a spin configuration file",
            defaultDir=defaultDir,
            defaultFile="",
            wildcard=wildcard,
            style=wx.OPEN | wx.CHANGE_DIR
            )

        # Show the dialog and retrieve the user response. If it is the OK response,
        # process the data.
        if dlg.ShowModal() == wx.ID_OK:
            # This returns a Python list of files that were selected.
            paths = dlg.GetPaths()
            #self.log.WriteText('You selected %d files:' % len(paths))
            self.spinfile=paths[0].encode('ascii')
            self.spinfilectrl.SetLabel("Spin File:%s"%(self.spinfile,))
            #wx.StaticText(FilePane, -1, "CellFile:%s"%(self.groupdata['cellfile'],))
        # Destroy the dialog. Don't do this until you are done with it!
        # BAD things can happen otherwise!
        dlg.Destroy()


    def OnOpenInt(self,event):
        # Create the dialog. In this case the current directory is forced as the starting
        # directory for the dialog, and no default file name is forced. This can easilly
        # be changed in your program. This is an 'open' dialog, and allows multitple
        # file selections as well.
        #
        # Finally, if the directory is changed in the process of getting files, this
        # dialog is set up to change the current working directory to the path chosen.

        #defaultDir=os.getcwd()
        #defaultDir=r'C:\polcorrecter\data'
        defaultDir=wx.ConfigBase.Get().GetPath()
        wildcard="files (*.txt)|*.txt|All files (*.*)|*.*"
        dlg = wx.FileDialog(
            self, message="Choose an Interaction file",
            defaultDir=defaultDir,
            defaultFile="",
            wildcard=wildcard,
            style=wx.OPEN | wx.CHANGE_DIR
            )

        # Show the dialog and retrieve the user response. If it is the OK response,
        # process the data.
        if dlg.ShowModal() == wx.ID_OK:
            # This returns a Python list of files that were selected.
            paths = dlg.GetPaths()
            #self.log.WriteText('You selected %d files:' % len(paths))
            self.interactionfile=paths[0].encode('ascii')
            self.interactionfilectrl.SetLabel("Interaction File:%s"%(self.interactionfile,))
            #wx.StaticText(FilePane, -1, "CellFile:%s"%(self.groupdata['cellfile'],))
        # Destroy the dialog. Don't do this until you are done with it!
        # BAD things can happen otherwise!
        dlg.Destroy()



if __name__=='__main__':
    app=MyApp()
    dlg=FormDialog(parent=None,id=-1)
    result=dlg.ShowModal()
    if result==wx.ID_OK:
        print "OK"
    else:
        print "Cancel"
    dlg.Destroy()
