import wx
import wxaddons.sized_controls as sc
import  wx.lib.intctrl
import  wx.grid as  gridlib
import numpy as N
import sys,os

   
    
class mFormValidator(wx.PyValidator):
    def __init__(self,data,key):
        wx.PyValidator.__init__(self)
        print 'FormValidator init', key
        self.data=data
        self.key=key
        #self.TransferToWindow()
    def Clone(self):
        return mFormValidator(self.data,self.key)

    def Validate(self,win):
        print 'validating'
        textCtrl=self.GetWindow()
        text=textCtrl.GetValue()
        
        try:
            value=float(text)
            textCtrl.SetBackgroundColour(wx.SystemSettings_GetColour(wx.SYS_COLOUR_WINDOW))
            textCtrl.Refresh()
            return True
            
        except ValueError:
            wx.MessageBox("This field must be a number","error")
            textCtrl.SetBackgroundColour("pink")
            textCtrl.SetFocus()
            textCtrl.Refresh()
            return False
        

    def TransferToWindow(self):
        print 'Form TransferToWindow',self.key
        
        textCtrl=self.GetWindow()
        ##print 'checkctrl',checkctrl
        print self.__dict__

        textCtrl.SetValue(self.data.get(self.key,""))
        return True

    def TransferFromWindow(self):
        print 'TransferFromWindow'
        textCtrl=self.GetWindow()
        self.data[self.key]=textCtrl.GetValue()
        #self.qfloat=float(textCtrl.GetValue())
        return True    


def WalkTree(parent):
        print 'walking', parent
        parent.SetExtraStyle(wx.WS_EX_VALIDATE_RECURSIVELY)
        for child in parent.GetChildren():
            if child==None:
                print 'child',child
            else:
                WalkTree(child)
                
class FormDialog(sc.SizedDialog):
    def __init__(self, parent, id):
        valstyle=wx.WS_EX_VALIDATE_RECURSIVELY
        sc.SizedDialog.__init__(self, None, -1, "Calculate Dispersion",
                        style=wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER)#| wx.WS_EX_VALIDATE_RECURSIVELY)
        self.SetExtraStyle(wx.WS_EX_VALIDATE_RECURSIVELY)
        pane = self.GetContentsPane()
        pane.SetSizerType("vertical")




        self.data={}
        self.data['kx']=str(1.0)
        self.data['ky']=str(0.0)
        self.data['kz']=str(0.0)

        DirectionsubPane = sc.SizedPanel(pane, -1)
        DirectionsubPane.SetSizerType("horizontal")
        DirectionsubPane.SetSizerProps(expand=True)
        DirectionsubPane.SetExtraStyle(wx.WS_EX_VALIDATE_RECURSIVELY)
        print 'pane'
        WalkTree(pane)
        wx.StaticText(DirectionsubPane, -1, "qx")     
        qx=wx.TextCtrl(DirectionsubPane, -1, self.data['kx'],validator=mFormValidator(self.data,'kx'))
        
        #wx.StaticText(DirectionsubPane, -1, "qy")
        #qy=wx.TextCtrl(DirectionsubPane, -1, self.data['ky'],validator=mFormValidator(self.data,'ky'))

        #wx.StaticText(DirectionsubPane, -1, "qz")
        #qz=wx.TextCtrl(DirectionsubPane, -1, self.data['kz'],validator=mFormValidator(self.data,'kz'))



        # add dialog buttons
        self.SetButtonSizer(self.CreateStdDialogButtonSizer(wx.OK | wx.CANCEL))
        # a little trick to make sure that you can't resize the dialog to
        # less screen space than the controls need
        self.Fit()
        self.SetMinSize(self.GetSize())




if __name__=='__main__':
    #app=MyApp()
    app=wx.PySimpleApp()
    frame=wx.Frame(None,-1,"A Frame",style=wx.DEFAULT_FRAME_STYLE,size=(200,100))
    frame.SetExtraStyle(wx.WS_EX_VALIDATE_RECURSIVELY)
    dlg=FormDialog(parent=frame,id=-1)
    frame.Show()
    result=dlg.ShowModal()
    if result==wx.ID_OK:
        frame.Validate()
        print "OK"
        #dlg.TransferFromWindow()
        print dlg.data
    else:
        print "Cancel"
    
    dlg.Destroy()
    app.MainLoop()
    