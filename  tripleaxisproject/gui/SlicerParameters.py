"""
    Panel class to show the slicer parameters 
"""

import wx
import wx.lib.newevent
from copy import deepcopy

(SlicerEvent, EVT_SLICER)   = wx.lib.newevent.NewEvent()
(SlicerParameterEvent, EVT_SLICER_PARS)   = wx.lib.newevent.NewEvent()

class SlicerParameterPanel(wx.Panel):
    #TODO: show units
    #TODO: order parameters properly
    
    def __init__(self, parent, *args, **kwargs):
        wx.Panel.__init__(self, parent, *args, **kwargs)
        self.params = {}
        self.parent = parent
        self.type = None
        self.listeners = []
        self.parameters = []
        self.bck = wx.GridBagSizer(5,5)
        self.SetSizer(self.bck)
               
        title = wx.StaticText(self, -1, "Right-click on 2D plot for slicer options", style=wx.ALIGN_LEFT)
        self.bck.Add(title, (0,0), (1,2), flag=wx.LEFT|wx.ALIGN_CENTER_VERTICAL, border=15)
        
        # Bindings
        self.parent.Bind(EVT_SLICER, self.onEVT_SLICER)
        self.parent.Bind(EVT_SLICER_PARS, self.onParamChange)

    def onEVT_SLICER(self, event):
        """
            Process EVT_SLICER events
            When the slicer changes, update the panel
            
            @param event: EVT_SLICER event
        """
        event.Skip()
        self.set_slicer(event.type, event.params)
        
    def set_slicer(self, type, params):
        """
            Rebuild the panel
        """
        self.bck.Clear(True)  
        self.type = type  
        
        if type==None:
            title = wx.StaticText(self, -1, "Right-click on 2D plot for slicer options", style=wx.ALIGN_LEFT)
            self.bck.Add(title, (0,0), (1,2), flag=wx.LEFT|wx.ALIGN_CENTER_VERTICAL, border=15)

        else:
            title = wx.StaticText(self, -1, "Slicer Parameters", style=wx.ALIGN_LEFT)
            self.bck.Add(title, (0,0), (1,2), flag=wx.LEFT|wx.ALIGN_CENTER_VERTICAL, border=15)
            
            n = 1
            self.parameters = []
            #params = slicer.get_params()
            keys = params.keys()
            keys.sort()
            
            for item in keys:
                n += 1
                text = wx.StaticText(self, -1, item, style=wx.ALIGN_LEFT)
                self.bck.Add(text, (n-1,0), flag = wx.LEFT|wx.ALIGN_CENTER_VERTICAL, border = 15)
                ctl = wx.TextCtrl(self, -1, size=(80,20), style=wx.TE_PROCESS_ENTER)
                
                ctl.SetToolTipString("Modify the value of %s to change the 2D slicer" % item)
                
                
                ctl.SetValue(str(params[item]))
                self.Bind(wx.EVT_TEXT_ENTER, self.onTextEnter)
                ctl.Bind(wx.EVT_KILL_FOCUS, self.onTextEnter)
                self.parameters.append([item, ctl])
                self.bck.Add(ctl, (n-1,1), flag=wx.TOP|wx.BOTTOM, border = 0)

        self.bck.Layout()
        self.bck.Fit(self)
        self.parent.GetSizer().Layout()

    def onParamChange(self, evt):
        evt.Skip()
        if evt.type == "UPDATE":
            for item in self.parameters:              
                if item[0] in evt.params:
                    item[1].SetValue("%-5.3g" %evt.params[item[0]])
                    item[1].Refresh()
        
    def onTextEnter(self, evt): 
        """
            Parameters have changed
        """ 
        params = {}
        for item in self.parameters:
            try:
                params[item[0]] = float(item[1].GetValue())
            except:
                #TODO: popup
                print "Bad value for parameter", item[0], item[1]

        # Post parameter event
        event = SlicerParameterEvent(type=self.type, params=params)
        wx.PostEvent(self.parent, event)
        