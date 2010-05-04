import wx
import  wx.lib.newevent

#class ClearTreeEvent(wx.PyCommandEvent):
#    def __init__(self,evtType,id):
#        wx.PyCommandEvent.__init__(self,evtType,id)
        

#myEVT_CLEAR_TREE=wx.NewEventType()
#EVT_CLEAR_TREE=wx.PyEventBinder(myEVT_CLEAR_TREE,1)



myEVT_CLEAR_TREE, EVT_CLEAR_TREE = wx.lib.newevent.NewCommandEvent()