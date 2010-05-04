import wx

class MyFrame(wx.Frame):
    def __init__(self):
        wx.Frame.__init__(self,None,-1,"Sub-menu Example")
        p=wx.Panel(self)
        menu=wx.Menu()
        
        submenu=wx.Menu()
        submenu.Append(-1,"sub-item-1")
        submenu.Append(-1,"sub-item-2")
        menu.AppendMenu(-1,'sub-menu',submenu)
        
        
        subsubmenu=wx.Menu()
        print subsubmenu
        print subsubmenu.__dict__
        #submenu.GetInvokingWindow()
        subsubevt=subsubmenu.Append(-1,"sub-sub-item-1")
        subsubmenu.Append(-1,"sub-sub-item-2")
        submenu.AppendMenu(-1,'subsub-menu',subsubmenu)
        self.Bind(wx.EVT_MENU,self.OnSubSub,subsubevt)
        
        
        menu.AppendSeparator()
        exit=menu.Append(-1,"Exit")
        submenu.Bind(wx.EVT_MENU, self.OnExit,exit)
        
        menuBar=wx.MenuBar()
        menuBar.Append(menu,"Menu")
        self.SetMenuBar(menuBar)
        
        
    def OnExit(self,event):
        self.Close()
        
    def OnSubSub(self,event):
        
        menubar=self.GetMenuBar()
        
        item=menubar.FindItemById(event.GetId())
        text=item.GetText()
        wx.MessageBox("You selected the '%s' item" %text)
        print 'bye'
        
        
        
if __name__=="__main__":
    app=wx.PySimpleApp()
    frame=MyFrame()
    frame.Show()
    app.MainLoop()
        
