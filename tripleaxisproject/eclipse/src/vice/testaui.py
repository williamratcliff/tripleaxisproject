import wx, wx.aui

class TestNotebook(wx.aui.AuiNotebook):
    def __init__(self, *args, **kwargs):
        super(self.__class__, self).__init__(*args, **kwargs)
        self.__auiManager = self.GetAuiManager()
        self.Bind(wx.aui.EVT_AUINOTEBOOK_PAGE_CHANGED, self.OnPageChanged)
        #self.Bind(wx.EVT_CONTEXT_MENU, self.OnContextMenu)
        self.Bind(wx.aui.EVT_AUINOTEBOOK_TAB_RIGHT_DOWN,self.OnTabContextMenu)
        self.tabs = None

    def setContextMenuOnTabs(self):
        # Have to call this out of an event handler, because the
        # AuiTabCtrl doesn't exist yet during __init__
        children = self.GetChildren()
        for widget in children:
            print widget
            if widget.__class__ == wx.aui.AuiTabCtrl:
                print "Found tabs!!!"
                self.tabs = widget
                break
        if self.tabs is not None:
            self.tabs.Bind(wx.EVT_CONTEXT_MENU, self.OnTabContextMenu)
            self.tabs.Bind(wx.EVT_BUTTON, self.OnTabContextMenu)

    def OnPageChanged(self, evt):
        if not self.tabs:
            self.setContextMenuOnTabs()

    def OnContextMenu(self, evt):
        print("Context menu over all notebook window")

    def OnTabContextMenu(self, event):
        print("Context menu over tabs")
        #pos=event.GetPosition()
        #pos=self.ScreenToClient(pos)
        print 'select', event.GetSelection()
        #currpage=self.GetSelection()
        currpage=event.GetSelection()
        print currpage
        
        self.PopupMenu(self.GetPage(currpage).popupmenu)
        #evt.Skip()


if __name__ == "__main__":
    def createPanel(parent, ContentClass, *args, **kwargs):
        panel = wx.Panel(parent)
        content = ContentClass(panel, *args, **kwargs)
        sizer = wx.BoxSizer()
        sizer.Add(content, flag=wx.EXPAND, proportion=1)
        panel.SetSizerAndFit(sizer)
        return panel

    app = wx.App(redirect=False)
    frame = wx.Frame(None)
    notebook = TestNotebook(frame,style=(wx.aui.AUI_NB_DEFAULT_STYLE|wx.aui.AUI_NB_WINDOWLIST_BUTTON)^
                            (wx.aui.AUI_NB_CLOSE_ON_ACTIVE_TAB
                            ))
    for index in range(5):
        page = createPanel(notebook, wx.TextCtrl,
            value='This is page %d'%index, style=wx.TE_MULTILINE)
        page.popupmenu = wx.Menu()
        item1 = page.popupmenu.Append(wx.ID_ANY, "page "+str(index))
        page.popupmenu.AppendSeparator()
        #item2 = page.popupmenu.Append(wx.ID_ANY, "item2 "+str(index))
        notebook.AddPage(page, 'Page %d'%index)
    frame.Show()
    app.MainLoop()