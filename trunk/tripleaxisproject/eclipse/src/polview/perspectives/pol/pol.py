"""
This software was developed by the University of Tennessee as part of the
Distributed Data Analysis of Neutron Scattering Experiments (DANSE)
project funded by the US National Science Foundation.

See the license text in license.txt

copyright 2008, University of Tennessee
"""


import wx
import sys
from sans.guicomm.events import EVT_NEW_PLOT
from polarization import FileTreeCtrl4 as FTC
from polarization import pol_catalog4 as polcatalog
if wx.Platform == '__WXMSW__':
    from wx.lib.pdfwin import PDFWindow
    
    

class HelpPanel(wx.Panel):
    def __init__(self, parent):
        wx.Panel.__init__(self, parent, id=-1)
        self.pdf = None

        sizer = wx.BoxSizer(wx.VERTICAL)
        btnSizer = wx.BoxSizer(wx.HORIZONTAL)

        self.pdf = PDFWindow(self, style=wx.SUNKEN_BORDER)
        self.pdf.LoadFile(r'c:\polcorrecter\polcorrecter.pdf')
        sizer.Add(self.pdf, proportion=1, flag=wx.EXPAND)
        
#        btn = wx.Button(self, wx.NewId(), "Open PDF File")
#        self.Bind(wx.EVT_BUTTON, self.OnOpenButton, btn)
#        btnSizer.Add(btn, proportion=1, flag=wx.EXPAND|wx.ALL, border=5)
#
#        btn = wx.Button(self, wx.NewId(), "Previous Page")
#        self.Bind(wx.EVT_BUTTON, self.OnPrevPageButton, btn)
#        btnSizer.Add(btn, proportion=1, flag=wx.EXPAND|wx.ALL, border=5)
#
#        btn = wx.Button(self, wx.NewId(), "Next Page")
#        self.Bind(wx.EVT_BUTTON, self.OnNextPageButton, btn)
#        btnSizer.Add(btn, proportion=1, flag=wx.EXPAND|wx.ALL, border=5)
#
#        btnSizer.Add((50,-1), proportion=2, flag=wx.EXPAND)
#        sizer.Add(btnSizer, proportion=0, flag=wx.EXPAND)

        self.SetSizer(sizer)
        self.SetAutoLayout(True)

    def OnOpenButton(self, event):
        # make sure you have PDF files available on your drive
        dlg = wx.FileDialog(self, wildcard="*.pdf")
        if dlg.ShowModal() == wx.ID_OK:
            wx.BeginBusyCursor()
            self.pdf.LoadFile(dlg.GetPath())
            wx.EndBusyCursor()
        dlg.Destroy()

    def OnPrevPageButton(self, event):
        self.pdf.gotoPreviousPage()

    def OnNextPageButton(self, event):
        self.pdf.gotoNextPage()
    



class Plugin:
    """
        Plug-in class to be instantiated by the GUI manager
    """

    def __init__(self):
        """
            Initialize the plug-in
        """
        ## Plug-in name
        self.sub_menu = "Polarization"

        ## Reference to the parent window
        self.parent = None

        ## List of panels for the simulation perspective (names)
        self.perspective = []

        ## Plot panels
        self.filetree_panel = None


    def populate_menu(self, id, parent):
        """
            Create a 'Plot' menu to list the panels
            available for displaying
            @param id: next available unique ID for wx events
            @param parent: parent window
        """
        self.menu = wx.Menu()
        id = wx.NewId()
        self.menu.Append(id, "Load files", "Load files")
        wx.EVT_MENU(self.parent, id, self.catalog_panel.OnOpen)

        return [(id, self.menu, "Polarization")]

    def help(self,event):
        print 'help!'
        if wx.Platform == '__WXMSW__':
            frame = wx.Frame(None, -1, "PDFWindow", size = (640, 480))
            # make an instance of the class
            HelpPanel(frame)
            # show the frame
            frame.Show(True)
        
    def get_panels(self, parent):
        """
            Create and return a list of panel objects
        """
        ## Save a reference to the parent
        self.parent = parent
        self.filetree_panel = FTC.FileTreePanel(parent, -1)
        self.catalog_panel=polcatalog.CatalogPanel(parent,-1)
        self.catalog_panel.filetree_panel=self.filetree_panel
        self.perspective = [self.filetree_panel.window_name,
                            self.catalog_panel.window_name]

        # We have no initial panels for this plug-in
        return [self.filetree_panel, self.catalog_panel]

    def get_perspective(self):
        """
            Get the list of panel names for this perspective
        """
        return self.perspective

    def on_perspective(self, event):
        """
            Call back function for the perspective menu item.
            We notify the parent window that the perspective
            has changed.
            @param event: menu event
        """
        self.parent.set_perspective(self.perspective)

    def post_init(self):
        """
            Post initialization call back to close the loose ends
            [Somehow openGL needs this call]
        """
        self.parent.set_perspective(self.perspective)

    def _on_show_panel(self, event):
        print "_on_show_panel"
