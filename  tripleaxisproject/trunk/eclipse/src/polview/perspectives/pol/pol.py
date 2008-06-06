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
