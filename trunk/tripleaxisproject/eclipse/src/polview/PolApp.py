"""
    Dummy application.
    Allows the user to set an external data manager
"""

#import gui_manager
from sans.guiframe import gui_manager

# For py2exe, import config here
import local_config
import sys

class SansView():
    
    def __init__(self):
        """
        
        """
        #from gui_manager import ViewApp
        self.gui = gui_manager.ViewApp(0)
        #self.gui = gui_manager.ViewApp()
        
        # Add perspectives to the basic application
        # Additional perspectives can still be loaded
        # dynamically
        import perspectives.pol as module
        plug = module.Plugin()
        self.gui.add_perspective(plug)
            
        # Build the GUI
        self.gui.build_gui()
        
        # Set the application manager for the GUI
        self.gui.set_manager(self)
        
        # Start the main loop
        self.gui.MainLoop()  

if __name__ == "__main__": 
    sansview = SansView()