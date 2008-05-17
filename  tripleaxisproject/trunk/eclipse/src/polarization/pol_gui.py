import wx
import  wx.grid as  gridlib
import sys,os
import classify_files
import numpy as N

class MyApp(wx.App):
    def __init__(self, redirect=False, filename=None, useBestVisual=False, clearSigInt=True):
        wx.App.__init__(self,redirect,filename,clearSigInt)


    def OnInit(self):
        return True




#---------------------------------------------------------------------------

class CustomDataTable(gridlib.PyGridTableBase):
    def __init__(self):
        gridlib.PyGridTableBase.__init__(self)
        self.colLabels = ['off off', 'Hsample', 'Vsample', 'off on','H', 'V', 'on off','H', 'V', 'on on','H', 'V']
        self.rowLabels=['Group 0']
        self.dataTypes = [gridlib.GRID_VALUE_STRING, #off off
                          gridlib.GRID_VALUE_STRING,
                          gridlib.GRID_VALUE_STRING,
                          gridlib.GRID_VALUE_STRING, #off on
                          gridlib.GRID_VALUE_STRING,
                          gridlib.GRID_VALUE_STRING,
                          gridlib.GRID_VALUE_STRING, #on off
                          gridlib.GRID_VALUE_STRING,
                          gridlib.GRID_VALUE_STRING,
                          gridlib.GRID_VALUE_STRING, #on on
                          gridlib.GRID_VALUE_STRING,
                          gridlib.GRID_VALUE_STRING
                          ]
        self.data = []
        self.data.append(['','','',\
                            '','','',\
                            '','','',\
                            '','',''\
                            ])
        return
            #[1010, "The foo doesn't bar", "major", 1, 'MSW', 1, 1, 1, 1.12],
            #[1011, "I've got a wicket in my wocket", "wish list", 2, 'other', 0, 0, 0, 1.50],
            #[1012, "Rectangle() returns a triangle", "critical", 5, 'all', 0, 0, 0, 1.56]
            #]
    #--------------------------------------------------
    # required methods for the wxPyGridTableBase interface
    def GetNumberRows(self):
        return len(self.data) + 0
    def GetNumberCols(self):
        return len(self.colLabels)
    def IsEmptyCell(self, row, col):
        try:
            return not self.data[row][col]
        except IndexError:
            return True
    # Get/Set values in the table.  The Python version of these
    # methods can handle any data-type, (as long as the Editor and
    # Renderer understands the type too,) not just strings as in the
    # C++ version.
    def GetValue(self, row, col):
        try:
            return self.data[row][col]
        except IndexError:
            return ''
    def SetValue(self, row, col, value):
        try:
            self.data[row][col] = value
        except IndexError:
            # add a new row
            self.data.append([''] * self.GetNumberCols())
            self.rowLabels.append('Group '+str(len(self.rowLabels)))
            #print 'setting row ',row,' col ',col, ' val ',value
            #print self.__dict__
            #self.SetValue(row, col, value)
            self.data[row][col] = value
            #print 'set'

            # tell the grid we've added a row
            msg = gridlib.GridTableMessage(self,            # The table
                    gridlib.GRIDTABLE_NOTIFY_ROWS_APPENDED, # what we did to it
                    1                                       # how many
                    )
            self.GetView().ProcessTableMessage(msg)
        return
    #--------------------------------------------------
    # Some optional methods
    # Called when the grid needs to display labels
    def GetColLabelValue(self, col):
        return self.colLabels[col]
    # Called when the grid needs to display labels
    def GetRowLabelValue(self, row):
        return self.rowLabels[row]
    # Called to determine the kind of editor/renderer to use by
    # default, doesn't necessarily have to be the same type used
    # natively by the editor/renderer if they know how to convert.
    def GetTypeName(self, row, col):
        return self.dataTypes[col]
    # Called to determine how the data can be fetched and stored by the
    # editor and renderer.  This allows you to enforce some type-safety
    # in the grid.
    def CanGetValueAs(self, row, col, typeName):
        colType = self.dataTypes[col].split(':')[0]
        if typeName == colType:
            return True
        else:
            return False
    def CanSetValueAs(self, row, col, typeName):
        return self.CanGetValueAs(row, col, typeName)
#---------------------------------------------------------------------------
class CustTableGrid(gridlib.Grid):
    def __init__(self, parent):
        gridlib.Grid.__init__(self, parent, -1)
        table = CustomDataTable()
        # The second parameter means that the grid is to take ownership of the
        # table and will destroy it when done.  Otherwise you would need to keep
        # a reference to it and call it's Destroy method later.
        self.SetTable(table, True)
        #self.SetRowLabelSize(0)
        self.SetMargins(0,0)
        self.AutoSize()
        #gridlib.Grid.SetSelectionMode(self,gridlib.Grid.SelectRows)
        gridlib.Grid.EnableEditing(self,False)
        gridlib.EVT_GRID_CELL_LEFT_DCLICK(self, self.OnLeftDClick)
    # I do this because I don't like the default behaviour of not starting the
    # cell editor on double clicks, but only a second click.
    def OnLeftDClick(self, evt):
        if self.CanEnableCellControl():
            self.EnableCellEditControl()


class CatalogFrame(wx.Frame):
    def __init__(self,parent,id):
        wx.Frame.__init__(self,parent,id,'File Catalog',size=(340,200))
        self.Bind(wx.EVT_CLOSE,self.OnCloseWindow)
        self.createMenuBar()
        p = wx.Panel(self, -1, style=0)
        grid = CustTableGrid(p)
        bs = wx.BoxSizer(wx.VERTICAL)
        bs.Add(grid, 1, wx.GROW|wx.ALL|wx.EXPAND, 5)
        p.SetSizer(bs)
        self.grid=grid

    def menuData(self):
        return(("&File",\
                ("&Open","Open",self.OnOpen),\
                ("&Quit","Quit",self.OnCloseWindow)),\
                )

    def createMenuBar(self):
        menuBar=wx.MenuBar()
        for eachMenuData in self.menuData():
            menuLabel=eachMenuData[0]
            menuItems=eachMenuData[1:]
            menuBar.Append(self.createMenu(menuItems),menuLabel)
            self.SetMenuBar(menuBar)

    def createMenu(self,menuData):
        menu=wx.Menu()
        for eachLabel, eachStatus, eachHandler in menuData:
            if not eachLabel:
                menu.AppendSeparator()
                continue
            menuItem=menu.Append(-1,eachLabel,eachStatus)
            self.Bind(wx.EVT_MENU,eachHandler,menuItem)
        return menu

    def OnCloseWindow(self,event):
        self.Destroy()

    def OnOpen(self,event):
        # Create the dialog. In this case the current directory is forced as the starting
        # directory for the dialog, and no default file name is forced. This can easilly
        # be changed in your program. This is an 'open' dialog, and allows multitple
        # file selections as well.
        #
        # Finally, if the directory is changed in the process of getting files, this
        # dialog is set up to change the current working directory to the path chosen.

        defaultDir=os.getcwd()
        defaultDir=r'c:\bifeo3xtal\jan8_2008\9175'
        wildcard="bt7 files (*.bt7)|*.bt7|All files (*.*)|*.*"
        dlg = wx.FileDialog(
            self, message="Choose a file",
            defaultDir=defaultDir,
            defaultFile="",
            wildcard=wildcard,
            style=wx.OPEN | wx.MULTIPLE | wx.CHANGE_DIR
            )

        # Show the dialog and retrieve the user response. If it is the OK response,
        # process the data.
        if dlg.ShowModal() == wx.ID_OK:
            # This returns a Python list of files that were selected.
            paths = dlg.GetPaths()

            #self.log.WriteText('You selected %d files:' % len(paths))
            self.files=paths
            #for path in paths:
            #    self.log.WriteText('           %s\n' % path)

        # Compare this with the debug above; did we change working dirs?
        #self.log.WriteText("CWD: %s\n" % os.getcwd())


        self.catalog=classify_files.readfiles(self.files)
        #print 'opening:'
        #print self.catalog.pm.files
        for row in range(len(self.catalog.mm.files)):
            gridlib.Grid.SetCellValue(self.grid,row,0,self.catalog.mm.files[row])
            gridlib.Grid.SetCellValue(self.grid,row,1,str(self.catalog.mm.data[row]['hsample']))
            gridlib.Grid.SetCellValue(self.grid,row,2,str(self.catalog.mm.data[row]['vsample']))

        for row in range(len(self.catalog.mp.files)):
            gridlib.Grid.SetCellValue(self.grid,row,3,self.catalog.mp.files[row])
            gridlib.Grid.SetCellValue(self.grid,row,4,str(self.catalog.mp.data[row]['hsample']))
            gridlib.Grid.SetCellValue(self.grid,row,5,str(self.catalog.mp.data[row]['vsample']))
        for row in range(len(self.catalog.pm.files)):
            gridlib.Grid.SetCellValue(self.grid,row,6,self.catalog.pm.files[row])
            gridlib.Grid.SetCellValue(self.grid,row,7,str(self.catalog.pm.data[row]['hsample']))
            gridlib.Grid.SetCellValue(self.grid,row,8,str(self.catalog.pm.data[row]['vsample']))

        for row in range(len(self.catalog.pp.files)):
            gridlib.Grid.SetCellValue(self.grid,row,9,self.catalog.pp.files[row])
            gridlib.Grid.SetCellValue(self.grid,row,10,str(self.catalog.pp.data[row]['hsample']))
            gridlib.Grid.SetCellValue(self.grid,row,11,str(self.catalog.pp.data[row]['vsample']))

        gridlib.Grid.AutoSize(self.grid)
        gridlib.Grid.ForceRefresh(self.grid)
        # Destroy the dialog. Don't do this until you are done with it!
        # BAD things can happen otherwise!
        dlg.Destroy()






if __name__=='__main__':
    app=MyApp()
    frame=CatalogFrame(parent=None,id=-1)
    frame.Show()
    app.MainLoop()