import wx
import  wx.grid as  gridlib
import sys,os
import classify_files2 as classify_files
import numpy as N
import gridbar

class MyApp(wx.App):
    def __init__(self, redirect=False, filename=None, useBestVisual=False, clearSigInt=True):
        wx.App.__init__(self,redirect,filename,clearSigInt)


    def OnInit(self):
        return True




#---------------------------------------------------------------------------

class CustomDataTable(gridlib.PyGridTableBase):
    def __init__(self):
        gridlib.PyGridTableBase.__init__(self)
        self.colLabels = ['Select?', 'filename', 'polarization state','hsample','vsample','h','k','l','e','a3','a4','temp','magfield']
        self.rowLabels=['Group 0']

        self.dataTypes = [gridlib.GRID_VALUE_BOOL, #selected
                          gridlib.GRID_VALUE_STRING,#filename
                          gridlib.GRID_VALUE_STRING, #polarization state
                          gridlib.GRID_VALUE_STRING, #hsample
                          gridlib.GRID_VALUE_STRING, #vsample
                          gridbar.GRID_VALUE_BAR, #h
                          gridbar.GRID_VALUE_BAR, #k
                          gridbar.GRID_VALUE_BAR, #l
                          gridbar.GRID_VALUE_BAR, #e
                          gridbar.GRID_VALUE_BAR, #a3
                          gridbar.GRID_VALUE_BAR, #a4
                          gridbar.GRID_VALUE_BAR, #temp
                          gridbar.GRID_VALUE_BAR, #magfield
                          #gridlib.GRID_VALUE_STRING,
                          ]
        self.data = []
        self.data.append([1, #selected
                        '', #filename
                        '', #polarization state
                        '', #hsample
                        '', #vsample
                        gridbar.Bar(5,6,(0,10)), #h
                        gridbar.Bar(5,7,(0,10)), #k
                        gridbar.Bar(5,8,(0,10)), #l
                        gridbar.Bar(5,9,(0,10)), #e
                        gridbar.Bar(5,10,(0,10)), #a3
                        gridbar.Bar(5,11,(0,10)), #a4
                        gridbar.Bar(5,12,(0,10)), #temp
                        gridbar.Bar(5,13,(0,10)), #magfield
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

        gridbar.register(self)
        # The second parameter means that the grid is to take ownership of the
        # table and will destroy it when done.  Otherwise you would need to keep
        # a reference to it and call it's Destroy method later.
        self.SetTable(table, True)
        #attr = gridlib.GridCellAttr()
        #attr.SetReadOnly(True)
        #attr.SetRenderer(gridbar.GridCellBarRenderer())
        #self.SetColAttr(13, attr)
        #self.SetCellValue(1,13,'q')
        #self.SetCellRenderer(1,13,gridbar.GridCellBarRenderer)
        #self.SetRowLabelSize(0)
        self.SetMargins(0,0)
        self.AutoSize()
        #gridlib.Grid.SetSelectionMode(self,gridlib.Grid.SelectRows)
        gridlib.Grid.EnableEditing(self,True)
        attr=gridlib.GridCellAttr()
        attr.SetReadOnly(False)
        self.SetColAttr(0,attr)
        for col in range(1,13):
            attr=gridlib.GridCellAttr()
            attr.SetReadOnly(True)
            self.SetColAttr(col,attr)
        gridlib.EVT_GRID_CELL_LEFT_DCLICK(self, self.OnLeftDClick)
        gridlib.EVT_GRID_LABEL_LEFT_DCLICK(self,self.onLeftDClickRowCell)

    # I do this because I don't like the default behaviour of not starting the
    # cell editor on double clicks, but only a second click.
    def OnLeftDClick(self, evt):
        if self.CanEnableCellControl():
            self.EnableCellEditControl()

    def onLeftDClickRowCell(self,evt):
        col=evt.GetCol()
        print 'col=',col


class CatalogFrame(wx.Frame):
    def __init__(self,parent,id):
        wx.Frame.__init__(self,parent,id,'File Catalog',size=(640,200))
        self.Bind(wx.EVT_CLOSE,self.OnCloseWindow)
        self.createMenuBar()
        p = wx.Panel(self, -1, style=0)
        grid = CustTableGrid(p)
        bs = wx.BoxSizer(wx.VERTICAL)
        bs.Add(grid, 1, wx.GROW|wx.ALL|wx.EXPAND, 5)
        p.SetSizer(bs)
        #p.Fit()
        self.grid=grid
        self.bs=bs
        self.tooltip = ''
        self.grid.GetGridWindow().Bind(wx.EVT_MOTION, self.onMouseOver)
        #self.grid.GetGridColLabelWindow().Bind(wx.EVT_MOTION, self.onLeftDclick)

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
            table=self.grid.GetTable()
            table.data=[]
            #table.Update()
            for row in range(len(self.catalog.files)):
                gridlib.Grid.SetCellValue(self.grid,row,1,self.catalog.files[row])
                gridlib.Grid.SetCellValue(self.grid,row,2,str(self.catalog.data[row]['polarization state']))
                gridlib.Grid.SetCellValue(self.grid,row,3,str(self.catalog.data[row]['hsample']))
                gridlib.Grid.SetCellValue(self.grid,row,4,str(self.catalog.data[row]['vsample']))
                if self.catalog.data[row].has_key('h'):
                    range_column=(self.catalog.h_range.min,self.catalog.h_range.max)
                    range_cell=(self.catalog.data[row]['h']['min'],self.catalog.data[row]['h']['max'])
                    currbar=gridbar.Bar(self.catalog.data[row]['h']['min'],self.catalog.data[row]['h']['max'],range_column)
                    table.SetValue(row,5,currbar)
                if self.catalog.data[row].has_key('k'):
                    range_column=(self.catalog.k_range.min,self.catalog.k_range.max)
                    range_cell=(self.catalog.data[row]['k']['min'],self.catalog.data[row]['k']['max'])
                    currbar=gridbar.Bar(self.catalog.data[row]['k']['min'],self.catalog.data[row]['k']['max'],range_column)
                    table.SetValue(row,6,currbar)
                if self.catalog.data[row].has_key('l'):
                    range_column=(self.catalog.l_range.min,self.catalog.l_range.max)
                    range_cell=(self.catalog.data[row]['l']['min'],self.catalog.data[row]['l']['max'])
                    currbar=gridbar.Bar(self.catalog.data[row]['l']['min'],self.catalog.data[row]['l']['max'],range_column)
                    table.SetValue(row,7,currbar)
                if self.catalog.data[row].has_key('e'):
                    range_column=(self.catalog.e_range.min,self.catalog.e_range.max)
                    range_cell=(self.catalog.data[row]['e']['min'],self.catalog.data[row]['e']['max'])
                    currbar=gridbar.Bar(self.catalog.data[row]['e']['min'],self.catalog.data[row]['e']['max'],range_column)
                    table.SetValue(row,8,currbar)
                if self.catalog.data[row].has_key('a3'):
                    range_column=(self.catalog.a3_range.min,self.catalog.a3_range.max)
                    range_cell=(self.catalog.data[row]['a3']['min'],self.catalog.data[row]['a3']['max'])
                    currbar=gridbar.Bar(self.catalog.data[row]['a3']['min'],self.catalog.data[row]['a3']['max'],range_column)
                    table.SetValue(row,9,currbar)
                if self.catalog.data[row].has_key('a4'):
                    range_column=(self.catalog.a4_range.min,self.catalog.a4_range.max)
                    range_cell=(self.catalog.data[row]['a4']['min'],self.catalog.data[row]['a4']['max'])
                    currbar=gridbar.Bar(self.catalog.data[row]['a4']['min'],self.catalog.data[row]['a4']['max'],range_column)
                    table.SetValue(row,10,currbar)
                if self.catalog.data[row].has_key('temp'):
                    #print 'temp'
                    range_column=(self.catalog.temp_range.min,self.catalog.temp_range.max)
                    range_cell=(self.catalog.data[row]['temp']['min'],self.catalog.data[row]['temp']['max'])
                    currbar=gridbar.Bar(self.catalog.data[row]['temp']['min'],self.catalog.data[row]['temp']['max'],range_column)
                    table.SetValue(row,11,currbar)
                if self.catalog.data[row].has_key('magfield'):
                    range_column=(self.catalog.magfield_range.min,self.catalog.magfield_range.max)
                    range_cell=(self.catalog.data[row]['magfield']['min'],self.catalog.data[row]['magfield']['max'])
                    currbar=gridbar.Bar(self.catalog.data[row]['magfield']['min'],self.catalog.data[row]['magfield']['max'],range_column)
                    table.SetValue(row,12,currbar)


            #table.Update()
            gridlib.Grid.AutoSize(self.grid)
            gridlib.Grid.ForceRefresh(self.grid)
        # Destroy the dialog. Don't do this until you are done with it!
        # BAD things can happen otherwise!
        dlg.Destroy()


    def onMouseOver(self, event):
        '''
        Method to calculate where the mouse is pointing and
        then set the tooltip dynamically.
        '''

        # Use CalcUnscrolledPosition() to get the mouse position within the
        # entire grid including what's offscreen
        x, y =self.grid.CalcUnscrolledPosition(event.GetX(),event.GetY())

        coords = self.grid.XYToCell(x, y)
        #coords = grid.XYToCell(x, y)
        col = coords[1]
        table=self.grid.GetTable()
        # Example colum limit to apply the custom tooltip to
        if col>4:
            row = coords[0]
            bar = table.GetValue(row, col)
            try:
                low=bar.low
                high=bar.high
                event.GetEventObject().SetToolTipString('range=(%4.3f,%4.3f)'%(low,high))
                self.tooltip='range=(%4.3f,%4.3f)' %(low,high)
            except:
                event.GetEventObject().SetToolTipString('')
                self.tooltip = ''
        else:
            event.GetEventObject().SetToolTipString('')
            self.tooltip = ''




if __name__=='__main__':
    app=MyApp()
    frame=CatalogFrame(parent=None,id=-1)
    frame.Show()
    app.MainLoop()