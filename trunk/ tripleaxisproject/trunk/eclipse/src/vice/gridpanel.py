import  wx
import  wx.grid as  gridlib
import sys,os
from polarization import classify_files2 as classify_files
from utilities import readncnr3 as readncnr
import numpy as N



keyMap = {
    wx.WXK_BACK : "WXK_BACK",
    wx.WXK_TAB : "WXK_TAB",
    wx.WXK_RETURN : "WXK_RETURN",
    wx.WXK_ESCAPE : "WXK_ESCAPE",
    wx.WXK_SPACE : "WXK_SPACE",
    wx.WXK_DELETE : "WXK_DELETE",
    wx.WXK_START : "WXK_START",
    wx.WXK_LBUTTON : "WXK_LBUTTON",
    wx.WXK_RBUTTON : "WXK_RBUTTON",
    wx.WXK_CANCEL : "WXK_CANCEL",
    wx.WXK_MBUTTON : "WXK_MBUTTON",
    wx.WXK_CLEAR : "WXK_CLEAR",
    wx.WXK_SHIFT : "WXK_SHIFT",
    wx.WXK_ALT : "WXK_ALT",
    wx.WXK_CONTROL : "WXK_CONTROL",
    wx.WXK_MENU : "WXK_MENU",
    wx.WXK_PAUSE : "WXK_PAUSE",
    wx.WXK_CAPITAL : "WXK_CAPITAL",
    wx.WXK_PRIOR : "WXK_PRIOR",
    wx.WXK_NEXT : "WXK_NEXT",
    wx.WXK_END : "WXK_END",
    wx.WXK_HOME : "WXK_HOME",
    wx.WXK_LEFT : "WXK_LEFT",
    wx.WXK_UP : "WXK_UP",
    wx.WXK_RIGHT : "WXK_RIGHT",
    wx.WXK_DOWN : "WXK_DOWN",
    wx.WXK_SELECT : "WXK_SELECT",
    wx.WXK_PRINT : "WXK_PRINT",
    wx.WXK_EXECUTE : "WXK_EXECUTE",
    wx.WXK_SNAPSHOT : "WXK_SNAPSHOT",
    wx.WXK_INSERT : "WXK_INSERT",
    wx.WXK_HELP : "WXK_HELP",
    wx.WXK_NUMPAD0 : "WXK_NUMPAD0",
    wx.WXK_NUMPAD1 : "WXK_NUMPAD1",
    wx.WXK_NUMPAD2 : "WXK_NUMPAD2",
    wx.WXK_NUMPAD3 : "WXK_NUMPAD3",
    wx.WXK_NUMPAD4 : "WXK_NUMPAD4",
    wx.WXK_NUMPAD5 : "WXK_NUMPAD5",
    wx.WXK_NUMPAD6 : "WXK_NUMPAD6",
    wx.WXK_NUMPAD7 : "WXK_NUMPAD7",
    wx.WXK_NUMPAD8 : "WXK_NUMPAD8",
    wx.WXK_NUMPAD9 : "WXK_NUMPAD9",
    wx.WXK_MULTIPLY : "WXK_MULTIPLY",
    wx.WXK_ADD : "WXK_ADD",
    wx.WXK_SEPARATOR : "WXK_SEPARATOR",
    wx.WXK_SUBTRACT : "WXK_SUBTRACT",
    wx.WXK_DECIMAL : "WXK_DECIMAL",
    wx.WXK_DIVIDE : "WXK_DIVIDE",
    wx.WXK_F1 : "WXK_F1",
    wx.WXK_F2 : "WXK_F2",
    wx.WXK_F3 : "WXK_F3",
    wx.WXK_F4 : "WXK_F4",
    wx.WXK_F5 : "WXK_F5",
    wx.WXK_F6 : "WXK_F6",
    wx.WXK_F7 : "WXK_F7",
    wx.WXK_F8 : "WXK_F8",
    wx.WXK_F9 : "WXK_F9",
    wx.WXK_F10 : "WXK_F10",
    wx.WXK_F11 : "WXK_F11",
    wx.WXK_F12 : "WXK_F12",
    wx.WXK_F13 : "WXK_F13",
    wx.WXK_F14 : "WXK_F14",
    wx.WXK_F15 : "WXK_F15",
    wx.WXK_F16 : "WXK_F16",
    wx.WXK_F17 : "WXK_F17",
    wx.WXK_F18 : "WXK_F18",
    wx.WXK_F19 : "WXK_F19",
    wx.WXK_F20 : "WXK_F20",
    wx.WXK_F21 : "WXK_F21",
    wx.WXK_F22 : "WXK_F22",
    wx.WXK_F23 : "WXK_F23",
    wx.WXK_F24 : "WXK_F24",
    wx.WXK_NUMLOCK : "WXK_NUMLOCK",
    wx.WXK_SCROLL : "WXK_SCROLL",
    wx.WXK_PAGEUP : "WXK_PAGEUP",
    wx.WXK_PAGEDOWN : "WXK_PAGEDOWN",
    wx.WXK_NUMPAD_SPACE : "WXK_NUMPAD_SPACE",
    wx.WXK_NUMPAD_TAB : "WXK_NUMPAD_TAB",
    wx.WXK_NUMPAD_ENTER : "WXK_NUMPAD_ENTER",
    wx.WXK_NUMPAD_F1 : "WXK_NUMPAD_F1",
    wx.WXK_NUMPAD_F2 : "WXK_NUMPAD_F2",
    wx.WXK_NUMPAD_F3 : "WXK_NUMPAD_F3",
    wx.WXK_NUMPAD_F4 : "WXK_NUMPAD_F4",
    wx.WXK_NUMPAD_HOME : "WXK_NUMPAD_HOME",
    wx.WXK_NUMPAD_LEFT : "WXK_NUMPAD_LEFT",
    wx.WXK_NUMPAD_UP : "WXK_NUMPAD_UP",
    wx.WXK_NUMPAD_RIGHT : "WXK_NUMPAD_RIGHT",
    wx.WXK_NUMPAD_DOWN : "WXK_NUMPAD_DOWN",
    wx.WXK_NUMPAD_PRIOR : "WXK_NUMPAD_PRIOR",
    wx.WXK_NUMPAD_PAGEUP : "WXK_NUMPAD_PAGEUP",
    wx.WXK_NUMPAD_NEXT : "WXK_NUMPAD_NEXT",
    wx.WXK_NUMPAD_PAGEDOWN : "WXK_NUMPAD_PAGEDOWN",
    wx.WXK_NUMPAD_END : "WXK_NUMPAD_END",
    wx.WXK_NUMPAD_BEGIN : "WXK_NUMPAD_BEGIN",
    wx.WXK_NUMPAD_INSERT : "WXK_NUMPAD_INSERT",
    wx.WXK_NUMPAD_DELETE : "WXK_NUMPAD_DELETE",
    wx.WXK_NUMPAD_EQUAL : "WXK_NUMPAD_EQUAL",
    wx.WXK_NUMPAD_MULTIPLY : "WXK_NUMPAD_MULTIPLY",
    wx.WXK_NUMPAD_ADD : "WXK_NUMPAD_ADD",
    wx.WXK_NUMPAD_SEPARATOR : "WXK_NUMPAD_SEPARATOR",
    wx.WXK_NUMPAD_SUBTRACT : "WXK_NUMPAD_SUBTRACT",
    wx.WXK_NUMPAD_DECIMAL : "WXK_NUMPAD_DECIMAL",
    wx.WXK_NUMPAD_DIVIDE : "WXK_NUMPAD_DIVIDE"
    }


class MyApp(wx.App):
    def __init__(self, redirect=False, filename=None, useBestVisual=False, clearSigInt=True):
        wx.App.__init__(self,redirect,filename,clearSigInt)


    def OnInit(self):
        return True




#---------------------------------------------------------------------------

class CustomDataTable(gridlib.PyGridTableBase):
    def __init__(self):
        gridlib.PyGridTableBase.__init__(self)
        self.colLabels = ['Select?', 'filename','seq #', 'polarization state','hsample','vsample','h','k','l','e','a3','a4','temp','magfield']
        self.rowLabels=['File 0']

        self.dataTypes = [gridlib.GRID_VALUE_STRING, #selected
                          gridlib.GRID_VALUE_STRING,#filename
                          gridlib.GRID_VALUE_STRING,#sequence number
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
        self.data.append(['', #selected
                        '', #filename
                        '', #sequence number
                        '', #polarization state
                        '', #hsample
                        '', #vsample
                        gridbar.Bar(0,0,(0,10)), #h  (lo,high, range)
                        gridbar.Bar(0,0,(0,10)), #k
                        gridbar.Bar(0,0,(0,10)), #l
                        gridbar.Bar(0,0,(0,10)), #e
                        gridbar.Bar(0,0,(0,10)), #a3
                        gridbar.Bar(0,0,(0,10)), #a4
                        gridbar.Bar(0,0,(0,10)), #temp
                        gridbar.Bar(0,0,(0,10)), #magfield
                        ])
        return
            #[1010, "The foo doesn't bar", "major", 1, 'MSW', 1, 1, 1, 1.12],
            #[1011, "I've got a wicket in my wocket", "wish list", 2, 'other', 0, 0, 0, 1.50],
            #[1012, "Rectangle() returns a triangle", "critical", 5, 'all', 0, 0, 0, 1.56]
            #]
    #--------------------------------------------------
    # required methods for the wxPyGridTableBase interface
    def GetNumberRows(self):
        return len(self.data)
        #return len(self.data)
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
            #print 'SetValue works',self.GetNumberRows(),self.data[row][1]
        except IndexError:
            # add a new row
            #print 'IndexError in SetValue',self.GetNumberRows()
            self.AppendRow()
            self.data[row][col]=value
            #print 'IndexError in SetValue after SetValue',self.GetNumberRows()
            #print 'setting row ',row,' col ',col, ' val ',value
            #print self.__dict__
            #self.SetValue(row, col, value)
        return

    def AppendRow(self):
            self.data.append([''] * self.GetNumberCols())
            #print 'After Append SetValue',self.GetNumberRows()
            #self.rowLabels[row]='File '+str(len(self.rowLabels))
            #self.rowLabels.append('File '+str(len(self.rowLabels)))

            # tell the grid we've added a row
            msg = gridlib.GridTableMessage(self,            # The table
                    gridlib.GRIDTABLE_NOTIFY_ROWS_APPENDED, # what we did to it
                    1                                       # how many
                    )
            #print 'size notified',self.GetNumberRows()
            self.GetView().ProcessTableMessage(msg)
            #print 'self.rowLabels', self.rowLabels
            #self.data[row][col] = value


    #--------------------------------------------------
    # Some optional methods
    # Called when the grid needs to display labels
    def GetColLabelValue(self, col):
        return self.colLabels[col]
    # Called when the grid needs to display labels
    def GetRowLabelValue(self, row):
        return 'File '+str(row)
        #return self.rowLabels[row]
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

    def DeleteRows(self,pos=0,numRows=1):
#        print 'Delete number',self.GetNumberRows()
#        print 'pos',pos
#        print 'numRows', numRows
        if numRows>=0 and numRows<=self.GetNumberRows():
#            print 'Delete',numRows
            #for i in range(numRows):
            #    self.data.pop()
            del self.data[pos:pos+numRows]
            msg = gridlib.GridTableMessage(self,            # The table
            gridlib.GRIDTABLE_NOTIFY_ROWS_DELETED, # what we did to it
            pos,numRows                                     # how many
            )
            #msg = wx.grid.GridTableMessage(self, 0, numRows)
            self.GetView().ProcessTableMessage(msg)
            
#            print 'Deleted'
            self.UpdateValues()
            return True
        else:
            return False
        
    def UpdateValues( self ):
            """Update all displayed values"""
            msg =gridlib.GridTableMessage(self, gridlib.GRIDTABLE_REQUEST_VIEW_GET_VALUES)
            self.GetView().ProcessTableMessage(msg)
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
        gridlib.Grid.EnableEditing(self,False)
        attr=gridlib.GridCellAttr()
        attr.SetReadOnly(True)
        self.SetColAttr(0,attr)
        for col in range(1,14):
            attr=gridlib.GridCellAttr()
            attr.SetReadOnly(True)
            #attr.SetBackgroundColour('grey' if col%2 else (139, 139, 122))
            #attr.SetTextColour((167,167,122) if col%2 else (139, 139, 122))
            self.SetColAttr(col,attr)
        #gridlib.EVT_GRID_CELL_LEFT_DCLICK(self, self.OnLeftDClick)
        gridlib.EVT_GRID_CELL_LEFT_CLICK(self,self.OnLeftClick)
        #gridlib.EVT_GRID_CELL_CHANGE(self,self.OnCellChange)
        gridlib.EVT_GRID_LABEL_LEFT_DCLICK(self,self.onLeftDClickRowCell)

    # I do this because I don't like the default behaviour of not starting the
    # cell editor on double clicks, but only a second click.
    def OnLeftClick(self, evt):
        print 'LeftClick'
        col=evt.GetCol()
        row=evt.GetRow()
        table=self.GetTable()
        if col<=0 and row >=0:
            currval=table.GetValue(row,0)
            if currval=='':
                table.SetValue(row,0,'x')
            else:
                table.SetValue(row,0,'')


        #if self.CanEnableCellControl():
        #    self.EnableCellEditControl()
        gridlib.Grid.ForceRefresh(self)

    def OnCellChange(self, evt):
#        print 'Changed'
        if self.CanEnableCellControl():
            self.EnableCellEditControl()
        gridlib.Grid.ForceRefresh(self)
        evt.Skip()



    def onLeftDClickRowCell(self,evt):
        col=evt.GetCol()
        table=self.GetTable()
        data=N.array(table.data)
#        print 'before ', data[:,0]
        col_to_sort=[(i,s) for i,s in enumerate(data[:,col])]
        col_to_sort.sort(lambda x,y: cmp(x[1],y[1]))
        g_col = [i for (i,s) in col_to_sort]
        #print col_to_sort
        if col >=0:
            if (N.diff(g_col)>0).all():
                g_col=g_col[::-1]

            #print 'col=',col
            #print 'sort '
            #print g
            for i in range(data.shape[1]):
                data[:,i]=data[g_col,i]
            table.data=data.tolist()
#            print 'after',data[:,0]
            gridlib.Grid.AutoSize(self)
            gridlib.Grid.ForceRefresh(self)
        #evt.Skip()
        
        
if __name__=='__main__':
    frame=