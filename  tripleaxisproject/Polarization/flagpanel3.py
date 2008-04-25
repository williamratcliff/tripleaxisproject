import wx
import wxaddons.sized_controls as sc
import  wx.lib.intctrl
import  wx.grid as  gridlib

class MyApp(wx.App):
    def __init__(self, redirect=False, filename=None, useBestVisual=False, clearSigInt=True):
        wx.App.__init__(self,redirect,filename,clearSigInt)


    def OnInit(self):
        return True


class ConstraintMatrixTable(gridlib.PyGridTableBase):
    def __init__(self):
        gridlib.PyGridTableBase.__init__(self)
        self.colLabels = ['off off', 'off on', 'on off','on on']
        self.rowLabels=['off off', 'off on', 'on off','on on']

        self.dataTypes = [gridlib.GRID_VALUE_FLOAT, #off off
                          gridlib.GRID_VALUE_FLOAT,#off on
                          gridlib.GRID_VALUE_FLOAT,#on off
                          gridlib.GRID_VALUE_FLOAT,#on on
                          ]
        self.data = []
#        self.data.append([1, #selected
#                        '', #filename
#                        '', #sequence number
#                        '', #polarization state
#                        '', #hsample
#                        '', #vsample
#                        ])
        return
            #[1010, "The foo doesn't bar", "major", 1, 'MSW', 1, 1, 1, 1.12],
            #[1011, "I've got a wicket in my wocket", "wish list", 2, 'other', 0, 0, 0, 1.50],
            #[1012, "Rectangle() returns a triangle", "critical", 5, 'all', 0, 0, 0, 1.56]
            #]
    #--------------------------------------------------
    # required methods for the wxPyGridTableBase interface
    def GetNumberRows(self):
        return len(self.data)+0
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
            #self.data[row][col] = value
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
class ConstraintMatrixGrid(gridlib.Grid):
    def __init__(self, parent):
        gridlib.Grid.__init__(self, parent, -1)
        table = ConstraintMatrixTable()
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
        self.sorted=False
        #gridlib.Grid.SetSelectionMode(self,gridlib.Grid.SelectRows)
        gridlib.Grid.EnableEditing(self,True)
        attr=gridlib.GridCellAttr()
        attr.SetReadOnly(False)
        self.SetColAttr(0,attr)
        for col in range(1,14):
            attr=gridlib.GridCellAttr()
            attr.SetReadOnly(True)
            #attr.SetBackgroundColour('grey' if col%2 else (139, 139, 122))
            #attr.SetTextColour((167,167,122) if col%2 else (139, 139, 122))
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
        table=self.GetTable()
        data=N.array(table.data)
        #print data.shape
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
            for i in range(data.shape[1]-1):
                data[:,i]=data[g_col,i]
            table.data=data
            gridlib.Grid.AutoSize(self)
            gridlib.Grid.ForceRefresh(self)











class FormDialog(sc.SizedDialog):
    def __init__(self, parent, id):
        sc.SizedDialog.__init__(self, None, -1, "Reduction Forms",
                        style=wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER)

        pane = self.GetContentsPane()
        pane.SetSizerType("vertical")



        monitorPane = sc.SizedPanel(pane, -1)
        monitorPane.SetSizerType("horizontal")
        monitorPane.SetSizerProps(expand=True)
        # row 1
        wx.StaticText(monitorPane, -1, "Monitor Position")
        wx.CheckBox(monitorPane, -1, "PrePolarizer")
        wx.CheckBox(monitorPane, -1, "PostPolarizer")
        #textCtrl = wx.TextCtrl(pane, -1, "Your name here")
        #textCtrl.SetSizerProps(expand=True)


        CountsEnablePane = sc.SizedPanel(pane, -1)
        CountsEnablePane.SetSizerType("horizontal")
        CountsEnablePane.SetSizerProps(expand=True)
        # row 1
        wx.StaticText(CountsEnablePane, -1, "Enabled Cross Sections")
        wx.CheckBox(CountsEnablePane, -1, "off off")
        wx.CheckBox(CountsEnablePane, -1, "off on")
        wx.CheckBox(CountsEnablePane, -1, "on off")
        wx.CheckBox(CountsEnablePane, -1, "off off")


        CountsAdd1Pane = sc.SizedPanel(pane, -1)
        CountsAdd1Pane.SetSizerType("horizontal")
        CountsAdd1Pane.SetSizerProps(expand=True)
        # row 1
        wx.StaticText(CountsAdd1Pane, -1, "Summed Cross Sections")
        #wx.StaticText(CountsAdd1Pane, -1, "C1->C1+C2")
        wx.CheckBox(CountsAdd1Pane, -1, "off off")
        wx.CheckBox(CountsAdd1Pane, -1, "off on")
        wx.CheckBox(CountsAdd1Pane, -1, "on off")
        wx.CheckBox(CountsAdd1Pane, -1, "off off")



        CountsAdd2Pane = sc.SizedPanel(pane, -1)
        CountsAdd2Pane.SetSizerType("horizontal")
        CountsAdd2Pane.SetSizerProps(expand=True)
        # row 1
        wx.StaticText(CountsAdd2Pane, -1, "Summed Cross Sections")
        #wx.StaticText(CountsAdd2Pane, -1, "Ca->Ca+Cb")
        wx.CheckBox(CountsAdd2Pane, -1, "off off")
        wx.CheckBox(CountsAdd2Pane, -1, "off on")
        wx.CheckBox(CountsAdd2Pane, -1, "on off")
        wx.CheckBox(CountsAdd2Pane, -1, "off off")



        ConstraintPane = sc.SizedPanel(pane, -1)
        ConstraintPane.SetSizerType("horizontal")
        ConstraintPane.SetSizerProps(expand=True)
        # row 1
        wx.StaticText(ConstraintPane, -1, "Constrained Cross Sections")
        #wx.StaticText(CountsAdd2Pane, -1, "Ca->Ca+Cb")
        wx.CheckBox(ConstraintPane, -1, "off off")
        wx.CheckBox(ConstraintPane, -1, "off on")
        wx.CheckBox(ConstraintPane, -1, "on off")
        wx.CheckBox(ConstraintPane, -1, "off off")


        mmConstraintPane = sc.SizedPanel(pane, -1)
        mmConstraintPane.SetSizerType("horizontal")
        mmConstraintPane.SetSizerProps(expand=True)
        # row 1
        wx.StaticText(mmConstraintPane, -1, "Constraint Coeff off off")
        #wx.StaticText(CountsAdd2Pane, -1, "Ca->Ca+Cb")
        wx.StaticText(mmConstraintPane, -1, "off off")
        wx.lib.intctrl.IntCtrl(mmConstraintPane, -1)
        wx.StaticText(mmConstraintPane, -1, "off on")
        wx.lib.intctrl.IntCtrl(mmConstraintPane, -1)
        wx.StaticText(mmConstraintPane, -1, "on off")
        wx.lib.intctrl.IntCtrl(mmConstraintPane, -1)
        wx.StaticText(mmConstraintPane, -1, "off off")
        wx.lib.intctrl.IntCtrl(mmConstraintPane, -1)


        mpConstraintPane = sc.SizedPanel(pane, -1)
        mpConstraintPane.SetSizerType("horizontal")
        mpConstraintPane.SetSizerProps(expand=True)
        # row 1
        wx.StaticText(mpConstraintPane, -1, "Constraint Coeff off on")
        #wx.StaticText(CountsAdd2Pane, -1, "Ca->Ca+Cb")
        wx.StaticText(mpConstraintPane, -1, "off off")
        wx.lib.intctrl.IntCtrl(mpConstraintPane, -1)
        wx.StaticText(mpConstraintPane, -1, "off on")
        wx.lib.intctrl.IntCtrl(mpConstraintPane, -1)
        wx.StaticText(mpConstraintPane, -1, "on off")
        wx.lib.intctrl.IntCtrl(mpConstraintPane, -1)
        wx.StaticText(mpConstraintPane, -1, "off off")
        wx.lib.intctrl.IntCtrl(mpConstraintPane, -1)


        pmConstraintPane = sc.SizedPanel(pane, -1)
        pmConstraintPane.SetSizerType("horizontal")
        pmConstraintPane.SetSizerProps(expand=True)
        # row 1
        wx.StaticText(pmConstraintPane, -1, "Constraint Coeff on off")
        #wx.StaticText(CountsAdd2Pane, -1, "Ca->Ca+Cb")
        wx.StaticText(pmConstraintPane, -1, "off off")
        wx.lib.intctrl.IntCtrl(pmConstraintPane, -1)
        wx.StaticText(pmConstraintPane, -1, "off on")
        wx.lib.intctrl.IntCtrl(pmConstraintPane, -1)
        wx.StaticText(pmConstraintPane, -1, "on off")
        wx.lib.intctrl.IntCtrl(pmConstraintPane, -1)
        wx.StaticText(pmConstraintPane, -1, "off off")
        wx.lib.intctrl.IntCtrl(pmConstraintPane, -1)


        ppConstraintPane = sc.SizedPanel(pane, -1)
        ppConstraintPane.SetSizerType("horizontal")
        ppConstraintPane.SetSizerProps(expand=True)
        # row 1
        wx.StaticText(ppConstraintPane, -1, "Constraint Coeff on on")
        #wx.StaticText(CountsAdd2Pane, -1, "Ca->Ca+Cb")
        wx.StaticText(ppConstraintPane, -1, "off off")
        wx.lib.intctrl.IntCtrl(ppConstraintPane, -1)
        wx.StaticText(ppConstraintPane, -1, "off on")
        wx.lib.intctrl.IntCtrl(ppConstraintPane, -1)
        wx.StaticText(ppConstraintPane, -1, "on off")
        wx.lib.intctrl.IntCtrl(ppConstraintPane, -1)
        wx.StaticText(ppConstraintPane, -1, "off off")
        wx.lib.intctrl.IntCtrl(ppConstraintPane, -1)



        # add dialog buttons
        self.SetButtonSizer(self.CreateStdDialogButtonSizer(wx.OK | wx.CANCEL))

        # a little trick to make sure that you can't resize the dialog to
        # less screen space than the controls need
        self.Fit()
        self.SetMinSize(self.GetSize())






if __name__=='__main__':
    app=MyApp()
    dlg=FormDialog(parent=None,id=-1)
    result=dlg.ShowModal()
    if result==wx.ID_OK:
        print "OK"
    else:
        print "Cancel"
    dlg.Destroy()
