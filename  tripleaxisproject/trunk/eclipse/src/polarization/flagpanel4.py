import wx
import wxaddons.sized_controls as sc
import  wx.lib.intctrl
import  wx.grid as  gridlib
import numpy as N
import sys,os

class MyApp(wx.App):
    def __init__(self, redirect=False, filename=None, useBestVisual=False, clearSigInt=True):
        wx.App.__init__(self,redirect,filename,clearSigInt)


    def OnInit(self):
        return True


class ConstraintMatrixTable(gridlib.PyGridTableBase):
    def __init__(self,pbflags):
        gridlib.PyGridTableBase.__init__(self)
        #pp mm pm mp
        self.colLabels = ['selected?','on on', 'off off', 'on off','off on']
        self.rowLabels=['on on', 'off off', 'on off','off on']

        self.dataTypes = [gridlib.GRID_VALUE_STRING,# selected?
                          gridlib.GRID_VALUE_STRING, #off off
                          gridlib.GRID_VALUE_STRING,#off on
                          gridlib.GRID_VALUE_STRING,#on off
                          gridlib.GRID_VALUE_STRING,#on on
                          #gridlib.GRID_VALUE_STRING,#row labels
                          ]
        #data = []
        data=N.zeros((4,5),'Float64')
        data=data.tolist()
        for col in range(1,5,1):
            #print 'col setup',col
            data[0][col]=str(pbflags.Spp[col-1]) #row,col
            data[1][col]=str(pbflags.Smm[col-1]) #row,col
            data[2][col]=str(pbflags.Spm[col-1]) #row,col
            data[3][col]=str(pbflags.Smp[col-1]) #row,col
        for row in range(4):
            data[row][row+1]=''
            data[row][0]=''

        #data[0][5]='off off'
        #data[1][5]='off on'
        #data[2][5]='on off'
        #data[3][5]='on on'
        self.data=data

        return

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

    def GetRowValues(self,row):
        rowvalues=[]
        for col in range(self.GetNumberCols()):
            try:
                rowvalues.append(float(self.GetValue(row,col)))
            except ValueError:
                rowvalues.append(float(0))
        return rowvalues

    def GetColValues(self,col):
        colvalues=[]
        for row in range(self.GetNumberRows()):
            if col >0:
                try:
                    colvalues.append(float(self.GetValue(row,col)))
                except ValueError:
                    colvalues.append(int(0))
            else:
                selected=self.GetValue(row,col)
                if selected=='x':
                    colvalues.append(int(1))
                else:
                    colvalues.append(int(0))
        return colvalues



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
    def __init__(self, parent,pbflags):
        gridlib.Grid.__init__(self, parent, -1)
        table = ConstraintMatrixTable(pbflags)
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
        attr.SetReadOnly(True)
        self.SetColAttr(0,attr)
        for col in range(1,5):
            attr=gridlib.GridCellAttr()
            attr.SetReadOnly(False)
            #attr.SetBackgroundColour('grey' if col%2 else (139, 139, 122))
            #attr.SetTextColour((167,167,122) if col%2 else (139, 139, 122))
            self.SetColAttr(col,attr)
        #keep diagonal blank
        for row in range(table.GetNumberRows()):
            self.SetReadOnly(row,row+1,True)
        #attr.SetReadOnly(True)
        #self.SetColAttr(table.GetNumberCols()-1,attr)
        for row in range(table.GetNumberRows()):
                         if pbflags.Sconstrain[row]==1:
                             table.SetValue(row,0,'x')
        gridlib.EVT_GRID_CELL_LEFT_CLICK(self,self.OnLeftClick)
        gridlib.Grid.ForceRefresh(self)
        #gridlib.EVT_GRID_CELL_LEFT_DCLICK(self, self.OnLeftDClick)
        #gridlib.EVT_GRID_LABEL_LEFT_DCLICK(self,self.onLeftDClickRowCell)

    # I do this because I don't like the default behaviour of not starting the
    # cell editor on double clicks, but only a second click.
    def OnLeftDClick(self, evt):
        if self.CanEnableCellControl():
            self.EnableCellEditControl()

    def onLeftDClickRowCell(self,evt):
        pass
        #col=evt.GetCol()
        #table=self.GetTable()

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
        else:
            evt.Skip()
        #if self.CanEnableCellControl():
        #    self.EnableCellEditControl()
        gridlib.Grid.ForceRefresh(self)


class mFormValidator(wx.PyValidator):
    def __init__(self,data,key):
        wx.PyValidator.__init__(self)
        self.data=data
        self.key=key
        #print 'FormValidator init', key
        #self.TransferToWindow()

    def Clone(self):
        return mFormValidator(self.data,self.key)

    def Validate(self,win):
        return True

    def TransferToWindow(self):
        #print 'Form TransferToWindow',self.key
        
        checkctrl=self.GetWindow()
        ##print 'checkctrl',checkctrl
        ##print self.__dict__
        #print 'key',self.key
        #print 'dict', self.data._array[self.key]
        #print 'data',self.data.__dict__[self.key]
        checkctrl.SetValue(self.data._array[self.key])
        return True

    def TransferFromWindow(self):
        #print 'TransferFromWindow'
        checkctrl=self.GetWindow()
        self.data._array[self.key]=int(checkctrl.GetValue())
        return True

class ListValidator(wx.PyValidator):
    def __init__(self,data,key,index):
        wx.PyValidator.__init__(self)
        #print 'ListValidator init', key
        self.data=data
        self.key=key
        self.index=index

    def Clone(self):
        return ListValidator(self.data,self.key,self.index)

    def Validate(self,win):
        return True

    def TransferToWindow(self):
        #print 'List Transfer2Window',self.key
        checkctrl=self.GetWindow()
        checkctrl.SetValue(self.data._array[self.key][self.index])
        return True

    def TransferFromWindow(self):
        checkctrl=self.GetWindow()
        self.data._array[self.key][self.index]=int(checkctrl.GetValue())
        return True



class FormDialog(sc.SizedDialog):
    def __init__(self, parent, id,individualdata=None,groupdata=None):
        sc.SizedDialog.__init__(self, None, -1, "Reduction Forms",
                        style=wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER)#| wx.WS_EX_VALIDATE_RECURSIVELY)

        pane = self.GetContentsPane()
        pane.SetSizerType("vertical")
        self.groupdata=groupdata
        self.cellfile=groupdata['cellfile']
        FilePane = sc.SizedPanel(pane, -1)
        FilePane.SetSizerType("vertical")
        FilePane.SetSizerProps(expand=True)
        self.cellfilectrl=wx.StaticText(FilePane, -1, "CellFile:%s"%(self.cellfile,))
        b = wx.Button(FilePane, -1, "Browse", (50,50))
        self.Bind(wx.EVT_BUTTON, self.OnOpen, b)


        FilterPane = sc.SizedPanel(pane, -1)
        FilterPane.SetSizerType("vertical")
        FilterPane.SetSizerProps(expand=True)
        filter=wx.CheckBox(FilterPane, -1, "Filter Before Polarizer?",validator=mFormValidator(groupdata['pbflags'],'MonoSelect'))
        self.Bind(wx.EVT_CHECKBOX, self.EvtFilter, filter)
        print 'Filtered'


        monitorPane = sc.SizedPanel(pane, -1)
        monitorPane.SetSizerType("horizontal")
        monitorPane.SetSizerProps(expand=True)
        # row 1
        wx.StaticText(monitorPane, -1, "Monitor Position")
        prepolarizer=wx.CheckBox(monitorPane, -1, "PrePolarizer",validator=mFormValidator(groupdata['pbflags'],'MonitorCorrect'))
        postpolarizer=wx.CheckBox(monitorPane, -1, "PostPolarizer",validator=mFormValidator(groupdata['pbflags'],'PolMonitorCorrect'))
        self.Bind(wx.EVT_CHECKBOX, self.EvtPrePolarizer, prepolarizer)
        self.Bind(wx.EVT_CHECKBOX, self.EvtPostPolarizer, postpolarizer)
        print 'Monitored'

        CountsEnablePane = sc.SizedPanel(pane, -1)
        CountsEnablePane.SetSizerType("horizontal")
        CountsEnablePane.SetSizerProps(expand=True)
        # row 1
        #pp mm pm mp
        wx.StaticText(CountsEnablePane, -1, "Enabled Measured Counts")
        ce1=wx.CheckBox(CountsEnablePane, -1, "off off",validator=ListValidator(groupdata['pbflags'],'CountsEnable',1))
        ce2=wx.CheckBox(CountsEnablePane, -1, "off on",validator=ListValidator(groupdata['pbflags'],'CountsEnable',3))
        ce3=wx.CheckBox(CountsEnablePane, -1, "on off",validator=ListValidator(groupdata['pbflags'],'CountsEnable',2))
        ce4=wx.CheckBox(CountsEnablePane, -1, "on on",validator=ListValidator(groupdata['pbflags'],'CountsEnable',0))
        self.Bind(wx.EVT_CHECKBOX, self.EvtCountsEnable_offoff, ce1)
        self.Bind(wx.EVT_CHECKBOX, self.EvtCountsEnable_offon, ce2)
        self.Bind(wx.EVT_CHECKBOX, self.EvtCountsEnable_onoff, ce3)
        self.Bind(wx.EVT_CHECKBOX, self.EvtCountsEnable_onon, ce4)
        print 'CountsEnabled'
        
        FilterPane.TransferDataToWindow()
        monitorPane.TransferDataToWindow()
        CountsEnablePane.TransferDataToWindow()
        #pane.TransferDataToWindow()


        CountsAddPane = sc.SizedPanel(pane, -1)
        CountsAddPane.SetSizerType("horizontal")
        CountsAddPane.SetSizerProps(expand=True)
        # row 1
        wx.StaticText(CountsAddPane, -1, "Combine Measured Counts")

        #wx.StaticText(CountsAdd1Pane, -1, "C1->C1+C2")
        nsf=wx.CheckBox(CountsAddPane, -1, "NSF")
        sf=wx.CheckBox(CountsAddPane, -1, "SF")

        self.Bind(wx.EVT_CHECKBOX, self.EvtNSF, nsf)
        self.Bind(wx.EVT_CHECKBOX, self.EvtSF, sf)


        ConstraintPane = sc.SizedPanel(pane, -1)
        ConstraintPane.SetSizerType("vertical")
        ConstraintPane.SetSizerProps(expand=True)
        # row 1
        wx.StaticText(ConstraintPane, -1, "Constraint Matrix")
        grid = ConstraintMatrixGrid(ConstraintPane,groupdata['pbflags'])
        self.grid=grid

        # add dialog buttons
        self.SetButtonSizer(self.CreateStdDialogButtonSizer(wx.OK | wx.CANCEL))
        self.TransferDataToWindow()
        # a little trick to make sure that you can't resize the dialog to
        # less screen space than the controls need
        self.Fit()
        self.SetMinSize(self.GetSize())


    def EvtFilter(self,evt):
        self.groupdata['pbflags'].MonoSelect=int(evt.IsChecked())
    def EvtPrePolarizer(self,evt):
        self.groupdata['pbflags'].MonitorCorrect=int(evt.IsChecked())

    def EvtPostPolarizer(self,evt):
        self.groupdata['pbflags'].PolMonitorCorrect=int(evt.IsChecked())

    def EvtCountsEnable_onon(self,evt):
        self.groupdata['pbflags'].CountsEnable[0]=int(evt.IsChecked())
    def EvtCountsEnable_offoff(self,evt):
        self.groupdata['pbflags'].CountsEnable[1]=int(evt.IsChecked())
    def EvtCountsEnable_onoff(self,evt):
        self.groupdata['pbflags'].CountsEnable[2]=int(evt.IsChecked())
    def EvtCountsEnable_offon(self,evt):
        self.groupdata['pbflags'].CountsEnable[3]=int(evt.IsChecked())

#pp mm pm mp
    def EvtNSF(self,evt):
        if evt.IsChecked():
            self.groupdata['pbflags'].CountsAdd1[0]=1
            self.groupdata['pbflags'].CountsAdd1[1]=2
            #self.groupdata['pbflags'].CountsEnable[1]=0
        else:
            self.groupdata['pbflags'].CountsAdd1[0]=0
            self.groupdata['pbflags'].CountsAdd1[1]=0

    def EvtSF(self,evt):
        if evt.IsChecked():
            self.groupdata['pbflags'].CountsAdd2[2]=3
            self.groupdata['pbflags'].CountsAdd2[3]=4
            #self.groupdata['pbflags'].CountsEnable[3]=0
        else:
            self.groupdata['pbflags'].CountsAdd1[2]=0
            self.groupdata['pbflags'].CountsAdd1[3]=0





    def OnOpen(self,event):
        # Create the dialog. In this case the current directory is forced as the starting
        # directory for the dialog, and no default file name is forced. This can easilly
        # be changed in your program. This is an 'open' dialog, and allows multitple
        # file selections as well.
        #
        # Finally, if the directory is changed in the process of getting files, this
        # dialog is set up to change the current working directory to the path chosen.

        #defaultDir=os.getcwd()
        #defaultDir=r'C:\polcorrecter\data'
        defaultDir=wx.ConfigBase.Get().GetPath()
        wildcard="cell files (*.txt)|*.txt|All files (*.*)|*.*"
        dlg = wx.FileDialog(
            self, message="Choose a file",
            defaultDir=defaultDir,
            defaultFile="",
            wildcard=wildcard,
            style=wx.OPEN | wx.CHANGE_DIR
            )

        # Show the dialog and retrieve the user response. If it is the OK response,
        # process the data.
        if dlg.ShowModal() == wx.ID_OK:
            # This returns a Python list of files that were selected.
            paths = dlg.GetPaths()
            #self.log.WriteText('You selected %d files:' % len(paths))
            self.groupdata['cellfile']=paths[0].encode('ascii')
            self.cellfilectrl.SetLabel("CellFile:%s"%(self.groupdata['cellfile'],))
            #wx.StaticText(FilePane, -1, "CellFile:%s"%(self.groupdata['cellfile'],))
        # Destroy the dialog. Don't do this until you are done with it!
        # BAD things can happen otherwise!
        dlg.Destroy()




if __name__=='__main__':
    app=MyApp()
    dlg=FormDialog(parent=None,id=-1)
    result=dlg.ShowModal()
    if result==wx.ID_OK:
        print "OK"
    else:
        print "Cancel"
    dlg.Destroy()
