import wx,wx.gizmos
import string
import os,sys
import wx.lib.colourselect as csel
import wx.lib.customtreectrl as CT
import flagpanel4 as flagpanel
import polcorrect3 as polcorrect
import numpy as N
from sans.guicomm.events import NewPlotEvent
from sans.guitools.plottables import Data1D
import copy
from myevents import *

#import images
#try:
#    import treemixin
#except ImportError:
#    from wx.lib.mixins import treemixin





penstyle = ["wx.SOLID", "wx.TRANSPARENT", "wx.DOT", "wx.LONG_DASH", "wx.DOT_DASH", "wx.USER_DASH",
           "wx.BDIAGONAL_HATCH", "wx.CROSSDIAG_HATCH", "wx.FDIAGONAL_HATCH", "wx.CROSS_HATCH",
           "wx.HORIZONTAL_HATCH", "wx.VERTICAL_HATCH"]

ArtIDs = [ "None",
           "wx.ART_ADD_BOOKMARK",
           "wx.ART_DEL_BOOKMARK",
           "wx.ART_HELP_SIDE_PANEL",
           "wx.ART_HELP_SETTINGS",
           "wx.ART_HELP_BOOK",
           "wx.ART_HELP_FOLDER",
           "wx.ART_HELP_PAGE",
           "wx.ART_GO_BACK",
           "wx.ART_GO_FORWARD",
           "wx.ART_GO_UP",
           "wx.ART_GO_DOWN",
           "wx.ART_GO_TO_PARENT",
           "wx.ART_GO_HOME",
           "wx.ART_FILE_OPEN",
           "wx.ART_PRINT",
           "wx.ART_HELP",
           "wx.ART_TIP",
           "wx.ART_REPORT_VIEW",
           "wx.ART_LIST_VIEW",
           "wx.ART_NEW_DIR",
           "wx.ART_HARDDISK",
           "wx.ART_FLOPPY",
           "wx.ART_CDROM",
           "wx.ART_REMOVABLE",
           "wx.ART_FOLDER",
           "wx.ART_FOLDER_OPEN",
           "wx.ART_GO_DIR_UP",
           "wx.ART_EXECUTABLE_FILE",
           "wx.ART_NORMAL_FILE",
           "wx.ART_TICK_MARK",
           "wx.ART_CROSS_MARK",
           "wx.ART_ERROR",
           "wx.ART_QUESTION",
           "wx.ART_WARNING",
           "wx.ART_INFORMATION",
           "wx.ART_MISSING_IMAGE",
           "SmileBitmap"
           ]

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
# CustomTreeCtrl Demo Implementation
#---------------------------------------------------------------------------
class CustomTreeCtrl(CT.CustomTreeCtrl):

    def __init__(self, parent, id=wx.ID_ANY, pos=wx.DefaultPosition,
                 size=wx.DefaultSize,
                 style=wx.SUNKEN_BORDER | CT.TR_HAS_BUTTONS | CT.TR_HAS_VARIABLE_ROW_HEIGHT|wx.WS_EX_VALIDATE_RECURSIVELY,
                 log=None):

        CT.CustomTreeCtrl.__init__(self, parent, id, pos, size, style)
        self.parent=parent
        alldata = dir(CT)

        treestyles = []
        events = []
        for data in alldata:
            if data.startswith("TR_"):
                treestyles.append(data)
            elif data.startswith("EVT_"):
                events.append(data)

        self.events = events
        self.styles = treestyles
        self.item = None

        il = wx.ImageList(16, 16)

        for items in ArtIDs[1:-1]:
            bmp = wx.ArtProvider_GetBitmap(eval(items), wx.ART_TOOLBAR, (16, 16))
            il.Add(bmp)

        #smileidx = il.Add(images.getSmilesBitmap())
        numicons = il.GetImageCount()

        self.AssignImageList(il)
        self.count = 0
        self.log = log

        # NOTE:  For some reason tree items have to have a data object in
        #        order to be sorted.  Since our compare just uses the labels
        #        we don't need any real data, so we'll just use None below for
        #        the item data.

        self.root = self.AddRoot("The Root Item")

        if not(self.GetTreeStyle() & CT.TR_HIDE_ROOT):
            self.SetPyData(self.root, None)
            self.SetItemImage(self.root, 24, CT.TreeItemIcon_Normal)
            self.SetItemImage(self.root, 13, CT.TreeItemIcon_Expanded)



        child = self.AppendItem(self.root, "Raw Data")
        self.SetItemBold(child, True)
        self.SetPyData(child, None)
        self.SetItemImage(child, 24, CT.TreeItemIcon_Normal)
        self.SetItemImage(child, 13, CT.TreeItemIcon_Expanded)

        child = self.AppendItem(self.root, "Reduced Data")
        self.SetItemBold(child, True)
        self.SetPyData(child, None)
        self.SetItemImage(child, 24, CT.TreeItemIcon_Normal)
        self.SetItemImage(child, 13, CT.TreeItemIcon_Expanded)
        #self.gauge = wx.Gauge(self, -1, 50, style=wx.GA_HORIZONTAL|wx.GA_SMOOTH)
        #self.gauge.SetValue(0)

        self.Bind(wx.EVT_LEFT_DCLICK, self.OnLeftDClick)
        #self.Bind(wx.EVT_IDLE, self.OnIdle)
        self.eventdict = {'EVT_TREE_BEGIN_DRAG': self.OnBeginDrag, 'EVT_TREE_BEGIN_LABEL_EDIT': self.OnBeginEdit,
                          'EVT_TREE_BEGIN_RDRAG': self.OnBeginRDrag, 'EVT_TREE_DELETE_ITEM': self.OnDeleteItem,
                          'EVT_TREE_END_DRAG': self.OnEndDrag, 'EVT_TREE_END_LABEL_EDIT': self.OnEndEdit,
                          'EVT_TREE_ITEM_ACTIVATED': self.OnActivate, 'EVT_TREE_ITEM_CHECKED': self.OnItemCheck,
                          'EVT_TREE_ITEM_CHECKING': self.OnItemChecking, 'EVT_TREE_ITEM_COLLAPSED': self.OnItemCollapsed,
                          'EVT_TREE_ITEM_COLLAPSING': self.OnItemCollapsing, 'EVT_TREE_ITEM_EXPANDED': self.OnItemExpanded,
                          'EVT_TREE_ITEM_EXPANDING': self.OnItemExpanding, 'EVT_TREE_ITEM_GETTOOLTIP': self.OnToolTip,
                          'EVT_TREE_ITEM_MENU': self.OnItemMenu, 'EVT_TREE_ITEM_RIGHT_CLICK': self.OnRightDown,
                          'EVT_TREE_KEY_DOWN': self.OnKey, 'EVT_TREE_SEL_CHANGED': self.OnSelChanged,
                          'EVT_TREE_SEL_CHANGING': self.OnSelChanging, "EVT_TREE_ITEM_HYPERLINK": self.OnHyperLink}

        #mainframe = wx.GetTopLevelParent(self)

        #if not hasattr(mainframe, "leftpanel"):
        self.Bind(CT.EVT_TREE_ITEM_EXPANDED, self.OnItemExpanded)
        self.Bind(CT.EVT_TREE_ITEM_COLLAPSED, self.OnItemCollapsed)
        self.Bind(CT.EVT_TREE_SEL_CHANGED, self.OnSelChanged)
        self.Bind(CT.EVT_TREE_SEL_CHANGING, self.OnSelChanging)
        self.Bind(wx.EVT_RIGHT_DOWN, self.OnRightDown)
        self.Bind(wx.EVT_RIGHT_UP, self.OnRightUp)
        self.Bind(EVT_CLEAR_TREE,self.OnClearTree)
        #self.Bind(wx.EVT_LEFT_DCLICK, self.OnRightUp)
        #else:
        #    for combos in mainframe.treeevents:
        #        self.BindEvents(combos)

        #if hasattr(mainframe, "leftpanel"):
        #    self.ChangeStyle(mainframe.treestyles)

        if not(self.GetTreeStyle() & CT.TR_HIDE_ROOT):
            self.SelectItem(self.root)
            self.Expand(self.root)


    def BindEvents(self, choice, recreate=False):

        value = choice.GetValue()
        text = choice.GetLabel()

        evt = "CT." + text
        binder = self.eventdict[text]

        if value == 1:
            if evt == "CT.EVT_TREE_BEGIN_RDRAG":
                self.Bind(wx.EVT_RIGHT_DOWN, None)
                self.Bind(wx.EVT_RIGHT_UP, None)
            self.Bind(eval(evt), binder)
        else:
            self.Bind(eval(evt), None)
            if evt == "CT.EVT_TREE_BEGIN_RDRAG":
                self.Bind(wx.EVT_RIGHT_DOWN, self.OnRightDown)
                self.Bind(wx.EVT_RIGHT_UP, self.OnRightUp)


    def ChangeStyle(self, combos):

        style = 0
        for combo in combos:
            if combo.GetValue() == 1:
                style = style | eval("CT." + combo.GetLabel())

        if self.GetTreeStyle() != style:
            self.SetTreeStyle(style)


    def OnCompareItems(self, item1, item2):

        t1 = self.GetItemText(item1)
        t2 = self.GetItemText(item2)

        self.log.write('compare: ' + t1 + ' <> ' + t2 + "\n")

        if t1 < t2:
            return -1
        if t1 == t2:
            return 0

        return 1

    def OnClearTree(self, event):
        root=self.GetRootItem()
        if root!=None:
            if self.HasChildren(root):
                raw_data,cookie=self.GetFirstChild(root)
                if raw_data!=None:
                    self.DeleteChildren(raw_data)
                    reduced_data,cookie=self.GetNextChild(root,cookie)
                    self.DeleteChildren(reduced_data)
                    
        
        
        print 'Clear Event caught'
        return
    
    def OnIdle(self, event):

        if self.gauge:
            try:
                if self.gauge.IsEnabled() and self.gauge.IsShown():
                    self.count = self.count + 1

                    if self.count >= 50:
                        self.count = 0

                    self.gauge.SetValue(self.count)

            except:
                self.gauge = None

        event.Skip()


    def OnRightDown(self, event):

        pt = event.GetPosition()
        item, flags = self.HitTest(pt)

        if item:
            self.item = item
            self.log.write("OnRightClick: %s, %s, %s" % (self.GetItemText(item), type(item), item.__class__) + "\n")
            self.SelectItem(item)


    def OnRightUp(self, event):

        item = self.item

        if not item:
            event.Skip()
            return



        #if not self.IsEnabled(item):
        #    event.Skip()
        #    return

        # Item Text Appearance
        ishtml = self.IsItemHyperText(item)
        back = self.GetItemBackgroundColour(item)
        fore = self.GetItemTextColour(item)
        isbold = self.IsBold(item)
        font = self.GetItemFont(item)

        # Icons On Item
        normal = self.GetItemImage(item, CT.TreeItemIcon_Normal)
        selected = self.GetItemImage(item, CT.TreeItemIcon_Selected)
        expanded = self.GetItemImage(item, CT.TreeItemIcon_Expanded)
        #selexp = self.GetItemImage(item, CT.TreeItemIcon_SelectedExpanded)


        # Generic Item's Info
        children = self.GetChildrenCount(item)
        itemtype = self.GetItemType(item)
        text = self.GetItemText(item)
        pydata = self.GetPyData(item)

        self.current = item
        self.itemdict = {"ishtml": ishtml, "back": back, "fore": fore, "isbold": isbold,
                         "font": font, "normal": normal, "selected": selected, "expanded": expanded,
                         "itemtype": itemtype, "text": text, "pydata": pydata}


        DataGroup=self.GetFirstChild(self.GetRootItem())[0]
        if self.GetItemParent(item)==DataGroup:
            menu = wx.Menu()
            item4 = menu.Append(wx.ID_ANY, "Load Reduction Preferences")
            menu.AppendSeparator()
            item1 = menu.Append(wx.ID_ANY, "Set Reduction Preferences")
            menu.AppendSeparator()
            item2 = menu.Append(wx.ID_ANY, "Reduce Group")
            menu.AppendSeparator()
            item3 = menu.Append(wx.ID_ANY, "Plot Group")
            

            #item10 = menu.Append(wx.ID_ANY, "Delete Item")
            #if item == self.GetRootItem():
            #    item10.Enable(False)
            #if item in self.GetFirstChild(self.GetRootItem()):
            #    item10.Enable(False)
            self.Bind(wx.EVT_MENU, self.OnItemSetReductionPreferences, item1)
            self.Bind(wx.EVT_MENU, self.OnItemReduceGroup, item2)
            self.Bind(wx.EVT_MENU, self.OnItemPlot, item3)
            self.Bind(wx.EVT_MENU, self.OnItemLoad, item4)
            #self.Bind(wx.EVT_MENU, self.OnItemForeground, item2)
            #self.Bind(wx.EVT_MENU, self.OnItemBold, item3)
            #self.Bind(wx.EVT_MENU, self.OnItemFont, item4)
            #self.Bind(wx.EVT_MENU, self.OnItemHyperText, item5)
            #self.Bind(wx.EVT_MENU, self.OnEnableWindow, item6)
            #self.Bind(wx.EVT_MENU, self.OnDisableItem, item7)
            #self.Bind(wx.EVT_MENU, self.OnItemIcons, item8)
            #self.Bind(wx.EVT_MENU, self.OnItemInfo, item9)
            #self.Bind(wx.EVT_MENU, self.OnItemDelete, item10)
            #self.Bind(wx.EVT_MENU, self.OnItemPrepend, item11)
            #self.Bind(wx.EVT_MENU, self.OnItemAppend, item12)

            self.PopupMenu(menu)
            menu.Destroy()

    def OnItemPlot(self,event):
        current_selected=self.current
        pydata=self.itemdict['pydata']
        children_data=self.GetChildData(self.current)
        i=0
        for data in children_data:
            #print data
            self.varying=data['data'].metadata['count_info']['varying']
            self.independent_variable=self.varying[0]
            x=data['data'].data[self.independent_variable]
            self.x=x
            y=data['data'].data['detector']
            dy=N.sqrt(y)
            print 'varying',self.independent_variable

            #from sans.guicomm.events import NewPlotEvent
            #from sans.guitools.plottables import Data1D
            new_plot = Data1D(x, y, dy=dy)
            new_plot.name =data['filename']+' '+data['polstate']
            new_plot.group_id='data'
            new_plot.xaxis(str(self.independent_variable), 'A^{-1}')
            new_plot.yaxis("\\rm{Intensity} ","Arb. units")
            wx.PostEvent(self.parent.parent, NewPlotEvent(plot=new_plot))

    def OnItemLoad(self,event):
        defaultDir=os.getcwd()
        defaultDir=r'C:\polcorrecter\data'
        wildcard="driver files (*.polcor)|*.polcor|All files (*.*)|*.*"
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
            myfilestr=paths[0].encode('ascii')
            myfile = open(myfilestr, 'r')
            mycount=0
            returnline=['']
            catalog=[]
            cellfiles=[]
            absolute=True
            inblock=False
            groupdata=self.itemdict['pydata']
            pbflags=groupdata['pbflags']
            pbflags.MonitorCorrect=0
            pbflags.PolMonitorCorrect=1
            pbflags.MonoSelect=1
            pbflags.Debug=0
            pbflags.SimFlux=0
            pbflags.SimDeviate=0
            pbflags.NoNegativeCS=0
            pbflags.HalfPolarized=0
            pbflags.CountsEnable=[0,0,0,0]
            pbflags.CountsAdd1=[0,0,0,0]
            pbflags.CountsAdd2=[0,0,0,0]
            pbflags.Sconstrain=[0,0,0,0]
            pbflags.Spp=[0,0,0,0]
            pbflags.Smm=[0,0,0,0]
            pbflags.Spm=[0,0,0,0]
            pbflags.Smp=[0,0,0,0]
        
            while 1:
                lineStr = myfile.readline()
                if not(lineStr):
                    break
                strippedLine=lineStr.rstrip().lower()
                tokenized=strippedLine.split()
                #print 'tokenized ', tokenized
                directive=tokenized[0]
                if tokenized[0]==[]:
                    pass
                elif tokenized[0]=='#absolute'.lower():
                    #print 'absolute'
                    absolute=True
                    mydirectory=os.curdir
                elif tokenized[0]=='#directory'.lower():
                    #print 'directory'
                    if os.path.isdir(mydirectory):
                        mydirectory=tokenized[1]
                    else:
                        print '%s is not an existing directory!!!'%(tokenized[1],)
                        break
                elif tokenized[0]=='#MonitorCorrect'.lower():
                    pbflags.MonitorCorrect=int(tokenized[1])
                    #print '##Monitor Correct ',tokenized
                elif tokenized[0]=='#PolMonitorCorrect'.lower():
                    pbflags.PolMonitorCorrect=int(tokenized[1])
                elif tokenized[0]=='#MonoSelect'.lower:
                    pbflags.MonoSelect=int(tokenized[1])
                elif tokenized[0]=='#CountsEnable'.lower():
                    if len(tokenized)!=5:
                        print 'Too Few Arguments for the %s pragma.  An example of proper use is:'%(directive,)
                        print '%s 1 0 1 3'%(directive,)
                        break
                    else:
                        pbflags.CountsEnable[0]=int(tokenized[1])
                        pbflags.CountsEnable[1]=int(tokenized[2])
                        pbflags.CountsEnable[2]=int(tokenized[3])
                        pbflags.CountsEnable[3]=int(tokenized[4])
                elif tokenized[0]=='#CountsAdd1'.lower():
                    if len(tokenized)!=5:
                        print 'Too Few Arguments for the %s pragma.  An example of proper use is:'%(directive,)
                        print '%s 1 0 1 3'%(directive,)
                        break
                    else:
                        pbflags.CountsAdd2[0]=int(tokenized[1])
                        pbflags.CountsAdd2[1]=int(tokenized[2])
                        pbflags.CountsAdd2[2]=int(tokenized[3])
                        pbflags.CountsAdd2[3]=int(tokenized[4])
                elif tokenized[0]=='#CountsAdd2'.lower():
                    if len(tokenized)!=5:
                        print 'Too Few Arguments for the %s pragma.  An example of proper use is:'%(directive,)
                        print '%s 1 0 1 3'%(directive,)
                        break
                    else:
                        pbflags.CountsAdd2[0]=int(tokenized[1])
                        pbflags.CountsAdd2[1]=int(tokenized[2])
                        pbflags.CountsAdd2[2]=int(tokenized[3])
                        pbflags.CountsAdd2[3]=int(tokenized[4])
                elif tokenized[0]=='#Sconstrain'.lower():
                    #print 'Sconstrain ',tokenized
                    if len(tokenized)!=5:
                        print 'Too Few Arguments for the %s pragma.  An example of proper use is:'%(directive,)
                        print '%s 1 0 1 3'%(directive,)
                        break
                    else:
                        pbflags.Sconstrain[0]=int(tokenized[1])
                        pbflags.Sconstrain[1]=int(tokenized[2])
                        pbflags.Sconstrain[2]=int(tokenized[3])
                        pbflags.Sconstrain[3]=int(tokenized[4])
                elif tokenized[0]=='#Spp'.lower():
                    if len(tokenized)!=5:
                        print 'Too Few Arguments for the %s pragma.  An example of proper use is:'%(directive,)
                        print '%s 1 0 1 3'%(directive,)
                        break
                    else:
                        pbflags.Spp[0]=float(tokenized[1])
                        pbflags.Spp[1]=float(tokenized[2])
                        pbflags.Spp[2]=float(tokenized[3])
                        pbflags.Spp[3]=float(tokenized[4])
                elif tokenized[0]=='#Smm'.lower:
                    if len(tokenized)!=5:
                        print 'Too Few Arguments for the %s pragma.  An example of proper use is:'%(directive,)
                        print '%s 1 0 1 3'%(directive,)
                        break
                    else:
                        pbflags.Smm[0]=float(tokenized[1])
                        pbflags.Smm[1]=float(tokenized[2])
                        pbflags.Smm[2]=float(tokenized[3])
                        pbflags.Smm[3]=float(tokenized[4])
                elif tokenized[0]=='#Spm'.lower():
                    if len(tokenized)!=5:
                        print 'Too Few Arguments for the %s pragma.  An example of proper use is:'%(directive,)
                        print '%s 1 0 1 3'%(directive,)
                        break
                    else:
                        pbflags.Spm[0]=float(tokenized[1])
                        pbflags.Spm[1]=float(tokenized[2])
                        pbflags.Spm[2]=float(tokenized[3])
                        pbflags.Spm[3]=float(tokenized[4])
                elif tokenized[0]=='#Smp'.lower():
                    if len(tokenized)!=5:
                        print 'Too Few Arguments for the %s pragma.  An example of proper use is:'%(directive,)
                        print '%s 1 0 1 3'%(directive,)
                        break
                    else:
                        pbflags.Smp[0]=float(tokenized[1])
                        pbflags.Smp[1]=float(tokenized[2])
                        pbflags.Smp[2]=float(tokenized[3])
                        pbflags.Smp[3]=float(tokenized[4])
        
#                elif tokenized[0]=='#begin'.lower():
#                    inblock=True
#                    #print 'begin '
#                    files={}
#                elif tokenized[0]=='#end'.lower():
#                    #print 'end'
#                    catalog.append(files)
#                    print 'correcting %s using cellfile %s'%(files,cellfile)
#                    mypolcor=polarization_correct(files,cellfile)
#                    corrected_counts=mypolcor.correct(pbflags)
#                    mypolcor.savefiles()
#                    print 'corrected'
#                    if 0:
#                        key='pm'.lower()
#                        #pylab.subplot(2,2,1+mycount)
#                        #pylab.title(key)
#                        mydatac=mypolcor.mydata
#                        #pylab.errorbar(mydatac[key].data['qx'],mydatac[key].data['detector'],N.sqrt(mydatac[key].data['detector']),
#                        #    marker='s',mfc='blue',linestyle='None')
#                        #pylab.errorbar(mydatac[key].data['qx'],corrected_counts['Spm'],corrected_counts['Epm'], marker='s',mfc='red',linestyle='None')
#                        #print 'pm'
#                        #print corrected_counts['Spm']
#                        key='mp'.lower()
#                        #pylab.subplot(2,2,2+mycount)
#                        #pylab.title(key)
#                        #print 'mp'
#                        #print corrected_counts['Smp']
#                        #pylab.errorbar(mydatac[key].data['qx'],mydatac[key].data['detector'],N.sqrt(mydatac[key].data['detector']),
#                        #    marker='s',mfc='blue',linestyle='None')
#                        #pylab.errorbar(mydatac[key].data['qx'],corrected_counts['Smp'],corrected_counts['Emp'], marker='s',mfc='red',linestyle='None')
#                        mycount=mycount+2
#                    inblock=False
#                elif inblock==True:
#                    #print 'inblock'
#                    toksplit=tokenized[0].split('=')
#                    #check to make sure that there actually is a file specified!
#                    if len(toksplit)==2:
#                        filetok=os.path.join(mydirectory,toksplit[1])
#                        if os.path.isfile(filetok):
#                            files[toksplit[0]]=filetok
#                        else:
#                            print filetok+' does not exist!!!'
#                            break
#                    else:
#                        pass
                elif tokenized[0]=='#cell'.lower():
                    #print 'acellfile'
                    cellfile=os.path.join(mydirectory,tokenized[1])
                    groupdata['cellfile']=cellfile
                    if os.path.isfile(cellfile):
                        cellfiles.append(cellfile)
                    else:
                        print '%s does not exist!'%(cellfile,)
                        break
                #elif tokenized[0]=='#Pcell'.lower():
                #    #print 'pcellfile'
                #    pcellfile=os.path.join(mydirectory,tokenized[1])
                #    if os.path.isfile(pcellfile):
                #        pcellfiles.append(pcellfile)
                #    else:
                #        print '%s does not exist!'%(acellfile,)
                #        break
            print 'done'
            myfile.close()
            print 'closed'

        # Destroy the dialog. Don't do this until you are done with it!
        # BAD things can happen otherwise!
        dlg.Destroy()

        return

    def OnItemSetReductionPreferences(self,event):
        current_selected=self.current
        pydata=self.itemdict['pydata']
        children_data=self.GetChildData(self.current)
        files={}
        for currdata in children_data:
            key=currdata['polstate']
            files[key]=currdata['filename']
        self.files=files 
        files=copy.deepcopy(self.files)
        #pp mm pm mp
        if files.has_key('pp'):
            pydata['pbflags'].CountsEnable[0]=1
        if files.has_key('mm'):
            pydata['pbflags'].CountsEnable[1]=1
        if files.has_key('pm'):
            pydata['pbflags'].CountsEnable[2]=1
        if files.has_key('mp'):
            pydata['pbflags'].CountsEnable[3]=1
        
        if pydata['pbflags'].CountsEnable[0]==0:
            pydata['pbflags'].Sconstrain[0]=1
        if pydata['pbflags'].CountsEnable[1]==0:
            pydata['pbflags'].Sconstrain[1]=1
        if pydata['pbflags'].CountsEnable[2]==0:
            pydata['pbflags'].Sconstrain[2]=1
        if pydata['pbflags'].CountsEnable[3]==0:
            pydata['pbflags'].Sconstrain[3]=1
        dlg=flagpanel.FormDialog(parent=self,id=-1,groupdata=pydata,individualdata=children_data)
        #self.TransferDataToWindow()
        result=dlg.ShowModal()
        if result==wx.ID_OK:
            print "OK"
            print 'cellfile',dlg.groupdata['cellfile']
            table=dlg.grid.GetTable()
#        self.colLabels = ['selected?','off off', 'off on', 'on off','on on']
#        self.rowLabels=['off off', 'off on', 'on off','on on']
#pp mm pm mp
            dlg.groupdata['pbflags'].Spp=table.GetRowValues(0)[1:] #row,col
            dlg.groupdata['pbflags'].Smm=table.GetRowValues(1)[1:] #row,col
            dlg.groupdata['pbflags'].Spm=table.GetRowValues(2)[1:] #row,col
            dlg.groupdata['pbflags'].Smp=table.GetRowValues(3)[1:] #row,col
            dlg.groupdata['pbflags'].Sconstrain=table.GetColValues(0) #row,col
            #If combining counts, disable the 2nd cross section
            if dlg.groupdata['pbflags'].CountsAdd1[0]==1:
                dlg.groupdata['pbflags'].CountsEnable[1]=0
            if dlg.groupdata['pbflags'].CountsAdd2[0]==3:
                dlg.groupdata['pbflags'].CountsEnable[3]=0

            print 'MonitorSelect',dlg.groupdata['pbflags'].MonoSelect
            print 'MonitorCorrect',dlg.groupdata['pbflags'].MonitorCorrect
            print 'PolMonitorCorrect',dlg.groupdata['pbflags'].PolMonitorCorrect
            print 'CountsEnable',dlg.groupdata['pbflags'].CountsEnable
            #print 'CountsEnable',dlg.groupdata['pbflags'].CountsEnable
            print 'CountsAdd1',dlg.groupdata['pbflags'].CountsAdd1
            print 'CountsAdd2',dlg.groupdata['pbflags'].CountsAdd2
            print 'Sconstrain', dlg.groupdata['pbflags'].Sconstrain
            print 'Spp', dlg.groupdata['pbflags'].Spp
            print 'Smm', dlg.groupdata['pbflags'].Smm
            print 'Spm', dlg.groupdata['pbflags'].Spm
            print 'Smp', dlg.groupdata['pbflags'].Smp
            print 'text',self.itemdict['text']
            files={}
            for currdata in children_data:
                key=currdata['polstate']
                files[key]=currdata['filename']
            print files
            self.files=files
            self.groupdata=dlg.groupdata

        else:
            print "Cancel"
        dlg.Destroy()

    def OnItemReduceGroup(self,event):
        print 'Reduce Group'
        children_data=self.GetChildData(self.current)
        polstates=[]
        for data in children_data:
            #print data
            self.varying=data['data'].metadata['count_info']['varying']
            self.independent_variable=self.varying[0]
            x=data['data'].data[self.independent_variable]
            polstates.append(data['polstate'])
        files=copy.deepcopy(self.files)
        for ckey,myfile in files.iteritems():
            myfile=os.path.join(os.getcwd(),myfile)+'.bt7'
            files[ckey]=myfile
        #files={}
        #files['pm']=r'C:\polcorrecter\data\fieldscansplusminusreset53630.bt7'
        #files['mp']=r'C:\polcorrecter\data\fieldscanminusplusreset53631.bt7'
        text=self.itemdict['text']+'.polcor'
        print 'driver file', text
        #cellfile='c:\polcorrecter\data\cells.txt'  
        cellfile=self.groupdata['cellfile']
        flags=self.groupdata['pbflags']
        if cellfile !='':
             f=open(text,'wt')
             s='#absolute'
             f.write(s+'\n')
             s='#cell %s'%(cellfile,)
             f.write(s+'\n')

             s='#MonoSelect %d'%(flags.MonoSelect)
             f.write(s+'\n')
             s='#MonitorCorrect %d'%(flags.MonitorCorrect)
             f.write(s+'\n')
             s='#PolMonitorCorrect %d'%(flags.PolMonitorCorrect)
             f.write(s+'\n')

             s='#CountsEnable %d %d %d %d'%(flags.CountsEnable[0],flags.CountsEnable[1],
                                        flags.CountsEnable[2],flags.CountsEnable[3])
             f.write(s+'\n')
             s='#CountsAdd1 %d %d %d %d'%(flags.CountsAdd1[0],flags.CountsAdd1[1],
                                        flags.CountsAdd1[2],flags.CountsAdd1[3])
             f.write(s+'\n')
             s='#CountsAdd2 %d %d %d %d'%(flags.CountsAdd2[0],flags.CountsAdd2[1],
                                        flags.CountsAdd2[2],flags.CountsAdd2[3])
             f.write(s+'\n')
             s='#Sconstrain %d %d %d %d'%(flags.Sconstrain[0],flags.Sconstrain[1],
                                        flags.Sconstrain[2],flags.Sconstrain[3])
             f.write(s+'\n')
             s='#Smm %2.3f %2.3f %2.3f %2.3f'%(flags.Smm[0],flags.Smm[1],
                                        flags.Smm[2],flags.Smm[3])
             f.write(s+'\n')
             s='#Spp %2.3f %2.3f %2.3f %2.3f'%(flags.Spp[0],flags.Spp[1],
                                        flags.Spp[2],flags.Spp[3])
             f.write(s+'\n')
             s='#Smp %2.3f %2.3f %2.3f %2.3f'%(flags.Smp[0],flags.Smp[1],
                                        flags.Smp[2],flags.Smp[3])
             f.write(s+'\n')
             s='#Spm %2.3f %2.3f %2.3f %2.3f'%(flags.Spm[0],flags.Spm[1],
                                        flags.Spm[2],flags.Spm[3])
             f.write(s+'\n')
             s='#directory %s'%(os.getcwd(),)
             s='#absolute'
             f.write(s+'\n')
             s='#begin'
             f.write(s+'\n')
             keylist=['pp','mm','pm','mp']
             s=''
             for key in keylist:
                 if files.has_key(key):
                     s=s+'%s=%s\n'%(key,files[key])
             f.write(s)
             s='#end'
             f.write(s)
             f.close()
             print 'files',files
             print 'cellfile',cellfile
             mypolcor=polcorrect.polarization_correct(files,cellfile)
             corrected_counts=mypolcor.correct(flags)
             mypolcor.savefiles()

             keys=['pp','mm','mp','pm']
             print 'corrected',corrected_counts.keys()
             #print corrected_counts['Spm']
             #print corrected_counts['Smp']
             for key in polstates:
                 if corrected_counts.has_key('S'+key):
                     #x=self.x
                     y=corrected_counts['S'+key]
                     dy=corrected_counts['E'+key]
                     print 'varying',self.independent_variable
                     new_plot = Data1D(x, y, dy=dy)
                     new_plot.name =key+'.corrected'
                     new_plot.group_id='reduced'
                     new_plot.xaxis(str(self.independent_variable), 'A^{-1}')
                     new_plot.yaxis("\\rm{Intensity} ","Arb. units")
                     wx.PostEvent(self.parent.parent, NewPlotEvent(plot=new_plot))

             self.OnItemPlot(event)




    def GetChildData(self,currgroup):
        if self.HasChildren(currgroup):
            treedata=[]
            (item,cookie)=self.GetFirstChild(self.current)
            while item:
                treedata.append(self.GetPyData(item))
                item, cookie=self.GetNextChild(self.current,cookie)
            return treedata
        else:
            return None

    def OnLeftDClick(self, event):

        pt = event.GetPosition()
        item, flags = self.HitTest(pt)
        if item and (flags & CT.TREE_HITTEST_ONITEMLABEL):
            if self.GetTreeStyle() & CT.TR_EDIT_LABELS:
                self.log.write("OnLeftDClick: %s (manually starting label edit)"% self.GetItemText(item) + "\n")
                self.EditLabel(item)
            else:
                self.log.write("OnLeftDClick: Cannot Start Manual Editing, Missing Style TR_EDIT_LABELS\n")

        event.Skip()



    def OnItemBackground(self, event):

        colourdata = wx.ColourData()
        colourdata.SetColour(self.itemdict["back"])
        dlg = wx.ColourDialog(self, colourdata)

        dlg.GetColourData().SetChooseFull(True)

        if dlg.ShowModal() == wx.ID_OK:
            data = dlg.GetColourData()
            col1 = data.GetColour().Get()
            self.SetItemBackgroundColour(self.current, col1)
        dlg.Destroy()


    def OnItemForeground(self, event):

        colourdata = wx.ColourData()
        colourdata.SetColour(self.itemdict["fore"])
        dlg = wx.ColourDialog(self, colourdata)

        dlg.GetColourData().SetChooseFull(True)

        if dlg.ShowModal() == wx.ID_OK:
            data = dlg.GetColourData()
            col1 = data.GetColour().Get()
            self.SetItemTextColour(self.current, col1)
        dlg.Destroy()


    def OnItemBold(self, event):

        self.SetItemBold(self.current, not self.itemdict["isbold"])


    def OnItemFont(self, event):

        data = wx.FontData()
        font = self.itemdict["font"]

        if font is None:
            font = wx.SystemSettings_GetFont(wx.SYS_DEFAULT_GUI_FONT)

        data.SetInitialFont(font)

        dlg = wx.FontDialog(self, data)

        if dlg.ShowModal() == wx.ID_OK:
            data = dlg.GetFontData()
            font = data.GetChosenFont()
            self.SetItemFont(self.current, font)

        dlg.Destroy()


    def OnItemHyperText(self, event):

        self.SetItemHyperText(self.current, not self.itemdict["ishtml"])


    def OnEnableWindow(self, event):

        enable = self.GetItemWindowEnabled(self.current)
        self.SetItemWindowEnabled(self.current, not enable)


    def OnDisableItem(self, event):

        self.EnableItem(self.current, False)


    def OnItemIcons(self, event):

        bitmaps = [self.itemdict["normal"], self.itemdict["selected"],
                   self.itemdict["expanded"], self.itemdict["selexp"]]

        wx.BeginBusyCursor()
        dlg = TreeIcons(self, -1, bitmaps=bitmaps)
        wx.EndBusyCursor()
        dlg.ShowModal()


    def SetNewIcons(self, bitmaps):

        self.SetItemImage(self.current, bitmaps[0], CT.TreeItemIcon_Normal)
        self.SetItemImage(self.current, bitmaps[1], CT.TreeItemIcon_Selected)
        self.SetItemImage(self.current, bitmaps[2], CT.TreeItemIcon_Expanded)
        self.SetItemImage(self.current, bitmaps[3], CT.TreeItemIcon_SelectedExpanded)


    def OnItemInfo(self, event):

        itemtext = self.itemdict["text"]
        numchildren = str(self.itemdict["children"])
        itemtype = self.itemdict["itemtype"]
        pydata = repr(type(self.itemdict["pydata"]))

        if itemtype == 0:
            itemtype = "Normal"
        elif itemtype == 1:
            itemtype = "CheckBox"
        else:
            itemtype = "RadioButton"

        strs = "Information On Selected Item:\n\n" + "Text: " + itemtext + "\n" \
               "Number Of Children: " + numchildren + "\n" \
               "Item Type: " + itemtype + "\n" \
               "Item Data Type: " + pydata + "\n"

        dlg = wx.MessageDialog(self, strs, "CustomTreeCtrlDemo Info", wx.OK | wx.ICON_INFORMATION)
        dlg.ShowModal()
        dlg.Destroy()



    def OnItemDelete(self, event):

        strs = "Are You Sure You Want To Delete Item " + self.GetItemText(self.current) + "?"
        dlg = wx.MessageDialog(None, strs, 'Deleting Item', wx.YES_NO | wx.NO_DEFAULT | wx.CANCEL | wx.ICON_QUESTION)

        if dlg.ShowModal() in [wx.ID_NO, wx.ID_CANCEL]:
            dlg.Destroy()
            return

        dlg.Destroy()

        self.DeleteChildren(self.current)
        self.Delete(self.current)
        self.current = None



    def OnItemPrepend(self, event):

        dlg = wx.TextEntryDialog(self, "Please Enter The New Item Name", 'Item Naming', 'Python')

        if dlg.ShowModal() == wx.ID_OK:
            newname = dlg.GetValue()
            newitem = self.PrependItem(self.current, newname)
            self.EnsureVisible(newitem)

        dlg.Destroy()


    def OnItemAppend(self, event):

        dlg = wx.TextEntryDialog(self, "Please Enter The New Item Name", 'Item Naming', 'Python')

        if dlg.ShowModal() == wx.ID_OK:
            newname = dlg.GetValue()
            newitem = self.AppendItem(self.current, newname)
            self.EnsureVisible(newitem)

        dlg.Destroy()


    def OnBeginEdit(self, event):

        self.log.write("OnBeginEdit" + "\n")
        # show how to prevent edit...
        item = event.GetItem()
        if item and self.GetItemText(item) == "The Root Item":
            wx.Bell()
            self.log.write("You can't edit this one..." + "\n")

            # Lets just see what's visible of its children
            cookie = 0
            root = event.GetItem()
            (child, cookie) = self.GetFirstChild(root)

            while child:
                self.log.write("Child [%s] visible = %d" % (self.GetItemText(child), self.IsVisible(child)) + "\n")
                (child, cookie) = self.GetNextChild(root, cookie)

            event.Veto()


    def OnEndEdit(self, event):

        self.log.write("OnEndEdit: %s %s" %(event.IsEditCancelled(), event.GetLabel()))
        # show how to reject edit, we'll not allow any digits
        for x in event.GetLabel():
            if x in string.digits:
                self.log.write(", You can't enter digits..." + "\n")
                event.Veto()
                return

        self.log.write("\n")



    def OnItemExpanded(self, event):

        item = event.GetItem()
        if item:
            self.log.write("OnItemExpanded: %s" % self.GetItemText(item) + "\n")


    def OnItemExpanding(self, event):

        item = event.GetItem()
        if item:
            self.log.write("OnItemExpanding: %s" % self.GetItemText(item) + "\n")

        event.Skip()


    def OnItemCollapsed(self, event):

        item = event.GetItem()
        if item:
            self.log.write("OnItemCollapsed: %s" % self.GetItemText(item) + "\n")


    def OnItemCollapsing(self, event):

        item = event.GetItem()
        if item:
            self.log.write("OnItemCollapsing: %s" % self.GetItemText(item) + "\n")

        event.Skip()


    def OnSelChanged(self, event):

        self.item = event.GetItem()
        if self.item:
            self.log.write("OnSelChanged: %s" % self.GetItemText(self.item))
            if wx.Platform == '__WXMSW__':
                self.log.write(", BoundingRect: %s" % self.GetBoundingRect(self.item, True) + "\n")
            else:
                self.log.write("\n")

        event.Skip()


    def OnSelChanging(self, event):

        item = event.GetItem()
        olditem = event.GetOldItem()

        if item:
            if not olditem:
                olditemtext = "None"
            else:
                olditemtext = self.GetItemText(olditem)
            self.log.write("OnSelChanging: From %s" % olditemtext + " To %s" % self.GetItemText(item) + "\n")

        event.Skip()


    def OnBeginDrag(self, event):

        self.item = event.GetItem()
        if self.item:
            self.log.write("Beginning Drag..." + "\n")

            event.Allow()


    def OnBeginRDrag(self, event):

        self.item = event.GetItem()
        if self.item:
            self.log.write("Beginning Right Drag..." + "\n")

            event.Allow()


    def OnEndDrag(self, event):

        self.item = event.GetItem()
        if self.item:
            self.log.write("Ending Drag!" + "\n")

        event.Skip()


    def OnDeleteItem(self, event):

        item = event.GetItem()

        if not item:
            return

        self.log.write("Deleting Item: %s" % self.GetItemText(item) + "\n")
        event.Skip()


    def OnItemCheck(self, event):

        item = event.GetItem()
        self.log.write("Item " + self.GetItemText(item) + " Has Been Checked!\n")
        event.Skip()


    def OnItemChecking(self, event):

        item = event.GetItem()
        self.log.write("Item " + self.GetItemText(item) + " Is Being Checked...\n")
        event.Skip()


    def OnToolTip(self, event):

        item = event.GetItem()
        if item:
            event.SetToolTip(wx.ToolTip(self.GetItemText(item)))


    def OnItemMenu(self, event):

        item = event.GetItem()
        if item:
            self.log.write("OnItemMenu: %s" % self.GetItemText(item) + "\n")

        event.Skip()


    def OnKey(self, event):

        keycode = event.GetKeyCode()
        keyname = keyMap.get(keycode, None)

        if keycode == wx.WXK_BACK:
            self.log.write("OnKeyDown: HAHAHAHA! I Vetoed Your Backspace! HAHAHAHA\n")
            return

        if keyname is None:
            if "unicode" in wx.PlatformInfo:
                keycode = event.GetUnicodeKey()
                if keycode <= 127:
                    keycode = event.GetKeyCode()
                keyname = "\"" + unichr(event.GetUnicodeKey()) + "\""
                if keycode < 27:
                    keyname = "Ctrl-%s" % chr(ord('A') + keycode-1)

            elif keycode < 256:
                if keycode == 0:
                    keyname = "NUL"
                elif keycode < 27:
                    keyname = "Ctrl-%s" % chr(ord('A') + keycode-1)
                else:
                    keyname = "\"%s\"" % chr(keycode)
            else:
                keyname = "unknown (%s)" % keycode

        self.log.write("OnKeyDown: You Pressed '" + keyname + "'\n")

        event.Skip()


    def OnActivate(self, event):

        if self.item:
            self.log.write("OnActivate: %s" % self.GetItemText(self.item) + "\n")

        event.Skip()


    def OnHyperLink(self, event):

        item = event.GetItem()
        if item:
            self.log.write("OnHyperLink: %s" % self.GetItemText(self.item) + "\n")


    def OnTextCtrl(self, event):

        char = chr(event.GetKeyCode())
        self.log.write("EDITING THE TEXTCTRL: You Wrote '" + char + \
                       "' (KeyCode = " + str(event.GetKeyCode()) + ")\n")
        event.Skip()


    def OnComboBox(self, event):

        selection = event.GetEventObject().GetValue()
        self.log.write("CHOICE FROM COMBOBOX: You Chose '" + selection + "'\n")
        event.Skip()

class FileTreePanel(wx.Panel):

    ## Internal name for the AUI manager
    window_name = "filetreepanel"
    ## Title to appear on top of the window
    window_caption = "File Tree Panel"

    def __init__(self,parent,id):
        self.parent=parent
        wx.Panel.__init__(self,parent,id,style=0)
        mytree=CustomTreeCtrl(self,-1,style=wx.SUNKEN_BORDER | CT.TR_HAS_BUTTONS | CT.TR_HAS_VARIABLE_ROW_HEIGHT\
        | CT.TR_HIDE_ROOT|CT.TR_TWIST_BUTTONS|CT.TR_EDIT_LABELS| wx.WS_EX_VALIDATE_RECURSIVELY\
        |CT.TR_MULTIPLE|CT.TR_EXTENDED,log=sys.stdout)
        bs = wx.BoxSizer(wx.VERTICAL)
        bs.Add(mytree, 1, wx.GROW|wx.ALL|wx.EXPAND, 5)
        self.SetSizer(bs)
        #bs.Fit()
        self.tooltip = ''
        self.tree=mytree



class FileTreeFrame(wx.Frame):
    def __init__(self,parent,id):
        wx.Frame.__init__(self,parent,id,'File Groups',size=(300,600),style=wx.DEFAULT_FRAME_STYLE^wx.CLOSE_BOX)
        self.Bind(wx.EVT_CLOSE,self.OnCloseWindow)
        self.filetree_panel = FileTreePanel(self, -1)
        #self.dataplot_frame=dataplot.TestFrame(self,-1)
        #self.dataplot_frame.Show()


    def OnCloseWindow(self,event):
        self.Destroy()



if __name__=='__main__':
    app=MyApp()
    frame=FileTreeFrame(parent=None,id=-1)
    frame.Show()
    app.MainLoop()