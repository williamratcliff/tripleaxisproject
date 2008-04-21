import wx
import wx.grid

class Bar:
    def __init__(self,low,high,range):
        self.low=low
        self.high=high
        self.range=range

GRID_VALUE_BAR=Bar.__name__

def register(grid):
    grid.RegisterDataType(GRID_VALUE_BAR, GridCellBarRenderer(), None)

class GridCellBarRenderer(wx.grid.PyGridCellRenderer):
    def __init__(self,*args, **kwargs):
        wx.grid.PyGridCellRenderer.__init__(self,*args,**kwargs)
        print 'hello'

    def Draw(self, grid, attr, dc, rect, row, col, isSelected):
        print 'Drawing'
        range = grid.GetCellValue(row, col)
        #dc=wx.DC
        #dc.Clear()
        #dc.BeginDrawing()
        dc.SetBrush(wx.Brush("navy",wx.SOLID))
        dc.SetPen(wx.Pen("red"))
        dc.DrawRectangleRect(rect)
        #dc.EndDrawing()
        print rect

    def GetBestSize(self, grid, attr, dc, row, col):
        print 'Getting Best Size'
        dc.SetFont(attr.GetFont())
        w, h = dc.GetTextExtent('mmmmm')
        return wx.Size(w, h)

    def Clone(self):
        return GridCellBarRenderer()



