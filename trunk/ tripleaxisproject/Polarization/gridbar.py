import wx
import wx.grid
import copy

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
        dc.Clear()
        #dc.BeginDrawing()
        cell_height=rect.height
        cell_width=rect.width
        cell_x=rect.x
        cell_y=rect.y
        rect.height=rect.height/4
        rect.width=rect.width/3
        rect.y=rect.y+2*rect.height
        rect.x=cell_x+cell_width/4
        dc.SetBrush(wx.Brush("navy",wx.SOLID))
        dc.SetPen(wx.Pen("navy"))
        dc.DrawRectangleRect(rect)
        dc.DrawText('h',cell_x,cell_y)
        print rect.x, rect.y
        #dc.EndDrawing()
        print rect

    def GetBestSize(self, grid, attr, dc, row, col):
        print 'Getting Best Size'
        dc.SetFont(attr.GetFont())
        w, h = dc.GetTextExtent('mmmmm')
        return wx.Size(w, h)

    def Clone(self):
        return GridCellBarRenderer()



