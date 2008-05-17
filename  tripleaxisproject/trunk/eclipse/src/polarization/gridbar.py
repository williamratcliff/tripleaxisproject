import wx
import wx.grid
import copy
import numpy as N
from pickle import loads,dumps
threshold=1e-3

class Range:
    def __init__(self,myrange):
        self.low=myrange[0]#900
        self.high=myrange[1]#-900
    def map(self,value):
        if self.low>=self.high:
            return 0.5
        else:
            return (value-self.low)/(self.high-self.low)


class Bar:
    def __init__(self,low,high,myrange):
        self.low=low
        self.high=high
        self.center=(self.low+self.high)/2
        self.range=myrange


    def __str__(self):
        return "Bar(%g,%g)"%(self.low,self.high)

    def __cmp__(self,other):
        if not isinstance(other,Bar): return -1
        if self.low <other.low:
            return -1
        elif self.low > other.low:
            return 1
        elif self.low==other.low:
            return cmp(self.high,other.high)


GRID_VALUE_BAR=Bar.__name__

def register(grid):
    grid.RegisterDataType(GRID_VALUE_BAR, GridCellBarRenderer(), None)

class GridCellBarRenderer(wx.grid.PyGridCellRenderer):
    def __init__(self,*args, **kwargs):
        wx.grid.PyGridCellRenderer.__init__(self,*args,**kwargs)
        #print 'hello'

    def DrawBar(self, grid, attr, dc, rect, row, col, isSelected):
        #print 'Drawing'
        #print str(grid.GetCellValue(row, col))
        #bar = loads(str(grid.GetCellValue(row, col)))
        #bar = Bar.getcell(int(grid.GetCellValue(row,col)))
        table=grid.GetTable()
        bar = table.GetValue(row,col)
        #print 'bar ',bar
        if bar=='':
            return
        low=bar.low
        high=bar.high
        bar_range=bar.range
        #print 'low',bar.low,'high',bar.high,' range', bar_range
        cell_height=rect.height
        cell_width=rect.width
        cell_x=rect.x
        cell_y=rect.y
        rect.height=rect.height/4
        rect.y=rect.y+2*rect.height
        rect_b=copy.deepcopy(rect)
        fill=attr.GetTextColour()
        background="white"
        outline=background#"Navy"
        dc.SetBrush(wx.Brush(background,wx.SOLID))
        dc.SetPen(wx.Pen(outline))
        dc.DrawRectangleRect(rect_b)
        rangefinder=Range(bar_range)
        low_px=rangefinder.map(low)*cell_width
        high_px=rangefinder.map(high)*cell_width
        rect.width=(high_px-low_px)
        if rect.width<1:
            rect.width=1
        rect.x=cell_x+low_px
        #print 'low_mapped',low_mapped,'high_mapped',high_mapped,'width',rect.width,'rect.x',rect.x
        dc.SetBrush(wx.Brush(fill,wx.SOLID))
        dc.SetPen(wx.Pen(fill))
        dc.DrawRectangleRect(rect)
        s1='%4.3f'%(low,)
        s2='%4.3f'%(high,)
        dc.SetFont(wx.Font(7, wx.DEFAULT, wx.NORMAL, wx.NORMAL))
        text_width, text_height = dc.GetTextExtent(s2)

        #dc.DrawText(s1,cell_x,cell_y)
        #dc.DrawText(s2,cell_x+cell_width-text_width,cell_y)

    def Draw(self, grid, attr, dc, rect, row, col, isSelected):
        colour=attr.GetBackgroundColour()

        dc.SetBrush(wx.Brush(colour,wx.SOLID))
        dc.SetPen(wx.Pen(colour))
        dc.DrawRectangleRect(rect)
        self.DrawBar(grid, attr, dc, rect, row, col, isSelected)


    def GetBestSize(self, grid, attr, dc, row, col):
        #print 'Getting Best Size'
        dc.SetFont(attr.GetFont())
        w, h = dc.GetTextExtent('999.999 9999.999')
        return wx.Size(w, h)

    def Clone(self):
        return GridCellBarRenderer()



