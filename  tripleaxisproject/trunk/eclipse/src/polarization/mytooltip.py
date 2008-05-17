self.tooltip = ''
self.myGrid.GetGridWindow().Bind(wx.EVT_MOTION, self.onMouseOver)

def onMouseOver(self, event):
    '''
    Method to calculate where the mouse is pointing and
    then set the tooltip dynamically.
    '''

    # Use CalcUnscrolledPosition() to get the mouse position within the
    # entire grid including what's offscreen
    x, y =
self.totals_sheet.CalcUnscrolledPosition(event.GetX(),event.GetY())

    coords = self.totals_sheet.XYToCell(x, y)
#coords = grid.XYToCell(x, y)
    col = coords[1]

    # Example colum limit to apply the custom tooltip to
    if col == 16:
        row = coords[0]
        val = float(self.totals_sheet.GetCellValue(row, col))
        if val != '':
            hrs = val * 10.0
            event.GetEventObject().SetToolTipString('Now the total = %s' %
hrs)
            self.tooltip = 'Now the total = %s' % hrs
    else:
        event.GetEventObject().SetToolTipString('')
        self.tooltip = ''