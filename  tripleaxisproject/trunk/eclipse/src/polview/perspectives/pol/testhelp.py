# read PDF files (.pdf) with wxPython
# using wx.lib.pdfwin.PDFWindow class ActiveX control
# from wxPython's new wx.activex module, this allows one
# to use an ActiveX control, as if it would be wx.Window
# it embeds the Adobe Acrobat Reader
# as far as HB knows this works only with Windows
# tested with Python24 and wxPython26 by HB

import  wx

if wx.Platform == '__WXMSW__':
    from wx.lib.pdfwin import PDFWindow


class MyPanel(wx.Panel):
    def __init__(self, parent):
        wx.Panel.__init__(self, parent, id=-1)
        self.pdf = None

        sizer = wx.BoxSizer(wx.VERTICAL)
        btnSizer = wx.BoxSizer(wx.HORIZONTAL)

        self.pdf = PDFWindow(self, style=wx.SUNKEN_BORDER)

        sizer.Add(self.pdf, proportion=1, flag=wx.EXPAND)

        btn = wx.Button(self, wx.NewId(), "Open PDF File")
        self.Bind(wx.EVT_BUTTON, self.OnOpenButton, btn)
        btnSizer.Add(btn, proportion=1, flag=wx.EXPAND|wx.ALL, border=5)

        btn = wx.Button(self, wx.NewId(), "Previous Page")
        self.Bind(wx.EVT_BUTTON, self.OnPrevPageButton, btn)
        btnSizer.Add(btn, proportion=1, flag=wx.EXPAND|wx.ALL, border=5)

        btn = wx.Button(self, wx.NewId(), "Next Page")
        self.Bind(wx.EVT_BUTTON, self.OnNextPageButton, btn)
        btnSizer.Add(btn, proportion=1, flag=wx.EXPAND|wx.ALL, border=5)

        btnSizer.Add((50,-1), proportion=2, flag=wx.EXPAND)
        sizer.Add(btnSizer, proportion=0, flag=wx.EXPAND)

        self.SetSizer(sizer)
        self.SetAutoLayout(True)

    def OnOpenButton(self, event):
        # make sure you have PDF files available on your drive
        dlg = wx.FileDialog(self, wildcard="*.pdf")
        if dlg.ShowModal() == wx.ID_OK:
            wx.BeginBusyCursor()
            self.pdf.LoadFile(dlg.GetPath())
            wx.EndBusyCursor()
        dlg.Destroy()

    def OnPrevPageButton(self, event):
        self.pdf.gotoPreviousPage()

    def OnNextPageButton(self, event):
        self.pdf.gotoNextPage()


app = wx.PySimpleApp()
# create window/frame, no parent, -1 is default ID, title, size
frame = wx.Frame(None, -1, "PDFWindow", size = (640, 480))
# make an instance of the class
MyPanel(frame)
# show the frame
frame.Show(True)
# start the event loop
app.MainLoop()