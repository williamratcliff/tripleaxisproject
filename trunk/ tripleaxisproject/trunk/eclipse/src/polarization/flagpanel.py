import wx
import wxaddons.sized_controls as sc


class MyApp(wx.App):
    def __init__(self, redirect=False, filename=None, useBestVisual=False, clearSigInt=True):
        wx.App.__init__(self,redirect,filename,clearSigInt)


    def OnInit(self):
        return True


class FormDialog(sc.SizedDialog):
    def __init__(self, parent, id):
        sc.SizedDialog.__init__(self, None, -1, "SizedForm Dialog",
                        style=wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER)

        pane = self.GetContentsPane()
        pane.SetSizerType("form")

        # row 1
        wx.StaticText(pane, -1, "Name")
        textCtrl = wx.TextCtrl(pane, -1, "Your name here")
        textCtrl.SetSizerProps(expand=True)

        # row 2
        wx.StaticText(pane, -1, "Email")
        emailCtrl = wx.TextCtrl(pane, -1, "")
        emailCtrl.SetSizerProps(expand=True)

        # row 3
        wx.StaticText(pane, -1, "Gender")
        wx.Choice(pane, -1, choices=["male", "female"])

        # row 4
        wx.StaticText(pane, -1, "State")
        wx.TextCtrl(pane, -1, size=(60, -1)) # two chars for state

        # row 5
        wx.StaticText(pane, -1, "Title")

        # here's how to add a 'nested sizer' using sized_controls
        radioPane = sc.SizedPanel(pane, -1)
        radioPane.SetSizerType("horizontal")
        radioPane.SetSizerProps(expand=True)

        # make these children of the radioPane to have them use
        # the horizontal layout
        wx.RadioButton(radioPane, -1, "Mr.")
        wx.RadioButton(radioPane, -1, "Mrs.")
        wx.RadioButton(radioPane, -1, "Dr.")
        # end row 5

        # add dialog buttons
        self.SetButtonSizer(self.CreateStdDialogButtonSizer(wx.OK | wx.CANCEL))

        # a little trick to make sure that you can't resize the dialog to
        # less screen space than the controls need
        self.Fit()
        self.SetMinSize(self.GetSize())





if __name__=='__main__':
    app=MyApp()
    frame=FormDialog(parent=None,id=-1)
    frame.Show()
    app.MainLoop()
