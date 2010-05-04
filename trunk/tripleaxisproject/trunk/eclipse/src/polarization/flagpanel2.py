import wx
import wxaddons.sized_controls as sc
import  wx.lib.intctrl

class MyApp(wx.App):
    def __init__(self, redirect=False, filename=None, useBestVisual=False, clearSigInt=True):
        wx.App.__init__(self,redirect,filename,clearSigInt)


    def OnInit(self):
        return True


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
        wx.CheckBox(CountsAdd1Pane, -1, "NSF")
        wx.CheckBox(CountsAdd1Pane, -1, "SF")



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
    frame=FormDialog(parent=None,id=-1)
    frame.Show()
    app.MainLoop()
