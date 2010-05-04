#!/usr/bin/env python
# Copyright (c) 2007-8 Qtrac Ltd. All rights reserved.
# This program or module is free software: you can redistribute it and/or
# modify it under the terms of the GNU General Public License as published
# by the Free Software Foundation, either version 2 of the License, or
# version 3 of the License, or (at your option) any later version. It is
# provided for educational purposes and is distributed in the hope that
# it will be useful, but WITHOUT ANY WARRANTY; without even the implied
# warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See
# the GNU General Public License for more details.

from __future__ import division
import gzip
import os
import platform
import sys
from PyQt4.QtCore import *
from PyQt4.QtGui import *
import genfiles
import numpy as N


(FILE_NUMBER, H, K, L, TTHETA, OMEGA) = range(6)


TIMESTAMPFORMAT = "yyyy-MM-dd hh:mm"

class Range(object):
    def __init__(self,myrange):
        self.low=myrange[0]#900
        self.high=myrange[1]#-900
    def map(self,value):
        if self.low>=self.high:
            return 0.5
        else:
            return (value-self.low)/(self.high-self.low)



class DataModel(QAbstractTableModel):

    def __init__(self, filename):
        super(DataModel, self).__init__()
        self.filename = filename
        self.dataset = None


    def load(self):
        myfilebase=str(19262)
        myend='dat'
        mydirectory=r'c:\tbmno3\aug25_2009_ill'
        self.dataset=genfiles.read_files(mydirectory,myfilebase,myend)
        self.reset()
        


    def data(self, index, role=Qt.DisplayRole):  
        try:
            mylen=len(self.dataset.hmin)
        except:
            mylen=0
        if not index.isValid() or \
           not (0 <= index.row() < mylen):
            return QVariant()
        column = index.column()
        #result = self.results[index.row()]
        row=index.row()
        if role == Qt.DisplayRole:
            #item = result[column]
            #if column == TIMESTAMP:
            #    item = item.toString(TIMESTAMPFORMAT)
            #else:
            #    item = QString("%1").arg(item, 0, "f", 2)
            if column==FILE_NUMBER:
                item=self.dataset.filenums[row]
                item = QString("%1").arg(int(item), 0)
            elif column==H:
                item=self.dataset.hmin[row]
                item = QString("%1").arg(item, 0, "f", 2)
            elif column==K:
                item=self.dataset.kmin[row]
                item = QString("%1").arg(item, 0, "f", 2)
            elif column==L:
                item=self.dataset.lmin[row]
                item = QString("%1").arg(item, 0, "f", 2)
            elif column==OMEGA:
                #item=(min(self.dataset.data[row].angle1)+max(self.dataset.data[row].angle1))/2
                #item = QString("%1").arg(item, 0, "f", 2)
                item=[self.dataset.omega_min[row],self.dataset.omega_max[row], self.dataset.omega_minimum,self.dataset.omega_maximum]
            elif column==TTHETA:
                item=self.dataset.tth[row]
                item = QString("%1").arg(item, 0, "f", 2)    
                #if not self.dataset.data[row].angle2==[]:
                #    item=(min(self.dataset.tth_min[row])+max(self.dataset.tth_max[row]))/2
                #    item = QString("%1").arg(item, 0, "f", 2)
                #else:
                #    item=self.dataset.tth[row]
                #    item = QString("%1").arg(item, 0, "f", 2)            
            return QVariant(item)
        elif role == Qt.TextAlignmentRole:
            return QVariant(int(Qt.AlignRight|Qt.AlignVCenter))
        elif role == Qt.BackgroundColorRole:
            return QVariant(QColor(Qt.white))
        
        elif role == Qt.ToolTipRole:
            if column==OMEGA:
                item1=self.dataset.omega_min[row]
                item2=self.dataset.omega_max[row]
                tip = QString("%1 %2").arg(item1, 0, "f", 2).arg(item2, 0, "f", 2)
                return(QVariant(tip))
            
            #if column != TIMESTAMP:
            #    return QVariant(int(Qt.AlignRight|Qt.AlignVCenter))
            #return QVariant(int(Qt.AlignLeft|Qt.AlignVCenter))
        #elif role == Qt.TextColorRole and column == INLETFLOW:
        #    if result[column] < 0:
        #        return QVariant(QColor(Qt.red))
        #elif role == Qt.TextColorRole and \
        #     column in (RAWPH, FLOCCULATEDPH):
        #    ph = result[column]
        #    if ph < 7:
        #        return QVariant(QColor(Qt.red))
        #    elif ph >= 8:
        #        return QVariant(QColor(Qt.blue))
        #    else:
        #        return QVariant(QColor(Qt.darkGreen))
        return QVariant()


    def headerData(self, section, orientation, role=Qt.DisplayRole):
        if role == Qt.TextAlignmentRole:
            if orientation == Qt.Horizontal:
                return QVariant(int(Qt.AlignCenter))
            return QVariant(int(Qt.AlignRight|Qt.AlignVCenter))
        if role != Qt.DisplayRole:
            return QVariant()
        if orientation == Qt.Horizontal:
            if section == FILE_NUMBER:
                return QVariant("File#")
            elif section == H:
                return QVariant("H")
            elif section == K:
                return QVariant("K")
            elif section == L:
                return QVariant("L")
            elif section == OMEGA:
                return QVariant(u"\u00B0" +"Omega")
            elif section == TTHETA:
                return QVariant(u"\u00B0" +"2Theta")

        return QVariant(int(section + 0))  #Start counting rows from 0


    def rowCount(self, index=QModelIndex()):
        try:
            return len(self.dataset.hmin)
        except:
            return 0


    def columnCount(self, index=QModelIndex()):
        return 6
    
    def sortByCountryOwner(self):
        def compare(a, b):
            if a.country != b.country:
                return QString.localeAwareCompare(a.country, b.country)
            if a.owner != b.owner:
                return QString.localeAwareCompare(a.owner, b.owner)
            return QString.localeAwareCompare(a.name, b.name)
        self.ships = sorted(self.ships, compare)
        self.reset()
        
    def sortByName(self,col):
        #self.ships = sorted(self.ships)
        coldict={}
        coldict[FILE_NUMBER]=self.dataset.filenums
        coldict[H]=self.dataset.hmin
        coldict[K]=self.dataset.kmin
        coldict[L]=self.dataset.lmin
        coldict[TTHETA]=self.dataset.tth
        coldict[OMEGA]=self.dataset.omega_min
        coldict[OMEGA+1]=self.dataset.omega_max
        
        
        col_to_sort=[(i,s) for i,s in enumerate(coldict[col])]
        col_to_sort.sort(lambda x,y: cmp(x[1],y[1]))
        g_col = [i for (i,s) in col_to_sort]
        #print col_to_sort
        if col >=0:
            if (N.diff(g_col)>0).all():
                g_col=g_col[::-1]

            #print 'col=',col
            #print 'sort '
            #print g
            for i in range(7):
                coldict[i][:]=N.array(coldict[i])[g_col]
                #data[:,i]=data[g_col,i]
            
        self.reset()
        
    def sortByRange(self,col):
        coldict={}
        coldict[FILE_NUMBER]=self.dataset.filenums
        coldict[H]=self.dataset.hmin
        coldict[K]=self.dataset.kmin
        coldict[L]=self.dataset.lmin
        coldict[TTHETA]=self.dataset.tth
        coldict[OMEGA]=self.dataset.omega_min
        coldict[OMEGA+1]=self.dataset.omega_max
        
        
        col_to_sort=[(i,s) for i,s in enumerate(coldict[col])]
        col_to_sort.sort(lambda x,y: cmp(x[1],y[1]))
        g_col = [i for (i,s) in col_to_sort]
        #print col_to_sort
        if col >=0:
            if (N.diff(g_col)>0).all():
                g_col=g_col[::-1]

            #print 'col=',col
            #print 'sort '
            #print g
            for i in range(7):
                coldict[i][:]=N.array(coldict[i])[g_col]
                #data[:,i]=data[g_col,i]
        self.reset()


class DataDelegate(QItemDelegate):

    def __init__(self, parent=None):
        super(DataDelegate, self).__init__(parent)
        
        
    def paint(self, painter, option, index):
        palette = QApplication.palette()
        if index.column() == FILE_NUMBER:
            QItemDelegate.paint(self, painter, option, index)
        elif index.column() in [OMEGA]:
            color = palette.highlight().color() \
                if option.state & QStyle.State_Selected \
                else QColor(index.model().data(index,
                        Qt.BackgroundColorRole))
            painter.save()
            painter.fillRect(option.rect, color)
            #painter.translate(option.rect.x(), option.rect.y())
            fm = option.fontMetrics
            size=fm.height()
            painter.setRenderHint(QPainter.Antialiasing)
            painter.setRenderHint(QPainter.TextAntialiasing)
            painter.setPen(Qt.NoPen)
            height=size
            color=QColor(0,0,0)
            painter.setBrush(color)
            cell_width=option.rect.width()
            low,high,om_min,om_max = index.model().data(index).toPyObject()
            #angle1_range= N.absolute( max(angle1)-min(angle1) )
            #low,high=min(angle1), max(angle1)
            #bar_range=[min(angle1),max(angle1)]
            bar_range=[om_min,om_max]
            rangefinder=Range(bar_range)
            low_px=rangefinder.map(low)*cell_width
            high_px=rangefinder.map(high)*cell_width
            
            new_width=high_px-low_px
            #set a minimum size of one character
            #if new_width<fm.width("9"):
            #    new_width=fm.width("9")
            if new_width < 1:
                new_width=1
            new_x=low_px#0.2/width
            x=option.rect.x()+new_x
            y=option.rect.y()+size
            painter.drawRect(x,y,new_width,height/3)
            painter.restore()
            
            
        elif index.column() in [H,K,L,TTHETA]:
            QItemDelegate.paint(self, painter, option, index)


    def paint2(self, painter, option, index):
        if index.column() == DESCRIPTION:
            text = index.model().data(index).toString()
            palette = QApplication.palette()
            document = QTextDocument()
            document.setDefaultFont(option.font)
            if option.state & QStyle.State_Selected:
                document.setHtml(QString("<font color=%1>%2</font>") \
                        .arg(palette.highlightedText().color().name())\
                        .arg(text))
            else:
                document.setHtml(text)
            color = palette.highlight().color() \
                if option.state & QStyle.State_Selected \
                else QColor(index.model().data(index,
                        Qt.BackgroundColorRole))
            painter.save()
            painter.fillRect(option.rect, color)
            painter.translate(option.rect.x(), option.rect.y())
            document.drawContents(painter)
            painter.restore()
        else:
            QItemDelegate.paint(self, painter, option, index)


    def sizeHint(self, option, index):
        fm = option.fontMetrics
        if index.column() == FILE_NUMBER:
            return QSize(fm.width("192612"), fm.height())
        else:
            return QSize(fm.width("192612"), fm.height())            
        #else:
            #text = index.model().data(index).toString()
            #document = QTextDocument()
            #document.setDefaultFont(option.font)
            #document.setHtml(text)
            #return QSize(document.idealWidth() + 5, fm.height())
            
        return QItemDelegate.sizeHint(self, option, index)
    
    
    



class MainForm(QDialog):

    def __init__(self, parent=None):
        super(MainForm, self).__init__(parent)

        self.model = DataModel(os.path.join(
                os.path.dirname(__file__), "waterdata.csv.gz"))
        self.tableView = QTableView()
        self.tableView.setAlternatingRowColors(True)
        self.tableView.setModel(self.model)
        self.tableView.setItemDelegate(DataDelegate(self))
        self.tableView.resizeColumnsToContents()
        header = self.tableView.horizontalHeader()
        self.connect(header, SIGNAL("sectionClicked(int)"),
                         self.sortTable)
        self.connect(self.tableView, SIGNAL("clicked(const QModelIndex&)"), 
                     self._onClick) 
        
        
        #self.waterView = DataView()
        #self.waterView.setModel(self.model)
        scrollArea = QScrollArea()
        scrollArea.setBackgroundRole(QPalette.Light)
        #scrollArea.setWidget(self.waterView)
        #self.waterView.scrollarea = scrollArea

        #splitter = QSplitter(Qt.Horizontal)
        #splitter.addWidget(self.tableView)
        #splitter.addWidget(scrollArea)
        #splitter.setSizes([600, 250])
        
        
        
        layout = QVBoxLayout()
        #layout.addWidget(splitter)
        layout.addWidget(self.tableView)
        #layout.addWidget(scrollArea)
        self.setLayout(layout)

        self.setWindowTitle("Data Table")
        QTimer.singleShot(0, self.initialLoad)


    def initialLoad(self):
        QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))
        splash = QLabel(self)
        pixmap = QPixmap(os.path.join(os.path.dirname(__file__),
                "iss013-e-14802.jpg"))
        splash.setPixmap(pixmap)
        splash.setWindowFlags(Qt.SplashScreen)
        splash.move(self.x() + ((self.width() - pixmap.width()) / 2),
                    self.y() + ((self.height() - pixmap.height()) / 2))
        splash.show()
        QApplication.processEvents()
        try:
            self.model.load()
            #self.tableView.reset()
        except IOError, e:
            QMessageBox.warning(self, "Water Quality - Error", e)
        else:
            self.tableView.resizeColumnsToContents()
        splash.close()
        QApplication.processEvents()
        QApplication.restoreOverrideCursor()
        
    def sortTable(self, section):
        if section in (H,K,L,TTHETA,FILE_NUMBER):
            self.model.sortByName(section)
        elif section == OMEGA:
            self.model.sortByRange(section)
        self.tableView.resizeColumnsToContents()
        
    def _onClick(self, index):
        print "_onClick", index
        column = index.column()
        if column==FILE_NUMBER:
            item=str(index.model().data(index).toString())+'.dat'
            print 'item',item
            self.emit(SIGNAL("fileClicked"),item)
            
        

if __name__=="__main__":
    app = QApplication(sys.argv)
    form = MainForm()
    form.resize(850, 620)
    form.show()
    app.exec_()

