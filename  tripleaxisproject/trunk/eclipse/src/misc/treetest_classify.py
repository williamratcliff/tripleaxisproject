"""***************************************************************************
**
** Copyright (C) 2005-2005 Trolltech AS. All rights reserved.
**
** This file is part of the example classes of the Qt Toolkit.
**
** This file may be used under the terms of the GNU General Public
** License version 2.0 as published by the Free Software Foundation
** and appearing in the file LICENSE.GPL included in the packaging of
** this file.  Please review the following information to ensure GNU
** General Public Licensing requirements will be met:
** http://www.trolltech.com/products/qt/opensource.html
**
** If you are unsure which license is appropriate for your use, please
** review the following information:
** http://www.trolltech.com/products/qt/licensing.html or contact the
** sales department at sa...@trolltech.com.
**
** This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
** WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
**
***************************************************************************"""

import sys
from PyQt4 import QtCore, QtGui
from PyQt4.examples.itemviews.simpletreemodel import simpletreemodel_rc
import utilities.readncnr3 as readncnr
import numpy as N
import utilities.scriptutil as SU
import re
from utilities.simple_combine import simple_combine
import copy
import pylab
from utilities.findpeak3 import findpeak
from openopt import NLP
import scipy.optimize
import scipy.odr
from scipy.optimize import curve_fit
pi=N.pi
from mpfit import mpfit
import rescalculator.rescalc as rescalc


nodetypes=set(['hkl','th','tth','q','other','leaf'])
#for hkl nodes, the itemdata is a string representation of hkl
#for other branches, the itemdata is a string of the branch type
#for the leaves, itemdata is the filename, measured_data will contain the actual data
#associated with the measurement

class TreeItem(object):
    def __init__(self, data, parent=None,nodetype='hkl',measured_data=None):
        self.parentItem = parent
        self.itemData = data
        self.childItems = []
        self.nodetype=nodetype
        self.measured_data=measured_data
        self._checkState=QtCore.Qt.Unchecked
        self.q=q
        self.Q=0.0
        self.mon0=1.0
        self.th_correction=1.0
        self.tth_correction=1.0
        self.q_correction=1.0
        self.th_integrated_intensity=0.0
        self.tth_integrated_intensity=0.0
        self.q_integrated_intensity=0.0
        

    def checkState(self):
        return self._checkState
    def setcheckState(self,checkState):
        self._checkState=checkState
    
    def toggleCheck(self):
        if self._checkState==QtCore.Qt.Unchecked:
            self._checkState=QtCore.Qt.Checked
        elif self._checkState==QtCore.Qt.Checked:
            self._checkState=QtCore.Qt.Unchecked
            

    def appendChild(self, item):
        self.childItems.append(item)

    def child(self, row):
        return self.childItems[row]

    def childCount(self):
        return len(self.childItems)

    def columnCount(self):
        return len(self.itemData)

    def data(self, column):
        return self.itemData[column]

    def parent(self):
        return self.parentItem

    def row(self):
        if self.parentItem:
            return self.parentItem.childItems.index(self)

        return 0


class TreeModel(QtCore.QAbstractItemModel):
    def __init__(self, filestrlist, parent=None,mon0=1.0):
        QtCore.QAbstractItemModel.__init__(self, parent)
        self.mon0=mon0
        self.idMap = {}
        self.hklmap={}
        rootData = []
        rootData.append(QtCore.QVariant("HKL"))
        rootData.append(QtCore.QVariant("Summary"))
        self.rootItem = TreeItem(rootData)
        self.idMap[id(self.rootItem)] = self.rootItem
        self.setupModelData(filestrlist, self.rootItem)

    def columnCount(self, parent):
        if parent.isValid():
            return self.idMap[parent.internalId()].columnCount()
        else:
            return self.rootItem.columnCount()

    def data(self, index, role):
        if not index.isValid():
            return QtCore.QVariant()

        if index.column()==0 and role==QtCore.Qt.CheckStateRole:
            item = self.idMap[index.internalId()]
            print 'checkstate', item.checkState()
            return QtCore.QVariant(item.checkState())
        if role != QtCore.Qt.DisplayRole:
            return QtCore.QVariant()

        try:
            item = self.idMap[index.internalId()]
            return QtCore.QVariant(item.data(index.column()))
        except KeyError:
            return QtCore.QVariant()
        
    def setData(self, index, value, role=QtCore.Qt.EditRole):
        if (role == QtCore.Qt.CheckStateRole and index.column() == 0):
            #self.checkstates[self.fileInfo(index).absoluteFilePath()] = QtCore.Qt.CheckState() 
            item = self.idMap[index.internalId()]
            #item.setcheckState(QtCore.Qt.CheckState()) 
            item.toggleCheck()
            print 'setting data',QtCore.Qt.CheckState()
            self.emit(QtCore.SIGNAL("dataChanged(QtCore.QModelIndex,QModelIndex)"), index, index)
            return True

        return QtCore.QAbstractItemModel.setData(self, index, value, role)

    def flags(self, index):
        if not index.isValid():
            return QtCore.Qt.ItemIsEnabled

        return QtCore.Qt.ItemIsEnabled | QtCore.Qt.ItemIsSelectable| QtCore.Qt.ItemIsUserCheckable

    def headerData(self, section, orientation, role):
        if orientation == QtCore.Qt.Horizontal and role == QtCore.Qt.DisplayRole:
            return self.rootItem.data(section)

        return QtCore.QVariant()

    def index(self, row, column, parent):
        if row < 0 or column < 0 or row >= self.rowCount(parent) or column >= self.columnCount(parent):
            return QtCore.QModelIndex()

        if not parent.isValid():
            parentItem = self.rootItem
        else:
            parentItem = self.idMap[parent.internalId()]

        childItem = parentItem.child(row)
        if childItem:
            index = self.createIndex(row, column, id(childItem))
            self.idMap.setdefault(index.internalId(), childItem)
            return index
        else:
            return QtCore.QModelIndex()

    def parent(self, index):
        if not index.isValid():
            return QtCore.QModelIndex()

        try:
            childItem = self.idMap[index.internalId()]
            parentItem = childItem.parent()

            if parentItem == self.rootItem:
                return QtCore.QModelIndex()

            return self.createIndex(parentItem.row(), 0, id(parentItem))
        except KeyError:
            return QtCore.QModelIndex()

    def rowCount(self, parent):
        if parent.column() > 0:
            return 0

        try:
            if not parent.isValid():
                parentItem = self.rootItem
            else:
                parentItem = self.idMap[parent.internalId()]

            return parentItem.childCount()
        except:
            return 0

    def addnode(self,data,parent,nodetype=None,measured_data=None):
        node=TreeItem(data, parent,nodetype=nodetype, measured_data=measured_data)
        if len(self.hklmap.keys())==0:
            node.mon0=mydata.metadata['count_info']['monitor']
        else:
            
        self.idMap[id(node)] = node
        parent.appendChild(node)
        return node
    
    
    
    def place_data(self,mydata,tol=1e-6):
        if mydata.metadata['file_info']['scantype']=='b':
            #print 'b'
            currfile=mydata.metadata['file_info']['filename']
            if N.abs(mydata.metadata['motor4']['step'])<tol and N.abs(mydata.metadata['motor3']['step'])>tol:
                #print currfile, 'a3 scan'
                target='th'
                #self.th.append(data_item(mydata))
                #print 'self.th',self.th
            elif N.abs(mydata.metadata['motor4']['step']-2*mydata.metadata['motor3']['step'])<tol and N.abs(mydata.metadata['motor3']['step'])>tol:
                #print currfile, 'th-2th scan'
                #self.th2th.append(data_item(mydata))
                target='tth'
            else:
                #print currfile, 'strange scan'
                #self.other.append(data_item(mydata))
                target='other'
        return target

    def setupModelData(self, filestrlist, parent):
        parents = []
        indentations = []
        parents.append(parent)
        indentations.append(0)

        #myfilestrlist=[r'C:\Ce2RhIn8\Mar10_2009\magsc035.bt9',r'C:\Ce2RhIn8\Mar10_2009\magsc034.bt9']
        
        
        
        for myfilestr in filestrlist:
            mydatareader=readncnr.datareader()
            mydata=mydatareader.readbuffer(myfilestr)
            filename=mydata.metadata['file_info']['filename']
            h=str(mydata.metadata['q_center']['h_center'])
            k=str(mydata.metadata['q_center']['k_center'])
            l=str(mydata.metadata['q_center']['l_center'])
            hkl=h+k+l
    
            print 'hkl',hkl
            #hkl=QtCore.QString(hkl)
    
            #nodetypes=set(['hkl','th','tth','q','other','leaf'])
    
    
            if hkl not in self.hklmap.keys():           
                hkl_data=[hkl,'']
                hklnode=self.addnode(hkl_data,parents[-1],nodetype='hkl')
                self.hklmap[hkl]=id(hklnode)
                #add the branches
                data=['theta','']
                thnode=self.addnode(data,hklnode,nodetype='th')
                data=['ttheta','']
                tthnode=self.addnode(data,hklnode,nodetype='tth')
                data=['q','']
                qnode=self.addnode(data,hklnode,nodetype='q')
                data=['other','']
                othernode=self.addnode(data,hklnode,nodetype='other')
                data=[filename,'']
                targetdict={}
                targetdict['th']=thnode
                targetdict['tth']=tthnode
                targetdict['other']=othernode
                targetdict['qnode']=qnode
                targetnode=targetdict[self.place_data(mydata)]
                leaf=self.addnode(data,targetnode,nodetype='leaf',measured_data=mydata)
                idx=self.index(0,0,QtCore.QModelIndex())
                idx.model().setData(idx,QtCore.QVariant(QtCore.Qt.Checked), QtCore.Qt.CheckStateRole) 


        number = 0
        if 0:
            while number < len(lines):
                position = 0
                while position < len(lines[number]):
                    if lines[number][position] != " ":
                        break
                    position += 1

                lineData = lines[number][position:].trimmed()

                if not lineData.isEmpty():
                    # Read the column data from the rest of the line.
                    columnStrings = lineData.split("\t", QtCore.QString.SkipEmptyParts)
                    columnData = []
                    for column in range(0, len(columnStrings)):
                        columnData.append(columnStrings[column])

                    if position > indentations[-1]:
                        # The last child of the current parent is now the new parent
                        # unless the current parent has no children.

                        if parents[-1].childCount() > 0:
                            parents.append(parents[-1].child(parents[-1].childCount() - 1))
                            indentations.append(position)

                    else:
                        while position < indentations[-1] and len(parents) > 0:
                            parents.pop()
                            indentations.pop()

                    # Append a new item to the current parent's list of children.
                    item = TreeItem(columnData, parents[-1])
                    self.idMap[id(item)] = item
                    parents[-1].appendChild(item)

                number += 1


class myTreeView(QtGui.QTreeView):
    def __init__(self, parent=None):
        super(myTreeView, self).__init__(parent)

        #self.myModel = myModel()
        #self.setModel(self.myModel)
        #item=self.currentItem()
        #item.setCheckState(0, Qt.Unchecked) # 0 is the column number
        
        #f = QtCore.QFile(":/default.txt")
        #f.open(QtCore.QIODevice.ReadOnly)
        #self.myModel=TreeModel(QtCore.QString(f.readAll()))
        #f.close()
        filestrlist=[r'C:\Ce2RhIn8\Mar10_2009\magsc035.bt9',r'C:\Ce2RhIn8\Mar10_2009\magsc034.bt9']
        self.myModel=TreeModel(filestrlist)
        self.setModel(self.myModel)
        self.dragEnabled()
        self.acceptDrops()
        self.showDropIndicator()
        self.setDragDropMode(QtGui.QAbstractItemView.InternalMove)
        self.connect(self.model(), QtCore.SIGNAL("dataChanged(QtCore.QModelIndex,QtCore.QModelIndex)"), self.change)
        self.expandAll()
        #QtCore.QObject.connect(self.selectionModel(),QtCore.SIGNAL("selectionChanged(QItemSelection, QItemSelection)"),
        # self.itemselected)
        QtCore.QObject.connect(self.selectionModel(),QtCore.SIGNAL("currentChanged(QModelIndex, QModelIndex)"),
         self.currentselected)
        

    def change(self, topLeftIndex, bottomRightIndex):
        self.update(topLeftIndex)
        self.expandAll()
        self.expanded()
    def expanded(self):
        for column in range(self.model().columnCount(QModelIndex())):
            self.resizeColumnToContents(column)               
    def itemselected(self,selected, deselected):
        print 'itemselected'
        print len(selected), "items selected"
        print len(deselected), "items deselected"
        print 'indexes', selected[0].indexes()
        idx=selected[0].indexes()[0]
        model=selected[0].model()
        print 'item',model.idMap[idx.internalId()].itemData
        idx2=selected[0].indexes()[1]
        print 'item2',model.idMap[idx2.internalId()].itemData
        #Why are there two???
    def currentselected(self,selected,deselected):
        print 'currentselected'
        print 'selected',selected.model().idMap[selected.internalId()].itemData
        if not deselected.model()==None:
            print 'deselected',deselected.model().idMap[deselected.internalId()].itemData
        

if __name__ == "__main__":
    app = QtGui.QApplication(sys.argv)

    #f = QtCore.QFile(":/default.txt")
    #f.open(QtCore.QIODevice.ReadOnly)
    #model = TreeModel(QtCore.QString(f.readAll()))
    #f.close()

    #view = QtGui.QTreeView()
    view=myTreeView()
    #view.setModel(model)
    view.setWindowTitle("Simple Tree Model")
    view.show()
    sys.exit(app.exec_())
