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
import pylab


def check_q(q1,q2,tol=1e-6):
    heq=False
    keq=False
    leq=False
    #print 'q1,q2',q1,q2
    if N.abs(q2['h_center']-q1['h_center'])< tol:
        heq=True
    if N.abs(q2['k_center']-q1['k_center'])< tol:
        keq=True
    if N.abs(q2['l_center']-q1['l_center'])< tol:
        leq=True
    #print 'heq',heq,'keq',keq,'leq',leq

    return (heq and keq and leq)



nodetypes=set(['hkl','th','tth','q','other','leaf'])
#for hkl nodes, the itemdata is a string representation of hkl
#for other branches, the itemdata is a string of the branch type
#for the leaves, itemdata is the filename, measured_data will contain the actual data
#associated with the measurement

class TreeItem(object):
    def __init__(self, data, parent=None,nodetype='hkl',measured_data=None,q=None):
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
            #print 'checkstate', item.checkState()
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
        
    def correct_data(self,mydata,qscan=None):
        mya=mydata.metadata['lattice']['a']
        myb=mydata.metadata['lattice']['b']
        myc=mydata.metadata['lattice']['c']
        myalpha=N.radians(mydata.metadata['lattice']['alpha'])
        mybeta=N.radians(mydata.metadata['lattice']['beta'])
        mygamma=N.radians(mydata.metadata['lattice']['gamma'])
        a=N.array([mya],dtype='float64') #for now, the code is broken if only one element in the array for indexing
        b=N.array([myb],dtype='float64')
        c=N.array([myc],dtype='float64')
        alpha=N.array([myalpha],dtype='float64')
        beta=N.array([mybeta],dtype='float64')
        gamma=N.array([mygamma],dtype='float64')
        #print mydata.metadata
        h=mydata.metadata['orient1']['h']
        k=mydata.metadata['orient1']['k']
        l=mydata.metadata['orient1']['l']
        orient1=N.array([[h,k,l]],dtype='float64')
        h=mydata.metadata['orient2']['h']
        k=mydata.metadata['orient2']['k']
        l=mydata.metadata['orient2']['l']
        orient2=N.array([[h,k,l]],dtype='float64')
        orientation=rescalc.lattice_calculator.Orientation(orient1,orient2)
        mylattice=rescalc.lattice_calculator.Lattice(a=a,b=b,c=c,alpha=alpha,beta=beta,gamma=gamma,orientation=orientation)
        #h=q['h_center'] #perhaps just take this from the file in the future
        #k=q['k_center']
        #l=q['l_center']
        h=mydata.metadata['q_center']['h_center']
        k=mydata.metadata['q_center']['k_center']
        l=mydata.metadata['q_center']['l_center']
        H=N.array([h,h],dtype='float64');
        K=N.array([k,k],dtype='float64');
        L=N.array([l,l],dtype='float64');
        W=N.array([0.0,0.0],dtype='float64')
        EXP={}
        EXP['ana']={}
        EXP['ana']['tau']='pg(002)'
        EXP['mono']={}
        EXP['mono']['tau']='pg(002)';
        EXP['ana']['mosaic']=25
        EXP['mono']['mosaic']=25
        EXP['sample']={}
        EXP['sample']['mosaic']=25
        EXP['sample']['vmosaic']=25
        coll1=mydata.metadata['collimations']['coll1']
        coll2=mydata.metadata['collimations']['coll2']
        coll3=mydata.metadata['collimations']['coll3']
        coll4=mydata.metadata['collimations']['coll4']	
        EXP['hcol']=N.array([coll1,coll2,coll3,coll4],dtype='float64')
        EXP['vcol']=N.array([120, 120, 120, 240],dtype='float64')

        EXP['infix']=-1 #positive for fixed incident energy
        EXP['efixed']=mydata.metadata['energy_info']['ef']
        EXP['method']=0
        setup=[EXP]
        myrescal=rescalc.rescalculator(mylattice)
        newinput=rescalc.lattice_calculator.CleanArgs(a=a,b=b,c=c,alpha=alpha,beta=beta,gamma=gamma,orient1=orient1,orient2=orient2,\
                                                      H=H,K=K,L=L,W=W,setup=setup)
        neworientation=rescalc.lattice_calculator.Orientation(newinput['orient1'],newinput['orient2'])
        mylattice=rescalc.lattice_calculator.Lattice(a=newinput['a'],b=newinput['b'],c=newinput['c'],alpha=newinput['alpha'],
                                                     beta=newinput['beta'],gamma=newinput['gamma'],orientation=neworientation)
        myrescal.__init__(mylattice)
        Q=myrescal.lattice_calculator.modvec(H,K,L,'latticestar')
        #print 'Q', Q
        R0,RM=myrescal.ResMat(Q,W,setup)
        #print 'RM '
        #print RM.transpose()
        #print 'R0 ',R0
        #exit()
        R0,RMS=myrescal.ResMatS(H,K,L,W,setup)
        #myrescal.ResPlot(H, K, L, W, setup)
        #print 'RMS'
        #print RMS.transpose()[0]
        #corrections=myrescal.calc_correction(H,K,L,W,setup,qscan=[[1,0,1],[1,0,1]])
        corrections=myrescal.calc_correction(H,K,L,W,setup,qscan=qscan)
        #print corrections
        return corrections,Q[0]
        
    def correct_node(self,qnode,nodetype='leaf'):
        #qnode=self.qlist[index]
        #for now, handle th, 2th corrections, otherwise, need to access corrections['q_correction']
        #we would also have to send in qscan as an array
        #Also, for now, we will have corrections for every scan, we should change this
        #to corrections for each h,k,l node
        if nodetype=='leaf':
            mydata=qnode.measured_data
            corrections,Q=self.correct_data(mydata,qscan=None)
            th_correction=corrections['th_correction'][0]
            qnode.th_correction=th_correction.flatten()[0]
            tth_correction=corrections['tth_correction'][0]
            qnode.tth_correction=tth_correction.flatten()[0]
            qnode.Q=Q
            print 'corrected', qnode.th_correction


        return
    
    def condense_node(self,qnode):
        print qnode.q
        #print qnode.th

        a3=[]
        counts=[]
        counts_err=[]
        monlist=[]
        for mydataitem in qnode.th:
            mydata=mydataitem.data
            monlist.append(mydata.metadata['count_info']['monitor'])
            counts_err.append(N.array(mydata.data['counts_err']))
            counts.append(N.array(mydata.data['counts']))
            a3.append(N.array(mydata.data['a3']))
        a3_out,counts_out,counts_err_out=simple_combine(a3,counts,counts_err,monlist)

        #print a3_out.shape
        #print counts_out.shape
        #print counts_err_out.shape
        qnode.th_condensed={}
        qnode.th_condensed['a3']=a3_out
        qnode.th_condensed['counts']=counts_out
        qnode.th_condensed['counts_err']=counts_err_out

        print qnode.th_condensed['counts'].std()
        print qnode.th_condensed['counts'].mean()
        print qnode.th_condensed['counts'].max()
        print qnode.th_condensed['counts'].min()
        if 0:
            pylab.errorbar(a3_out,counts_out,counts_err_out,marker='s',linestyle='None',mfc='black',mec='black',ecolor='black')
            pylab.show()       
        return 

    def addnode(self,data,parent,nodetype=None,measured_data=None):
        node=TreeItem(data, parent,nodetype=nodetype, measured_data=copy.deepcopy(measured_data))
        if (not measured_data==None) and nodetype=='leaf':
            print 'adding leaf len',len(self.hklmap.keys())
            if len(self.hklmap.keys())==1:
                print 'first leaf'
                self.mon0=measured_data.metadata['count_info']['monitor']
        
            mon=node.measured_data.metadata['count_info']['monitor']
            node.measured_data.data['counts']=N.array(node.measured_data.data['counts'],'Float64')
            #node.measured_data.data['counts_err']=N.array(node.measured_data.data['counts_err'],'Float64')
            counts_new=node.measured_data.data['counts']*self.mon0/mon
            counts_new_err=N.sqrt(node.measured_data.data['counts'])*self.mon0/mon
            node.measured_data.data['counts']=counts_new
            node.measured_data.data['counts_err']=counts_new_err
            node.measured_data.metadata['count_info']['monitor']=self.mon0
            node.mon0=self.mon0
            print 'mon0 in add self',self.mon0,'node',node.mon0
            node.q=copy.copy(measured_data.metadata['q_center'])
            self.correct_node(node)
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
            hkl='('+h+' '+k+' '+l+')'
    
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
        self.resizeColumnToContents(0)
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
        currselected=selected.model().idMap[selected.internalId()].measured_data
        print 'selected',selected.model().idMap[selected.internalId()].itemData
        node=selected.model().idMap[selected.internalId()]
        print 'nodetype',node.nodetype
        if node.nodetype=='leaf':
            measured_data=node.measured_data
            print 'leaf'
            #scantype=node.parentItem.itemData[0]
            scantype=node.parentItem.nodetype
            print 'scantype',scantype
            if scantype in ['tth','th','q']:
                plotdict={}
                print 'valid scantype'
                plotdict['scantype']=scantype
                plotdict['data']={}
                if scantype=='th':
                    plotdict['data']['x']=node.measured_data.data['a3']
                    plotdict['xlabel']='th (degrees)'
                plotdict['data']['y']=node.measured_data.data['counts']
                plotdict['data']['yerr']=node.measured_data.data['counts_err']
                plotdict['ylabel']='Counts/%g5.1'%(node.mon0)  #assume counting by neutrons for now    
                self.emit(QtCore.SIGNAL("clearplot"),'clear')
                self.emit(QtCore.SIGNAL("plot"),plotdict)
                print 'emitted signal'
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
