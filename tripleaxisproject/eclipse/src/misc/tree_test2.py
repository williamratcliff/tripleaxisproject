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


class TreeItem(object):
    def __init__(self, data, parent=None):
        self.parentItem = parent
        self.itemData = data
        self.childItems = []

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
    def __init__(self, data, parent=None):
        QtCore.QAbstractItemModel.__init__(self, parent)

        self.idMap = {}

        rootData = []
        rootData.append(QtCore.QVariant("Title"))
        rootData.append(QtCore.QVariant("Summary"))
        self.rootItem = TreeItem(rootData)
        self.idMap[id(self.rootItem)] = self.rootItem
        self.setupModelData(data.split("\n"), self.rootItem)

    def columnCount(self, parent):
        if parent.isValid():
            return self.idMap[parent.internalId()].columnCount()
        else:
            return self.rootItem.columnCount()

    def data(self, index, role):
        if not index.isValid():
            return QtCore.QVariant()

        if role != QtCore.Qt.DisplayRole:
            return QtCore.QVariant()

        try:
            item = self.idMap[index.internalId()]
            return QtCore.QVariant(item.data(index.column()))
        except KeyError:
            return QtCore.QVariant()

    def flags(self, index):
        if not index.isValid():
            return QtCore.Qt.ItemIsEnabled

        return QtCore.Qt.ItemIsEnabled | QtCore.Qt.ItemIsSelectable

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

    def setupModelData(self, lines, parent):
        parents = []
        indentations = []
        parents.append(parent)
        indentations.append(0)

        number = 0

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


if __name__ == "__main__":
    app = QtGui.QApplication(sys.argv)

    f = QtCore.QFile(":/default.txt")
    f.open(QtCore.QIODevice.ReadOnly)
    model = TreeModel(QtCore.QString(f.readAll()))
    f.close()

    view = QtGui.QTreeView()
    view.setModel(model)
    view.setWindowTitle("Simple Tree Model")
    view.show()
    sys.exit(app.exec_())