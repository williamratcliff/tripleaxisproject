import sys
from PyQt4 import QtCore, QtGui


class TreeModel(QtCore.QAbstractItemModel):
    NAME = 0
    FILEID = QtCore.Qt.UserRole + 1
    horizontalHeaderLabels = ["File Name",]
    inventory = None

    def set_tree(self, inventory, root_item):
        self.emit(QtCore.SIGNAL("layoutAboutToBeChanged()"))
        self.inventory = inventory
        self.id2fileid = []
        self.fileid2id = {}
        self.dir_children_ids = {}
        self.parent_ids = []

        # Create internal ids for all items in the tree for use in
        # ModelIndex's.
        root_fileid = root_item.file_id
        self.append_fileid(root_fileid, None)
        remaining_dirs = [root_fileid,]
        while remaining_dirs:
            dir_fileid = remaining_dirs.pop(0)
            dir_id = self.fileid2id[dir_fileid]
            dir_children_ids = []
            for child in inventory[dir_fileid].children:
                id = self.append_fileid(child.file_id, dir_id)
                dir_children_ids.append(id)
                if child.children:
                    remaining_dirs.append(child.file_id)

                if len(self.id2fileid) % 100 == 0:
                    QtCore.QCoreApplication.processEvents()
            self.dir_children_ids[dir_id] = dir_children_ids

        self.emit(QtCore.SIGNAL("layoutChanged()"))

    def append_fileid(self, fileid, parent_id):
        ix = len(self.id2fileid)
        self.id2fileid.append(fileid)
        self.parent_ids.append(parent_id)
        self.fileid2id[fileid] = ix
        return ix

    def columnCount(self, parent):
        if parent.isValid():
            return 0
        return len(self.horizontalHeaderLabels)

    def rowCount(self, parent):
        if self.inventory is None:
            return 0
        parent_id = parent.internalId()
        if parent_id not in self.dir_children_ids:
            return 0
        return len(self.dir_children_ids[parent_id])

    def _index(self, row, column, parent_id):
        item_id = self.dir_children_ids[parent_id][row]
        return self.createIndex(row, column, item_id)

    def index(self, row, column, parent = QtCore.QModelIndex()):
        if self.inventory is None:
            return self.createIndex(row, column, 0)
        parent_id = parent.internalId()
        return self._index(row, column, parent_id)

    def sibling(self, row, column, index):
        sibling_id = child.internalId()
        if sibling_id == 0:
            return QtCore.QModelIndex()
        parent_id = self.parent_ids[child_id]
        return self._index(row, column, parent_id)

    def parent(self, child):
        child_id = child.internalId()
        if child_id == 0:
            return QtCore.QModelIndex()
        item_id = self.parent_ids[child_id]
        if item_id == 0 :
            return self.createIndex(0, 0, item_id)

        parent_id = self.parent_ids[item_id]
        row = self.dir_children_ids[parent_id].index(item_id)
        return QtCore.QModelIndex()
        #return self.createIndex(row, 0, item_id)

    def hasChildren(self, parent):
        if self.inventory is None:
            return False

        parent_id = parent.internalId()
        return parent_id in self.dir_children_ids

    def data(self, index, role):
        if not index.isValid():
            return QtCore.QVariant()

        fileid = self.id2fileid[index.internalId()]

        if role == self.FILEID:
            return QtCore.QVariant(fileid)

        item = self.inventory[fileid]

        column = index.column()
        if column == self.NAME:
            if role == QtCore.Qt.DisplayRole:
                return QtCore.QVariant(item.file_name)

        return QtCore.QVariant()

    def flags(self, index):
        if not index.isValid():
            return QtCore.Qt.ItemIsEnabled

        return QtCore.Qt.ItemIsEnabled | QtCore.Qt.ItemIsSelectable

    def headerData(self, section, orientation, role):
        if orientation == QtCore.Qt.Horizontal and role == QtCore.Qt.DisplayRole:
            return QtCore.QVariant(self.horizontalHeaderLabels[section])
        return QtCore.QVariant()

inventory = {}

class InventoryItem():

    def __init__(self, file_id, file_name, children=[]):
        self.file_id = file_id
        self.file_name = file_name
        self.children = children
        global inventory
        inventory[file_id] = self

root_item = InventoryItem("root-id", "", [
    InventoryItem("dir1-id", "dir1", [
        InventoryItem("file1-id", "file1")
    ]),
    InventoryItem("file1-id", "file1")
])

app = QtGui.QApplication(sys.argv)

model = TreeModel()
model.set_tree(inventory, root_item)

tree_view = QtGui.QTreeView()
tree_view.setModel(model)

tree_view.show()
app.exec_()