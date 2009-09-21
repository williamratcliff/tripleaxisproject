import sys
from PyQt4 import QtCore, QtGui

app = QtGui.QApplication(sys.argv)

model = QtGui.QDirModel()
tree = QtGui.QTreeView()
tree.setModel(model)

tree.setWindowTitle(tree.tr("Dir View"))
tree.resize(640, 480)
tree.show()

sys.exit(app.exec_())