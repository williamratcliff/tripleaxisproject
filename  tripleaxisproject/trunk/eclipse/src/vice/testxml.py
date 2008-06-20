import sys,os
from xml.dom import minidom, Node


class SimpleScanner:
    def __init__(self,doc):
        root=doc.childNodes[0]
        print 'root',root
        for child in root.childNodes:
            if child.nodeType==Node.ELEMENT_NODE and \
            child.tagName=='menu':
                self.handleMenu(child)
                
    def gettext(self, nodelist):
        """Given a list of one or more nodes, revursively finds all text nodes
        in that list, concats them, returns the result."""
        retlist=[]
        return retlist
        
        
                
    def handleMenu(self,node):
        """Process the menu tag  """
        
        for child in node.childNodes:
            if child.nodeType!=Node.ELEMENT_NODE:
                continue
            if child.tagName=='title':
                #print child.__dict__
                text=child.childNodes[0].data
                #print text.data
                print 'menu title is', text# self.gettext(child.childNodes)
            if child.tagName=='submenu':
                self.handleSubMenu(child)
                continue
            if child.tagName=='variables':
                text=child.childNodes[0].data
                split_text=text.lower().encode('ascii').replace('\n','').replace(' ','').split(',')
                print 'variables',split_text


    def handleSubMenu(self,node):
        """Process the menu tag  """
        
        for child in node.childNodes:
            if child.nodeType!=Node.ELEMENT_NODE:
                continue
            if child.tagName=='title':
                #print child.__dict__
                text=child.childNodes[0].data
                #print text.data
                print 'submenu title is', text# self.gettext(child.childNodes)
            if child.tagName=='variables':
                text=child.childNodes[0].data
                split_text=text.lower().encode('ascii').replace('\n','').replace(' ','').split(',')
                print 'submenu variables',split_text
            



if __name__=='__main__':
    mydirectory=r'C:\tripleaxisproject2\trunk\eclipse\src\vice'
    #mydirectory=r'C:\mytripleaxisproject\trunk\eclipse\src\vice'
    xmlfile='test2.xml'
    myfilestr=os.path.join(mydirectory,xmlfile)
    #print myfilestr
    doc=minidom.parse(myfilestr)
    myscanner=SimpleScanner(doc)
    #print doc.childNodes