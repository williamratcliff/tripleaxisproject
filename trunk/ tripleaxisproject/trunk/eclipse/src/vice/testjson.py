import demjson
import sys,os


if __name__=='__main__':
    mydirectory=r'C:\tripleaxisproject2\trunk\eclipse\src\vice'
    mydirectory=r'C:\mytripleaxisproject\trunk\eclipse\src\vice'
    jsonfile='test.json'
    myfilestr=os.path.join(mydirectory,jsonfile)
    #print myfilestr
    myfile=open(myfilestr)
    jsonstr=myfile.read()
    #print jsonstr
    json_obj=demjson.decode(jsonstr)
    print json_obj['menu']
    print json_obj['menu'][-1]
