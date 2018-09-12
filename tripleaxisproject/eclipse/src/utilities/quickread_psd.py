from . import readncnr3 as readncnr


myfilestr='c:\psd_test\psd72615.bt7'
mydatareader=readncnr.datareader()
mydata=mydatareader.readbuffer(myfilestr,lines=91)
print(list(mydata.data.keys()))
print(mydata.data['psdc11'])