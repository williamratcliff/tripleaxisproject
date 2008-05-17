from bbfreeze import Freezer
f = Freezer("polcorrector", includes=("numpy",))
f.addScript("polcorrect2.py")
#f.addScript("hello-version.py")
f()    # starts the freezing process
print 'done'