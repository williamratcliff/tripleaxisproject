from bbfreeze import Freezer
f = Freezer("polcorrecter", includes=("writebt7","readncnr2"))
f.addScript("polcorrect2.py")
f.include_py=False
#f.addScript("hello-version.py")
f()    # starts the freezing process