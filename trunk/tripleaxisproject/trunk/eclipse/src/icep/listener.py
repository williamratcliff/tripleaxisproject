from ice.event.communication import DeviceValueChangedListener
from ice.event.communication import DeviceValueChangedEvent
from ice.jython.management import DeviceListenerProxy

class A3Listener(DeviceValueChangedListener):
	values = ""

	def attachToA3(self):
		controller = Controller.getReference()
		a3Device = controller.getFirstDevice("A3")
		a3Device.addValueChangeListener(self)

	#
	# please remember to remove the listener after you are done with it
	# you are only aloud to have so many.
	#
	def detachFromA3(self):
		controller = Controller.getReference()
		a3Device = controller.getFirstDevice("A3")
		a3Device.removeValueChangeListener(self)

	#
	# don't complain I already am not happy. For some reason
	# print does not work from inside a listener call back. Not sure
	# why I tried System.out.print and it threw an exception. This should
	# make debuging challenging. 
	# I had to concat the values on a string to prove it works for you. I already
	# spent six hours of my day figuring out that it was working, but I just couldn't print
	#
	def actionPerformed(self,dvce):
		theProperty = dvce.getPropertyName()
		o = dvce.getOldValue()
		n = dvce.getNewValue()
		if theProperty == "hardwareValue":
			self.values += " oldValue=%s,newValue=%s" % (o,n)

	def getValues(self):
		return self.values
		
		

class test():

	def __init__(self):
		print "test"

	def moveA3(self,a3Value):
		a3Listener = A3Listener()

		a3Listener.attachToA3()
		m = MoveCommand("A3",a3Value)
		m.run()
		a3Listener.detachFromA3()
		print a3Listener.getValues()
