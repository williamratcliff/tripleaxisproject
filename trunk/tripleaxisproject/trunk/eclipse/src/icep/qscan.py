class test:
	def run(self):
		controller = Controller.getReference()
		experiment = controller.getExperiment()

		exprId = experiment.getID()
		qlist = [[1,0,3],[2,0,4],[5,0,4]]
		qscan = QScan(qlist,exprId,"Detector","Time",0,1,10,-1.0,1.0)
		qscan.scan()


class QScan:
	qList = [[]]
	baseName = "A3_HKL_MOTORSCAN_"
	exprId = ""
	detectorType = "Detector"
	countType = "Time"
	prefactor = 1
	counts = 1
	numberOfPoints = 1
	rangeStart = 0
	rangeStop = 1
	

	def __init__(self,qList,exprId,detectorType,countType,prefactor,counts,numberOfPoints,rangeStart,rangeStop):
		self.setQList(qList)
		self.setExprId(exprId)
		self.setDetectorType(detectorType)
		self.setCountType(countType)
		self.setPrefactor(prefactor)
		self.setCounts(counts)
		self.setNumberOfPoints(numberOfPoints)
		self.setRange(rangeStart,rangeStop)

	def getQList(self):
		return self.qList

	def setQList(self,qList):
		self.qList = qList

	def getBaseName(self):
		return self.baseName

	def setBaseName(self,baseName):
		self.baseName = baseName

	def getExprId(self):
		return self.exprId

	def setExprId(self,exprId):
		self.exprId = exprId 

	def getDetectorType(self):
		return self.detectorType

	def setDetectorType(self,detectorType):
		self.detectorType = detectorType

	def getCountType(self):
		return self.countType

	def setCountType(self,countType):
		self.countType = countType

	def getPrefactor(self):
		return self.prefactor

	def setPrefactor(self,prefactor):
		self.prefactor = prefactor

	def getCounts(self):
		return self.counts

	def setCounts(self,counts):
		self.counts = counts

	def getNumberOfPoints(self):
		return self.numberOfPoints
	
	def setNumberOfPoints(self,numberOfPoints):
		self.numberOfPoints = numberOfPoints

	def setRangeStart(self,rangeStart):
		self.rangeStart = rangeStart

	def setRangeStop(self,rangeStop):
		self.rangeStop = rangeStop

	def setRange(self,rangeStart,rangeStop):
		self.setRangeStart(rangeStart)
		self.setRangeStop(rangeStop)
	
	def getRangeStart(self):
		return self.rangeStart

	def getRangeStop(self):
		return self.rangeStop

	def scan(self):
		baseName = self.getBaseName()
		qList = self.getQList()
		exprId = self.getExprId()
		detectorType = self.getDetectorType()
		countType = self.getCountType()
		prefactor = self.getPrefactor()
		counts = self.getCounts()
		numberOfPoints = self.getNumberOfPoints()
		rangeStart = self.getRangeStart()
		rangeStop = self.getRangeStop()

		mScan = MotorScan("A3",baseName,exprId,numberOfPoints,counts,prefactor,detectorType,countType,rangeStart,rangeStop)
		for q in qList:
			newName = baseName + "%i-%i-%i" % (q[0],q[1],q[2])
			hklDestination = "[%f,%f,%f]" % (q[0],q[1],q[2])  
			moveCommand = MoveCommand("HKL",hklDestination)
			moveCommand.run()
			mScan.setName(newName)
			mScan.runScan(True)
			

class MotorScan:
	numberOfPoints = 0
	counts = 0;
	prefactor = 0;
	detectorType = "Detector"
	countType = "Time"
	rangeStart = 0.0
	rangeStop = 0.0
	motor = ""
	name = ""
	exptID = ""

	def __init__(self):
		print "motor scan created"

	def __init__(self,motor,name,exptId,numberOfPoints,counts,prefactor,detectorType,countType,rangeStart,rangeStop):
		self.initializeScan(motor,name,exptId,numberOfPoints,counts,prefactor,detectorType,countType,rangeStart,rangeStop)


	def setName(self,name):
		self.name = name

	def getName(self):
		return self.name

	def setMotor(self,motor):
		self.motor = motor

	def getMotor(self):
		return self.motor

	def getExptId(self):
		return self.exptId

	def setExptId(self,exptId):
		self.exptId = exptId

	def getNumberOfPoints(self):
		return self.numberOfPoints
	
	def setNumberOfPoints(self,numberOfPoints):
		self.numberOfPoints = numberOfPoints

	def getCounts(self):
		return self.counts

	def setCounts(self,counts):
		self.counts = counts

	def getPrefactor(self):
		return self.prefactor

	def setPrefactor(self,prefactor):
		self.prefactor = prefactor

	def getDetectorType(self):
		return self.detectorType

	def setDetectorType(self,detectorType):
		self.detectorType = detectorType

	def setCountType(self,countType):
		self.countType = countType

	def getCountType(self):
		return self.countType

	def setRangeStart(self,rangeStart):
		self.rangeStart = rangeStart

	def setRangeStop(self,rangeStop):
		self.rangeStop = rangeStop

	def setRange(self,rangeStart,rangeStop):
		self.setRangeStart(rangeStart)
		self.setRangeStop(rangeStop)
	
	def getRangeStart(self):
		return self.rangeStart

	def getRangeStop(self):
		return self.rangeStop

	def deleteScanFromScanList(self):
		controller = Controller.getReference()
		sendManager = controller.getSendManager()
		scanName = self.getName()
		deleteMessage = "scan delete %s" % (scanName)

		sendManager.addMessage(deleteMessage)

	def initializeScan(self,motor,name,exptId,numberOfPoints,counts,prefactor,detectorType,countType,rangeStart,rangeStop):
		self.setMotor(motor)
		self.setName(name)
		self.setExptId(exptId)
		self.setNumberOfPoints(numberOfPoints)
		self.setCounts(counts)
		self.setPrefactor(prefactor)
		self.setDetectorType(detectorType)
		self.setCountType(countType)
		self.setRange(rangeStart,rangeStop)

	def buildScanDescription(self):
		description = ""
		description += "Scan:SubID=%s:" % self.getExptId()
		description += "JType=MOTOR:"
		description += "Npts=%i:" % self.getNumberOfPoints()
		description += "Counts=%d:" % self.getCounts()
		description += "Prefac=%d:" % self.getPrefactor()
		description += "DetectorType=" + self.getDetectorType() + ":"
		description += "CountType=" + self.getCountType() + ":"
		description += "Range=" + self.getMotor() + "="
		description += "%f %f" % (self.getRangeStart(),  self.getRangeStop())
		description += " s"
		
		return description


	def runScan(self,deleteScanAfterWards):
		scanDescription = ""
		scanName = self.getName();

		scanDescription =  self.buildScanDescription()

		scanDesrToListCommand = ScanDescrToListCommand(scanName,scanDescription)
		scanDesrToListCommand.run()

		scanRunCommand = ScanRunCommand(scanName)
		scanRunCommand.run()

		if(deleteScanAfterWards == True):
			self.deleteScanFromScanList()
		
