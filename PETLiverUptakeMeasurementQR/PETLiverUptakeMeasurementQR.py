import os
import unittest
import vtk, qt, ctk, slicer
from slicer.ScriptedLoadableModule import *
import logging
import getpass
from datetime import datetime
import json
import shutil

#
# PETLiverUptakeMeasurementQR
#

class PETLiverUptakeMeasurementQR(ScriptedLoadableModule):
  """Uses ScriptedLoadableModule base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def __init__(self, parent):
    ScriptedLoadableModule.__init__(self, parent)
    self.parent.title = "PET Liver Uptake Measurement"
    self.parent.categories = ["Quantification"]
    self.parent.dependencies = ["DCMQI"]
    self.parent.contributors = ["Christian Bauer (Univeristy of Iowa)"]
    self.parent.helpText = """
    Measurement of uptake in a liver reference region in an SUVbw normalized FDG-18 whole-body PET scan with export to DICOM. \
    <a href="https://www.slicer.org/slicerWiki/index.php/Documentation/Nightly/Modules/PETLiverUptakeMeasurement">Documentation.</a>
    """
    self.parent.acknowledgementText = """
    This work was partially funded by NIH grants U01-CA140206 and U24-CA180918.
""" # replace with organization, grant and thanks.

#
# PETLiverUptakeMeasurementQRWidget
#

class PETLiverUptakeMeasurementQRWidget(ScriptedLoadableModuleWidget):
  """Uses ScriptedLoadableModuleWidget base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """
  
  def __init__(self,parent=None):
    ScriptedLoadableModuleWidget.__init__(self, parent)
    self.slicerTempDir = slicer.util.tempDirectory()
    
  def setup(self):
    self.measurementsLogic = PETLiverUptakeMeasurementQRLogic()

    ScriptedLoadableModuleWidget.setup(self)
    
    # Instantiate and connect widgets ...

    #
    # Parameters Area
    #
    parametersCollapsibleButton = ctk.ctkCollapsibleButton()
    parametersCollapsibleButton.text = "Liver Uptake Region Segmentation"
    self.layout.addWidget(parametersCollapsibleButton)
    parametersFormLayout = qt.QFormLayout(parametersCollapsibleButton)
    #
    # input volume selector
    #
    self.inputSelector = slicer.qMRMLNodeComboBox()
    self.inputSelector.nodeTypes = ["vtkMRMLScalarVolumeNode"]
    self.inputSelector.addAttribute("vtkMRMLScalarVolumeNode", "DICOM.instanceUIDs", None)
    self.inputSelector.addAttribute("vtkMRMLScalarVolumeNode", "DICOM.MeasurementUnitsCodeValue", "{SUVbw}g/ml") # ensures that the input is a SUV normalized PET scan
    self.inputSelector.selectNodeUponCreation = True
    self.inputSelector.addEnabled = False
    self.inputSelector.removeEnabled = False
    self.inputSelector.noneEnabled = False
    self.inputSelector.showHidden = False
    self.inputSelector.showChildNodeTypes = False
    self.inputSelector.setMRMLScene( slicer.mrmlScene )
    self.inputSelector.setToolTip( "Input SUVbw normalized DICOM PET volume")
    self.inputSelector.connect("currentNodeChanged(bool)",self.refreshUIElements)
    parametersFormLayout.addRow("Input Volume", self.inputSelector)    
    
    self.segmentationSelector = slicer.qMRMLNodeComboBox()
    self.segmentationSelector.nodeTypes = ["vtkMRMLLabelMapVolumeNode"]
    self.segmentationSelector.selectNodeUponCreation = True
    self.segmentationSelector.addEnabled = True
    self.segmentationSelector.removeEnabled = True
    self.segmentationSelector.noneEnabled = False
    self.segmentationSelector.showHidden = False
    self.segmentationSelector.showChildNodeTypes = False
    self.segmentationSelector.setMRMLScene( slicer.mrmlScene )
    self.segmentationSelector.setToolTip( "Output liver reference region volume")
    self.segmentationSelector.connect("currentNodeChanged(bool)",self.refreshUIElements)
    parametersFormLayout.addRow("Output Volume", self.segmentationSelector)
    
    self.regionSelector = slicer.qMRMLNodeComboBox()
    self.regionSelector.nodeTypes = ["vtkMRMLAnnotationROINode"]
    self.regionSelector.selectNodeUponCreation = True
    self.regionSelector.addEnabled = True
    self.regionSelector.removeEnabled = True
    self.regionSelector.noneEnabled = False
    self.regionSelector.showHidden = False
    self.regionSelector.showChildNodeTypes = False
    self.regionSelector.setMRMLScene( slicer.mrmlScene )
    self.regionSelector.setToolTip( "Search Region")
    self.regionSelector.connect("currentNodeChanged(bool)",self.refreshUIElements)
    parametersFormLayout.addRow("Region", self.regionSelector)
    
    self.thresholdRangeSlider = ctk.ctkRangeWidget()
    self.thresholdRangeSlider.minimum=0.0
    self.thresholdRangeSlider.maximum=20.0
    self.thresholdRangeSlider.singleStep = 0.1
    self.thresholdRangeSlider.minimumValue = 1.0
    self.thresholdRangeSlider.maximumValue = 10.0
    self.thresholdRangeSlider.connect("minimumValueChanged(double)",self.refreshUIElements)
    self.thresholdRangeSlider.connect("maximumValueChanged(double)",self.refreshUIElements)
    parametersFormLayout.addRow("Thresholds", self.thresholdRangeSlider)    
    
    self.erosionSlider = ctk.ctkSliderWidget()
    self.erosionSlider.minimum=0.0
    self.erosionSlider.maximum=20.0
    self.erosionSlider.singleStep=0.1
    self.erosionSlider.value = 3.0
    self.erosionSlider.connect("valueChanged(double)",self.refreshUIElements)
    parametersFormLayout.addRow("Erosion", self.erosionSlider)
    
    self.segmentButton = qt.QPushButton("Segment Reference Region")
    self.segmentButton.toolTip = "Segment reference region and measure uptake"
    self.segmentButton.connect('clicked(bool)', self.onSegmentButton)
    parametersFormLayout.addRow("Segment",self.segmentButton)
    
    # Measurement results
    measurementsCollapsibleButton = ctk.ctkCollapsibleButton()
    measurementsCollapsibleButton.text = "Measurements"
    self.layout.addWidget(measurementsCollapsibleButton)
    parametersFormLayout = qt.QFormLayout(measurementsCollapsibleButton)
    
    self.meanValueLineEdit = qt.QLineEdit("-")
    self.meanValueLineEdit.setReadOnly(True)
    self.meanValueLineEdit.setToolTip( "Mean value in reference region")
    parametersFormLayout.addRow("Mean",self.meanValueLineEdit)
    self.stdValueLineEdit = qt.QLineEdit("-")
    self.stdValueLineEdit.setReadOnly(True)
    self.stdValueLineEdit.setToolTip( "Standard deviation in reference region")
    parametersFormLayout.addRow("Standard Deviation",self.stdValueLineEdit)
    self.medianValueLineEdit = qt.QLineEdit("-")
    self.medianValueLineEdit.setReadOnly(True)
    self.medianValueLineEdit.setToolTip( "Mean value in reference region")
    parametersFormLayout.addRow("Median",self.medianValueLineEdit)

    # Reporting Section
    dicomReportingCollapsibleButton = ctk.ctkCollapsibleButton()
    dicomReportingCollapsibleButton.text = "DICOM Reporting"
    self.layout.addWidget(dicomReportingCollapsibleButton)
    reportingFormLayout = qt.QFormLayout(dicomReportingCollapsibleButton)
    
    self.readerValueLineEdit = qt.QLineEdit(getpass.getuser())
    #self.readerValueLineEdit.setReadOnly(True)
    self.readerValueLineEdit.setToolTip( "Name of reader reporting the measurement results")
    reportingFormLayout.addRow("Reader Name",self.readerValueLineEdit)
        
    self.reportBoxesFrame = qt.QFrame()
    self.reportBoxesFrame.setLayout(qt.QHBoxLayout())
    self.reportBoxesFrame.layout().setSpacing(0)
    self.reportBoxesFrame.layout().setMargin(0)
    self.saveReportButton = qt.QPushButton("Save Report")
    self.saveReportButton.toolTip = "Create partially completed DICOM Structured Report which could be continued at a later time (work in progress)"
    self.saveReportButton.connect('clicked(bool)', self.onSaveReportButtonClicked)
    # note: We don't have the functionality to reload and edit a temporary report, so this option doesn't make sense for us
    #self.reportBoxesFrame.layout().addWidget(self.saveReportButton)
    self.completeReportButton = qt.QPushButton("Complete Report")
    self.completeReportButton.toolTip = "Create the completed DICOM Structured Report and save to database"
    self.completeReportButton.connect('clicked(bool)', self.onCompleteReportButtonClicked)
    self.reportBoxesFrame.layout().addWidget(self.completeReportButton)
    reportingFormLayout.addRow("Report",self.reportBoxesFrame)

    # Add vertical spacer
    self.layout.addStretch(1)
    
    self.refreshUIElements()
    
  def enter(self):
    pass

  def exit(self):
    pass
  
  def refreshUIElements(self):
    self.saveReportButton.enabled = False
    self.completeReportButton.enabled = False
    self.meanValueLineEdit.text = "-"
    self.stdValueLineEdit.text = "-"
    self.medianValueLineEdit.text = "-"
    self.segmentButton.enabled = self.inputSelector.currentNode()!=None and self.segmentationSelector.currentNode()!=None
  
  def onSegmentButton(self):
    inputVolume = self.inputSelector.currentNode()
    outputVolume = self.segmentationSelector.currentNode()
    region = self.regionSelector.currentNode()
    
    cliParams = {'inputVolume': inputVolume.GetID(), \
                 'outputVolume': outputVolume.GetID(), \
                 'lowerThreshold': self.thresholdRangeSlider.minimumValue, \
                 'upperThreshold': self.thresholdRangeSlider.maximumValue, \
                 'erosion': self.erosionSlider.value \
                 }
    if region: cliParams['region'] = region.GetID()

    pd = qt.QProgressDialog('Running PET Liver Uptake Measurement Algorithm...', 'Cancel', 0, 100, slicer.util.mainWindow())
    pd.setModal(True)
    pd.setMinimumDuration(0)
    pd.show()
    pd.setValue(30)
    cliNode = None
    cliNode = slicer.cli.run(slicer.modules.petliveruptakemeasurement, cliNode, cliParams, wait_for_completion=False)
    while cliNode.IsBusy():
      slicer.app.processEvents()
      if pd.wasCanceled:
        cliNode.Cancel()
    pd.setValue(100)
    if pd.wasCanceled:
      return

    if cliNode.GetStatusString() != 'Completed':
      qt.QMessageBox().warning(None,"Warning","Segmentation of a liver reference region for uptake measurement was not successful. Please adjust parameters and try again.")
      return
    
    for i in range(cliNode.GetNumberOfParametersInGroup(1)):
      name = cliNode.GetParameterName(1,i)
      value = cliNode.GetParameterDefault(1,i)
      if name=='Mean_s': self.meanValueLineEdit.setText(value)
      if name=='Std_s': self.stdValueLineEdit.setText(value)
      if name=='Median_s': self.medianValueLineEdit.setText(value)    
    
    if self.meanValueLineEdit.text=='--':
      qt.QMessageBox().warning(None,"Warning","Segmentation of a liver reference region for uptake measurement was not successful. Please adjust parameters and try again.")
      return
    
    self.saveReportButton.enabled = True
    self.completeReportButton.enabled = True    
  
  # from https://github.com/QIICR/QuantitativeReporting/blob/master/Py/QuantitativeReporting.py
  def onSaveReportButtonClicked(self):
    success = self.saveReport()
    self.saveReportButton.enabled = not success
    if success:
      slicer.util.infoDisplay("Report successfully saved into SlicerDICOMDatabase")
  
  # from https://github.com/QIICR/QuantitativeReporting/blob/master/Py/QuantitativeReporting.py
  def onCompleteReportButtonClicked(self):
    success = self.saveReport(completed=True)
    self.saveReportButton.enabled = not success
    self.completeReportButton.enabled = not success
    if success:
      slicer.util.infoDisplay("Report successfully completed and saved into SlicerDICOMDatabase")
  
  # from https://github.com/QIICR/QuantitativeReporting/blob/master/Py/QuantitativeReporting.py
  def saveReport(self, completed=False):
    try:
      dcmSegmentationPath = self.createSEG()
      self.createDICOMSR(dcmSegmentationPath, completed)
      self.addProducedDataToDICOMDatabase()
    except (RuntimeError, ValueError, AttributeError) as exc:
      slicer.util.warningDisplay(exc.message)
      return False
    finally:
      self.cleanupTemporaryData()
    return True
    
  # base on https://github.com/QIICR/QuantitativeReporting/blob/master/Py/QuantitativeReporting.py
  def createSEG(self):   
    data = dict()
    data.update(self._getSeriesAttributes())
    description = "Automatic Liver Reference Region Segmentation"
    if not self.isAlgorithmAutomatic():
      description = "Semiautomatic Liver Reference Region Segmentation"
    data["SeriesDescription"] = description
    data.update(self._getAdditionalSeriesAttributes())
    data["segmentAttributes"] = self.generateJSON4DcmSEGExport()

    logging.debug("DICOM SEG Metadata output:")
    logging.debug(data)

    self.currentDateTime = datetime.now().strftime('%Y-%m-%d_%H%M%S')
    self.tempDir = os.path.join(self.slicerTempDir, self.currentDateTime)
    os.mkdir(self.tempDir)
    
    self.tempDicomDir = os.path.join(self.slicerTempDir, self.currentDateTime+'_dicoms')
    os.mkdir(self.tempDicomDir)
    for dicomFile in self.getDICOMFileList(self.inputSelector.currentNode(), absolutePaths=True):
      shutil.copy(dicomFile, self.tempDicomDir)    
    
    segmentFiles = []
    slicer.util.saveNode(self.segmentationSelector.currentNode(), os.path.join(self.tempDir, "seg.nrrd")), 
    segmentFiles.append(os.path.join(self.tempDir, "seg.nrrd"))       

    metaFilePath = self.saveJSON(data, os.path.join(self.tempDir, "seg_meta.json"))
    outputSegmentationPath = os.path.join(self.tempDir, "seg.dcm")

    params = {"dicomDirectory": self.tempDicomDir,
              "segImageFiles": ', '.join(segmentFiles).replace(', ', ","),
              "metaDataFileName": metaFilePath,
              "outputSEGFileName": outputSegmentationPath}

    logging.debug(params)

    cliNode = None
    cliNode = slicer.cli.run(slicer.modules.itkimage2segimage, cliNode, params, wait_for_completion=True)
    waitCount = 0
    while cliNode.IsBusy() and waitCount < 20:
      slicer.util.delayDisplay("Running SEG Encoding... %d" % waitCount, 1000)
      waitCount += 1
   
    if cliNode.GetStatusString() != 'Completed':
      raise RuntimeError("itkimage2segimage CLI did not complete cleanly")

    if not os.path.exists(outputSegmentationPath):
      raise RuntimeError("DICOM Segmentation was not created. Check Error Log for further information.")

    logging.debug("Saved DICOM Segmentation to {}".format(outputSegmentationPath))
    return outputSegmentationPath

  # from https://github.com/QIICR/QuantitativeReporting/blob/master/Py/QuantitativeReporting.py
  def createDICOMSR(self, referencedSegmentation, completed):
    data = self._getSeriesAttributes()
    data["SeriesDescription"] = "Liver Reference Region Measurement Report"

    compositeContextDataDir, data["compositeContext"] = os.path.dirname(referencedSegmentation), [os.path.basename(referencedSegmentation)]
    imageLibraryDataDir, data["imageLibrary"] = self.getDICOMFileList(self.inputSelector.currentNode())
    data.update(self._getAdditionalSRInformation(completed))

    data["Measurements"] = self.generateJSON4DcmSR(referencedSegmentation,
                                                   self.inputSelector.currentNode())
    
    logging.debug("DICOM SR Metadata output:")
    logging.debug(json.dumps(data, indent=2, separators=(',', ': ')))

    metaFilePath = self.saveJSON(data, os.path.join(self.tempDir, "sr_meta.json"))
    outputSRPath = os.path.join(self.tempDir, "sr.dcm")

    params = {"metaDataFileName": metaFilePath,
              "compositeContextDataDir": compositeContextDataDir,
              "imageLibraryDataDir": self.tempDicomDir,#imageLibraryDataDir,
              "outputFileName": outputSRPath}

    logging.debug(params)
    cliNode = None
    cliNode = slicer.cli.run(slicer.modules.tid1500writer, cliNode, params, wait_for_completion=True)
    waitCount = 0
    while cliNode.IsBusy() and waitCount < 20:
      slicer.util.delayDisplay("Running SR Encoding... %d" % waitCount, 1000)
      waitCount += 1

    if cliNode.GetStatusString() != 'Completed':
      raise Exception("tid1500writer CLI did not complete cleanly")

  # from https://github.com/QIICR/QuantitativeReporting/blob/master/Py/QuantitativeReporting.py
  def addProducedDataToDICOMDatabase(self):
    indexer = ctk.ctkDICOMIndexer()
    indexer.addFile(slicer.dicomDatabase, os.path.join(self.tempDir,'seg.dcm'), "copy")  # Note: doesn't really expect a destination dir
    indexer.addFile(slicer.dicomDatabase, os.path.join(self.tempDir,'sr.dcm'), "copy")  # Note: doesn't really expect a destination dir

  # from https://github.com/QIICR/QuantitativeReporting/blob/master/Py/QuantitativeReporting.py
  def cleanupTemporaryData(self):
    try:
      import shutil
      logging.debug("Cleaning up temporarily created directory {}".format(self.tempDir))
      logging.debug("Cleaning up temporarily created dicom image directory {}".format(self.tempDicomDir))
      shutil.rmtree(self.tempDir)
      shutil.rmtree(self.tempDicomDir)
      self.tempDir = None
      self.tempDicomDir = None
    except AttributeError:
      pass
  
  # from https://github.com/QIICR/QuantitativeReporting/blob/master/Py/QuantitativeReporting.py
  def _getSeriesAttributes(self):
    attributes = dict()
    sourceFileName = self.inputSelector.currentNode().GetStorageNode().GetFileName()
    DICOMTAGS_SERIES_NUMBER = '0020,0011'
    self.seriesNumber = self._getDICOMValue(sourceFileName, DICOMTAGS_SERIES_NUMBER)
    try: 
      self.seriesNumber = "100" if self.seriesNumber in [None,''] else str(int(self.seriesNumber)+100)
    except ValueError: # series "number" is not an integer
      self.seriesNumber = "100"
    attributes["SeriesNumber"] = self.seriesNumber
    attributes["InstanceNumber"] = "1"
    return attributes
  
  # from https://github.com/SlicerProstate/SlicerProstate/blob/master/SlicerProstate/SlicerProstateUtils/mixins.py
  def _getDICOMValue(self, currentFile, tag, default=""):
    try:
      return slicer.dicomDatabase.fileValue(str(currentFile), str(tag))
    except RuntimeError:
      logging.info("There are problems with accessing DICOM value %s from file %s" % (tag, currentFile))
    return default

  # from https://github.com/QIICR/QuantitativeReporting/blob/master/Py/QuantitativeReporting.py
  def _getAdditionalSeriesAttributes(self):
    attributes = { "ClinicalTrialSeriesID": "1",
                   "ClinicalTrialTimePointID": "1",
                   "ClinicalTrialCoordinatingCenterName": "QIICR UIowa"}
    if self.isAlgorithmAutomatic():
      attributes["ContentCreatorName"] = "Automated Liver Reference Region Segmentation"
    else:
      attributes["ContentCreatorName"] = str(self.readerValueLineEdit.text)
    return attributes
  
  def isAlgorithmAutomatic(self):
    return ( self.regionSelector.currentNode()==None and
             self.thresholdRangeSlider.minimumValue==1.0 and
             self.thresholdRangeSlider.maximumValue==10.0 and
             self.erosionSlider.value==3.0 )

  # from https://github.com/QIICR/QuantitativeReporting/blob/master/Py/QuantitativeReporting.py
  def _getAdditionalSRInformation(self, completed=False):
    data = dict()
    if self.isAlgorithmAutomatic():
      # deviceUID according to: https://github.com/QIICR/Iowa2DICOM/blob/master/Scripts/encode_all_objects.py#L50
      data["observerContext"] = {"ObserverType": "DEVICE", "DeviceObserverUID": "2.25.309174254785529415496603545181867685804"} 
    else:
      data["observerContext"] = {"ObserverType": "PERSON", "PersonObserverName": str(self.readerValueLineEdit.text)}
    data["VerificationFlag"] = "VERIFIED" if completed else "UNVERIFIED"
    data["CompletionFlag"] = "COMPLETE" if completed else "PARTIAL"
    data["activitySession"] = "1"
    data["timePoint"] = "1"
    return data

  # from https://github.com/QIICR/QuantitativeReporting/blob/master/Py/QuantitativeReporting.py
  def saveJSON(self, data, destination):
    with open(os.path.join(destination), 'w') as outfile:
      json.dump(data, outfile, indent=2)
    return destination

  # from https://github.com/QIICR/QuantitativeReporting/blob/master/Py/QuantitativeReporting.py
  def getDICOMFileList(self, volumeNode, absolutePaths=False):
    attributeName = "DICOM.instanceUIDs"
    instanceUIDs = volumeNode.GetAttribute(attributeName)
    if not instanceUIDs:
      raise ValueError("VolumeNode {0} has no attribute {1}".format(volumeNode.GetName(), attributeName))
    fileList = []
    rootDir = None
    for uid in instanceUIDs.split():
      rootDir, filename = self.getInstanceUIDDirectoryAndFileName(uid)
      fileList.append(str(filename if not absolutePaths else os.path.join(rootDir, filename)))
    if not absolutePaths:
      return rootDir, fileList
    return fileList

  # from https://github.com/QIICR/QuantitativeReporting/blob/master/Py/QuantitativeReporting.py
  def getInstanceUIDDirectoryAndFileName(self, uid):
    path = slicer.dicomDatabase.fileForInstance(uid)
    return os.path.dirname(path), os.path.basename(path)

  # based on https://github.com/QIICR/QuantitativeReporting/blob/af7001e2492ba168a50740a4470f502b172a2a34/Py/QuantitativeReporting.py#L933
  def generateJSON4DcmSEGExport(self):
    segmentsData = []
    segmentData = dict()
    segmentData["labelID"] = 1
    segmentData["SegmentDescription"] = "Liver Reference Region"
    segmentData["SegmentAlgorithmType"] = "AUTOMATIC" if self.isAlgorithmAutomatic() else "SEMIAUTOMATIC"
    segmentData["SegmentAlgorithmName"] = "Iowa QIN Liver Segmentation"
    segmentData["recommendedDisplayRGBValue"] = [255, 0, 0]
    segmentData["SegmentedPropertyCategoryCodeSequence"] = { \
      "CodeValue": "R-42018", \
      "CodingSchemeDesignator": "SRT", \
      "CodeMeaning": "Spatial and Relational Concept"
      }
    segmentData["SegmentedPropertyTypeCodeSequence"] = { \
        "CodeValue": "C94970",
        "CodingSchemeDesignator": "NCIt",
        "CodeMeaning": "Reference Region"
      }
    segmentData["AnatomicRegionSequence"] = { \
        "CodeValue": "T-6200",
        "CodingSchemeDesignator": "SRT",
        "CodeMeaning": "Liver"
      }
    segmentsData.append([segmentData]) # note: the extra brackets are necessary, because the segment attributes are a nested array
    return segmentsData
    
  # based on: https://github.com/QIICR/QuantitativeReporting/blob/master/Py/QuantitativeReporting.py#L999
  def generateJSON4DcmSR(self, dcmSegmentationFile, sourceVolumeNode):
    measurements = []

    sourceImageSeriesUID = self._getDICOMValue(sourceVolumeNode.GetStorageNode().GetFileName(), '0020,000E')
    logging.debug("SourceImageSeriesUID: {}".format(sourceImageSeriesUID))
    segmentationSOPInstanceUID = self._getDICOMValue(dcmSegmentationFile, "0008,0018")
    logging.debug("SegmentationSOPInstanceUID: {}".format(segmentationSOPInstanceUID))
    
    data = dict()
    data["TrackingIdentifier"] = "Liver reference region"
    data["ReferencedSegment"] = 1
    data["SourceSeriesForImageSegmentation"] = sourceImageSeriesUID
    data["segmentationSOPInstanceUID"] = segmentationSOPInstanceUID
    rwvid = sourceVolumeNode.GetAttribute("DICOM.RWV.instanceUID")
    if rwvid:
      data["rwvmMapUsedForMeasurement"] = rwvid
    data["Finding"] = { \
      "CodeValue": "C94970",
      "CodingSchemeDesignator": "NCIt",
      "CodeMeaning": "Reference Region"
      }
    data["FindingSite"] =  { \
      "CodeValue": "T-6200",
      "CodingSchemeDesignator": "SRT",
      "CodeMeaning": "Liver"
      }
    data["measurementItems"] = self.createMeasurementItems(sourceVolumeNode)
    measurements.append(data)

    return measurements
  
  # based on: https://github.com/QIICR/QuantitativeReporting/blob/master/Py/QuantitativeReporting.py#L1030
  def createMeasurementItems(self, sourceVolumeNode):
    measurementItems = []
    #sourceVolumeNode.GetAttribute("DICOM.MeasurementUnitsCodeValue") # should return "{SUVbw}g/ml" for our case
    #note: imaging modality and units have already been verified by input volume selector
    
    # mean
    item = dict()
    item["value"] = self.meanValueLineEdit.text
    item["quantity"] = {"CodeValue": "126401", "CodingSchemeDesignator": "DCM", "CodeMeaning": "SUVbw"}
    item["units"] = {"CodeValue": "{SUVbw}g/ml", "CodingSchemeDesignator": "UCUM", "CodeMeaning": "Standardized Uptake Value body weight"}
    item["derivationModifier"] = {"CodeValue": "R-00317", "CodingSchemeDesignator": "SRT", "CodeMeaning": "Mean" }
    measurementItems.append(item)
    
    # standard deviation
    item = dict()
    item["value"] = self.stdValueLineEdit.text
    item["quantity"] = {"CodeValue": "126401", "CodingSchemeDesignator": "DCM", "CodeMeaning": "SUVbw"}
    item["units"] = {"CodeValue": "{SUVbw}g/ml", "CodingSchemeDesignator": "UCUM", "CodeMeaning": "Standardized Uptake Value body weight"}
    item["derivationModifier"] = {"CodeValue": "R-10047", "CodingSchemeDesignator": "SRT", "CodeMeaning": "Standard Deviation" }
    measurementItems.append(item)
    
    # median
    item = dict()
    item["value"] = self.medianValueLineEdit.text
    item["quantity"] = {"CodeValue": "126401", "CodingSchemeDesignator": "DCM", "CodeMeaning": "SUVbw"}
    item["units"] = {"CodeValue": "{SUVbw}g/ml", "CodingSchemeDesignator": "UCUM", "CodeMeaning": "Standardized Uptake Value body weight"}
    item["derivationModifier"] = {"CodeValue": "R-00319", "CodingSchemeDesignator": "SRT", "CodeMeaning": "Median" }
    measurementItems.append(item)
    
    return measurementItems

#
# PETLiverUptakeMeasurementQRLogic
#

class PETLiverUptakeMeasurementQRLogic(ScriptedLoadableModuleLogic):

  def __init__(self):
    pass

  def __del__(self):
    pass

#
# PETLiverUptakeMeasurementQRTest
#

class PETLiverUptakeMeasurementQRTest(ScriptedLoadableModuleTest):
  """
  This is the test case for your scripted module.
  Uses ScriptedLoadableModuleTest base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def setUp(self):
    """ Do whatever is needed to reset the state - typically a scene clear will be enough.
    """
    slicer.mrmlScene.Clear(0)

  def runTest(self):
    """Run as few or as many tests as needed here.
    """
    self.setUp()
    self.test_PETLiverUptakeMeasurementQR1()

  def test_PETLiverUptakeMeasurementQR1(self):
    pass
