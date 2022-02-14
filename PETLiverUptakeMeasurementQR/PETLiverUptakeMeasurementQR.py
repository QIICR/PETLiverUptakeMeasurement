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
    self.parent.dependencies = []
    self.parent.contributors = ["Christian Bauer (University of Iowa)"]
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
    #self.inputSelector.addAttribute("vtkMRMLScalarVolumeNode", "DICOM.MeasurementUnitsCodeValue", "{SUVbw}g/ml") # ensures that the input is a SUV normalized PET scan
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
    if self.inputSelector.currentNode().GetStorageNode():
      sourceFileName = self.inputSelector.currentNode().GetStorageNode().GetFileName()
    else:
      sourceFileName = slicer.dicomDatabase.fileForInstance(self.inputSelector.currentNode().GetAttribute("DICOM.instanceUIDs").split(" ")[0])
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
    if sourceVolumeNode.GetStorageNode():
      sourceFileName = sourceVolumeNode.GetStorageNode().GetFileName()
    else:
      sourceFileName = slicer.dicomDatabase.fileForInstance(sourceVolumeNode.GetAttribute("DICOM.instanceUIDs").split(" ")[0])
    sourceImageSeriesUID = self._getDICOMValue(sourceFileName, '0020,000E')
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

    quantity = {
      "CodeValue": sourceVolumeNode.GetVoxelValueQuantity().GetCodeValue(),
      "CodingSchemeDesignator": sourceVolumeNode.GetVoxelValueQuantity().GetCodingSchemeDesignator(),
      "CodeMeaning": sourceVolumeNode.GetVoxelValueQuantity().GetCodeMeaning()
      }

    units = {
      "CodeValue": sourceVolumeNode.GetVoxelValueUnits().GetCodeValue(),
      "CodingSchemeDesignator": sourceVolumeNode.GetVoxelValueUnits().GetCodingSchemeDesignator(),
      "CodeMeaning": sourceVolumeNode.GetVoxelValueUnits().GetCodeMeaning()
      }

    # mean
    item = dict()
    item["value"] = self.meanValueLineEdit.text
    item["quantity"] = quantity
    item["units"] = units
    item["derivationModifier"] = {"CodeValue": "R-00317", "CodingSchemeDesignator": "SRT", "CodeMeaning": "Mean" }
    measurementItems.append(item)

    # standard deviation
    item = dict()
    item["value"] = self.stdValueLineEdit.text
    item["quantity"] = quantity
    item["units"] = units
    item["derivationModifier"] = {"CodeValue": "R-10047", "CodingSchemeDesignator": "SRT", "CodeMeaning": "Standard Deviation" }
    measurementItems.append(item)

    # median
    item = dict()
    item["value"] = self.medianValueLineEdit.text
    item["quantity"] = quantity
    item["units"] =units
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

import pydicom
from DICOMLib import DICOMUtils
class PETLiverUptakeMeasurementQRTest(ScriptedLoadableModuleTest):
  """
  This is the test case for your scripted module.
  Uses ScriptedLoadableModuleTest base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def runTest(self):
    """Run as few or as many tests as needed here.
    """
    self.setUp()
    self.test_PETLiverUptakeMeasurementQR1()
    slicer.mrmlScene.Clear(0)
    self.test_PETLiverUptakeMeasurementQR2()
    self.tearDown()

  def setUp(self):
    """ Open temporary DICOM database
    """
    slicer.mrmlScene.Clear(0)
    self.delayMs = 700
    self.tempDataDir = os.path.join(slicer.app.temporaryPath,'PETTest')
    self.tempDicomDatabaseDir = os.path.join(slicer.app.temporaryPath,'PETTestDicom')

  def doCleanups(self):
    """ cleanup temporary data in case an exception occurs
    """
    self.tearDown()

  def tearDown(self):
    """ Close temporary DICOM database and remove temporary data
    """
    try:
      import shutil
      if os.path.exists(self.tempDataDir):
        shutil.rmtree(self.tempDataDir)
    except Exception as e:
      import traceback
      traceback.print_exc()
      self.delayDisplay('Test caused exception!\n' + str(e),self.delayMs*2)

  def loadTestData(self):
    self.patienName = 'QIN-HEADNECK-01-0139'
    #download data and add to dicom database
    zipFileUrl = 'http://github.com/QIICR/PETTumorSegmentation/releases/download/4.10.2/QIN-HEADNECK-01-0139-PET.zip'
    zipFilePath = self.tempDataDir+'/dicom.zip'
    zipFileData = self.tempDataDir+'/dicom'
    expectedNumOfFiles = 545
    if not os.access(self.tempDataDir, os.F_OK):
      os.mkdir(self.tempDataDir)
    if not os.access(zipFileData, os.F_OK):
      os.mkdir(zipFileData)
      slicer.util.downloadAndExtractArchive( zipFileUrl, zipFilePath, zipFileData, expectedNumOfFiles)
    DICOMUtils.importDicom(zipFileData)
    
    # load dataset
    dicomFiles = slicer.util.getFilesInDirectory(zipFileData)
    loadablesByPlugin, loadEnabled = DICOMUtils.getLoadablesFromFileLists([dicomFiles],['DICOMScalarVolumePlugin'])
    loadedNodeIDs = DICOMUtils.loadLoadables(loadablesByPlugin)
    imageNode = slicer.mrmlScene.GetNodeByID(loadedNodeIDs[0])
    imageNode.SetSpacing(3.3940266832237, 3.3940266832237, 2.02490234375) # mimic spacing as produced by Slicer 4.10 for which the test was originally developed
    imageNode.SetOrigin(285.367523193359375,494.58682250976556816,-1873.3819580078125) # mimic origin as produced by Slicer 4.10 for which the test was originally developed
    
    # apply the SUVbw conversion factor and set units and quantity
    suvNormalizationFactor = 0.00040166400000000007
    quantity = slicer.vtkCodedEntry()
    quantity.SetFromString('CodeValue:126400|CodingSchemeDesignator:DCM|CodeMeaning:Standardized Uptake Value')
    units = slicer.vtkCodedEntry()
    units.SetFromString('CodeValue:{SUVbw}g/ml|CodingSchemeDesignator:UCUM|CodeMeaning:Standardized Uptake Value body weight')
    multiplier = vtk.vtkImageMathematics()
    multiplier.SetOperationToMultiplyByK()
    multiplier.SetConstantK(suvNormalizationFactor)
    multiplier.SetInput1Data(imageNode.GetImageData())
    multiplier.Update()
    imageNode.GetImageData().DeepCopy(multiplier.GetOutput())
    imageNode.GetVolumeDisplayNode().SetWindowLevel(6,3)
    imageNode.GetVolumeDisplayNode().SetAndObserveColorNodeID('vtkMRMLColorTableNodeInvertedGrey')
    imageNode.SetVoxelValueQuantity(quantity)
    imageNode.SetVoxelValueUnits(units)

    return imageNode

  def test_PETLiverUptakeMeasurementQR1(self):
    """ test standard segmentation and report generation
    """
    try:
      self.assertIsNotNone( slicer.modules.petliveruptakemeasurement )
      with DICOMUtils.TemporaryDICOMDatabase(self.tempDicomDatabaseDir) as db:
        self.assertTrue(db.isOpen)
        self.assertEqual(slicer.dicomDatabase, db)

        self.delayDisplay('Loading PET DICOM dataset (including download if necessary)')
        petNode = self.loadTestData()

        self.delayDisplay('Running segmentation')
        m = slicer.util.mainWindow()
        m.moduleSelector().selectModule('PETLiverUptakeMeasurementQR')
        qrWidget = slicer.modules.PETLiverUptakeMeasurementQRWidget
        qrWidget.inputSelector.setCurrentNode(petNode)
        segmentationNode = qrWidget.segmentationSelector.addNode()
        qrWidget.segmentButton.click()

        self.assertTrue(abs(float(qrWidget.meanValueLineEdit.text)-2.36253)<0.01)
        self.assertTrue(abs(float(qrWidget.stdValueLineEdit.text)-0.402997)<0.01)
        self.assertTrue(abs(float(qrWidget.medianValueLineEdit.text)-2.335)<0.01)

        self.delayDisplay('Completing and writing DICOM report')
        qrWidget.readerValueLineEdit.text = 'autotest'
        self.assertTrue(qrWidget.saveReport(completed=True))

        self.delayDisplay('Checking for DICOM SEG and SR')
        patientUID = DICOMUtils.getDatabasePatientUIDByPatientName(self.patienName)
        studies = slicer.dicomDatabase.studiesForPatient(patientUID)
        series = slicer.dicomDatabase.seriesForStudy(studies[0])
        SRSeries = None
        SEGSeries = None
        for serie in series:
          description = slicer.dicomDatabase.descriptionForSeries(serie)
          if description=='Automatic Liver Reference Region Segmentation':
            SEGSeries = serie
          if description=='Liver Reference Region Measurement Report':
            SRSeries = serie
        self.assertIsNotNone(SRSeries)
        self.assertIsNotNone(SEGSeries)
        SRFile = slicer.dicomDatabase.filesForSeries(SRSeries)[0]

        self.delayDisplay('Loading DICOM SR and verifying stored measurements')
        sr = pydicom.read_file(SRFile)
        dicomMean = self._getMeasuredValue(sr,'Mean')
        self.assertIsNotNone(dicomMean)
        self.assertEqual(dicomMean.MeasuredValueSequence[0].NumericValue, 2.36253)
        dicomStandardDeviation = self._getMeasuredValue(sr,'Standard Deviation')
        self.assertIsNotNone(dicomStandardDeviation)
        self.assertEqual(dicomStandardDeviation.MeasuredValueSequence[0].NumericValue, 0.402997)
        dicomMedian = self._getMeasuredValue(sr,'Median')
        self.assertIsNotNone(dicomMedian)
        self.assertEqual(dicomMedian.MeasuredValueSequence[0].NumericValue, 2.335)

        # clean up data from DICOM database
        db.removePatient(patientUID)

        self.delayDisplay('Test passed!')

    except Exception as e:
      import traceback
      traceback.print_exc()
      self.delayDisplay('Test caused exception!\n' + str(e),self.delayMs*2)

  def _getMeasuredValue(self, item, ConceptCodeMeaning):
    """ traverse pyDicom object hierarchy to search for measurement sequence with matching ConceptCodeMeaning
    """
    if 'MeasuredValueSequence' in item and \
      'NumericValue' in item.MeasuredValueSequence[0] and \
      'ContentSequence' in item and \
      'ConceptCodeSequence' in item.ContentSequence[0] and \
      'CodeMeaning' in item.ContentSequence[0].ConceptCodeSequence[0] and \
      item.ContentSequence[0].ConceptCodeSequence[0].CodeMeaning==ConceptCodeMeaning:
      return item
    if 'ContentSequence' in item:
     for csItem in item.ContentSequence:
       r = self._getMeasuredValue(csItem, ConceptCodeMeaning)
       if r: return r
    return None

  def test_PETLiverUptakeMeasurementQR2(self):
    """ test segmentation options
    """
    try:
      self.assertIsNotNone( slicer.modules.petliveruptakemeasurement )
      with DICOMUtils.TemporaryDICOMDatabase(self.tempDicomDatabaseDir) as db:
        self.assertTrue(db.isOpen)
        self.assertEqual(slicer.dicomDatabase, db)
        self.delayDisplay('Loading PET DICOM dataset (including download if necessary)')
        petNode = self.loadTestData()

        qrWidget = slicer.modules.PETLiverUptakeMeasurementQRWidget
        qrWidget.inputSelector.setCurrentNode(petNode)
        segmentationNode = qrWidget.segmentationSelector.addNode()

        self.delayDisplay('Running segmentation with standard settings')
        qrWidget.segmentButton.click()
        self.assertTrue(abs(float(qrWidget.meanValueLineEdit.text)-2.36253)<0.01)

        self.delayDisplay('Specifying annotation ROI')
        roi = slicer.mrmlScene.AddNewNodeByClass('vtkMRMLAnnotationROINode')
        roi.SetXYZ([-34,243,-1168])
        roi.SetRadiusXYZ([85,102,82])
        roi.SetName('ROI')
        qrWidget.regionSelector.setCurrentNode(roi)
        qrWidget.segmentButton.click()
        self.assertTrue(abs(float(qrWidget.meanValueLineEdit.text)-2.91891)<0.01)

        self.delayDisplay('Changing erosion range')
        originalErosion = qrWidget.erosionSlider.value
        qrWidget.erosionSlider.value = 0
        qrWidget.segmentButton.click()
        self.assertTrue(abs(float(qrWidget.meanValueLineEdit.text)-2.71982)<0.01)

        self.delayDisplay('Changing thresholds')
        originalMinimValue = qrWidget.thresholdRangeSlider.minimumValue
        originalMaximumValue = qrWidget.thresholdRangeSlider.maximumValue
        qrWidget.thresholdRangeSlider.minimumValue = 2
        qrWidget.thresholdRangeSlider.maximumValue = 20
        qrWidget.segmentButton.click()
        self.assertTrue(abs(float(qrWidget.meanValueLineEdit.text)-3.72669)<0.01)

        self.delayDisplay('Completing and writing DICOM report')
        qrWidget.readerValueLineEdit.text = 'semiautotest'
        self.assertTrue(qrWidget.saveReport(completed=True))

        self.delayDisplay('Testing that report was saved as semiautomatic result')
        patientUID = DICOMUtils.getDatabasePatientUIDByPatientName('QIN-HEADNECK-01-0139')
        studies = slicer.dicomDatabase.studiesForPatient(patientUID)
        series = slicer.dicomDatabase.seriesForStudy(studies[0])
        SEGSeries = None
        for serie in series:
          description = slicer.dicomDatabase.descriptionForSeries(serie)
          if description=='Semiautomatic Liver Reference Region Segmentation':
            SEGSeries = serie
        self.assertIsNotNone(SEGSeries)

        # reset values
        qrWidget.regionSelector.removeCurrentNode()
        qrWidget.erosionSlider.value = originalErosion
        qrWidget.thresholdRangeSlider.minimumValue = originalMinimValue
        qrWidget.thresholdRangeSlider.maximumValue = originalMaximumValue

        # clean up data from DICOM database
        db.removePatient(patientUID)

        self.delayDisplay('Test passed!')

    except Exception as e:
      import traceback
      traceback.print_exc()
      self.delayDisplay('Test caused exception!\n' + str(e),self.delayMs*2)

  def _loadTestData_WithPETDICOMExtension(self):
    """ load SUV normalized PET scan from DICOM fileassuming Slicer-PETDICOM extension is installed and enabled
    """
    from DICOMLib import DICOMUtils
    import urllib
    data = {
      'PETVolume': {
        'UID': '1.3.6.1.4.1.14519.5.2.1.2744.7002.886851941687931416391879144903',
        'url': 'http://github.com/QIICR/PETTumorSegmentation/releases/download/4.10.2/QIN-HEADNECK-01-0139-PET.zip',
        'zipFile': 'QIN-HEADNECK-01-0139-PET.zip',
        'SUVNormalizationFactor': 0.00040166400000000007
      }
    }
    destinationDirectory = self.tempDataDir
    for key, value in data.iteritems(): # download data if necessary
      UID = value['UID']
      if not len(slicer.dicomDatabase.filesForSeries(UID)):
        url = value['url']
        zipFile = value['zipFile']
        filePath = os.path.join(destinationDirectory, zipFile)
        if not os.path.exists(os.path.dirname(filePath)):
          os.makedirs(os.path.dirname(filePath))
        logging.debug('Saving download %s to %s ' % (url, filePath))
        if not os.path.exists(filePath) or os.stat(filePath).st_size == 0:
          slicer.util.delayDisplay('Requesting download of %s...\n' % url, 1000)
          urllib.urlretrieve(url, filePath)
        if os.path.exists(filePath) and os.path.splitext(filePath)[1]=='.zip':
          success = slicer.app.applicationLogic().Unzip(filePath, destinationDirectory)
          if not success:
            logging.error("Archive %s was NOT unzipped successfully." %  filePath)
      indexer = ctk.ctkDICOMIndexer()
      indexer.addDirectory(slicer.dicomDatabase, destinationDirectory, None)
      indexer.waitForImportFinished()
    # load dataset
    success=DICOMUtils.loadSeriesByUID([data['PETVolume']['UID']])
    if not success:
      logging.error("Unable to load dicom data %s\n." %  data['PETVolume']['UID'])

    return slicer.mrmlScene.GetFirstNodeByClass('vtkMRMLScalarVolumeNode')

  def _test_PETLiverUptakeMeasurementQR2_WithQRExtension(self):
    """ verify SR measurements assuming QuantiativeReporting Extension is installed and enabled
    """
    from DICOMLib import DICOMUtils
    self.delayDisplay('Loading DICOM SEG and SR')

    # note: introduces dependency on QuantiativeReporting and DICOMPETExtension

    # with QuantiativeReporting installed, this should load the DICOM SR and SEG together with the PET scan
    # disable all dicom plugins except those necessary
    dicomWidget = slicer.modules.dicom.widgetRepresentation().self()
    dicomPluginCheckbox =  dicomWidget.detailsPopup.pluginSelector.checkBoxByPlugin
    dicomPluginStates = {(key,value.checked) for key,value in dicomPluginCheckbox.iteritems()}
    for cb in dicomPluginCheckbox.itervalues():
      cb.checked=False
    dicomPluginCheckbox['DICOMRWVMPlugin'].checked = True
    dicomPluginCheckbox['DICOMPETSUVPlugin'].checked = True
    dicomPluginCheckbox['DICOMSegmentationPlugin'].checked = True
    dicomPluginCheckbox['DICOMTID1500Plugin'].checked = True
    DICOMUtils.loadPatientByName(self.patienName)
    for key,value in dicomPluginStates:
      dicomPluginCheckbox[key].checked=value

    # with QuantiativeReporting installed, this should load the DICOM SR and SEG together with the PET scan
    DICOMUtils.loadPatientByName(self.patienName)
    tableNode = slicer.mrmlScene.GetFirstNodeByClass("vtkMRMLTableNode")
    self.assertIsNotNone(tableNode, "Loading SR into mrmlScene failed. No vtkMRMLTableNodes were found within the scene.")
    self.assertEqual(tableNode.GetNumberOfRows(),1)
    self.assertTrue(tableNode.GetNumberOfColumns()>3)
    self.assertEqual(tableNode.GetColumnName(0),'Segment')
    self.assertEqual(tableNode.GetColumnName(1),'Mean [{SUVbw}g/ml]')
    self.assertEqual(tableNode.GetColumnName(2),'Standard Deviation [{SUVbw}g/ml]')
    self.assertEqual(tableNode.GetColumnName(3),'Median [{SUVbw}g/ml]')
    self.assertEqual(tableNode.GetCellText(0,0),'Liver reference region')
    self.assertTrue(abs(float(tableNode.GetCellText(0,1))-2.36253)<0.01)
    self.assertTrue(abs(float(tableNode.GetCellText(0,2))-0.402997)<0.01)
    self.assertTrue(abs(float(tableNode.GetCellText(0,3))-2.335)<0.01)
    self.delayDisplay('Test passed!')

    for key,value in dicomPluginStates:
      dicomPluginCheckbox[key].checked=value
