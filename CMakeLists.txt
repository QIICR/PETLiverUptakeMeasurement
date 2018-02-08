cmake_minimum_required(VERSION 2.8.9)

project(PETLiverUptakeMeasurement)
set(MODULE_NAME PETLiverUptakeMeasurement)

#-----------------------------------------------------------------------------
# Extension meta-information
set(EXTENSION_HOMEPAGE "https://www.slicer.org/slicerWiki/index.php/Documentation/Nightly/Extensions/PETLiverUptakeMeasurement")
set(EXTENSION_CATEGORY "Quantification")
set(EXTENSION_CONTRIBUTORS "Christian Bauer (University of Iowa), Markus van Tol (University of Iowa), Andrey Fedorov (SPL), Ethan Ulrich (University of Iowa), Reinhard Beichel (University of Iowa), John Buatti (University of Iowa)")
set(EXTENSION_DESCRIPTION "This extensions contains a module for automated measurement of liver uptake in whole-body FDG-18 PET scans.")
set(EXTENSION_ICONURL "https://raw.githubusercontent.com/QIICR/PETLiverUptakeMeasurement/master/PETLiverUptakeMeasurement.png")
set(EXTENSION_SCREENSHOTURLS "https://www.slicer.org/slicerWiki/images/6/6d/PETLiverUptakeMeasurementWithModel.png")
set(EXTENSION_DEPENDS "DCMQI")

#-----------------------------------------------------------------------------
# Extension dependencies
find_package(Slicer REQUIRED)
include(${Slicer_USE_FILE})

#find_package(DCMQI REQUIRED)
set(DEPENDENCIES_ADDITIONAL_MODULE_PATHS
  ${DCMQI_DIR}/bin
  ${DCMQI_DIR}/${Slicer_CLIMODULES_LIB_DIR}
  )
  
#-----------------------------------------------------------------------------
# Extension modules
add_subdirectory(PETLiverUptakeMeasurement)
add_subdirectory(PETLiverUptakeMeasurementQR)
## NEXT_MODULE

#-----------------------------------------------------------------------------
include(${Slicer_EXTENSION_GENERATE_CONFIG})
include(${Slicer_EXTENSION_CPACK})
