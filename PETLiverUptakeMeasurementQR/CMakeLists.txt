cmake_minimum_required(VERSION 2.8.6)

#-----------------------------------------------------------------------------
set(MODULE_NAME PETLiverUptakeMeasurementQR)

#-----------------------------------------------------------------------------
set(MODULE_PYTHON_SCRIPTS
  PETLiverUptakeMeasurementQR.py
  )

set(MODULE_PYTHON_RESOURCES
  )

#-----------------------------------------------------------------------------
SlicerMacroBuildScriptedModule(
  NAME PETTumorSegmentationEffect
  SCRIPTS "${MODULE_PYTHON_SCRIPTS}"
  RESOURCES "${MODULE_PYTHON_RESOURCES}"
  WITH_GENERIC_TESTS
  )

#-----------------------------------------------------------------------------
if(BUILD_TESTING)
  add_subdirectory(Testing)
endif()

