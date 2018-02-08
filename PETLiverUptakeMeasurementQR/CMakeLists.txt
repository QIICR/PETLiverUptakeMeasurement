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
  )

#-----------------------------------------------------------------------------
if(BUILD_TESTING)

  # Register the unittest subclass in the main script as a ctest.
  # Note that the test will also be available at runtime.
  slicer_add_python_unittest(SCRIPT ${MODULE_NAME}.py
  SLICER_ARGS --additional-module-paths
    ${CMAKE_BINARY_DIR}/${Slicer_QTSCRIPTEDMODULES_LIB_DIR}
    ${DEPENDENCIES_ADDITIONAL_MODULE_PATHS})

  add_subdirectory(Testing)
endif()

