#----------------------------------------------------------------
# Generated CMake target import file for configuration "MinSizeRel".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "embree" for configuration "MinSizeRel"
set_property(TARGET embree APPEND PROPERTY IMPORTED_CONFIGURATIONS MINSIZEREL)
set_target_properties(embree PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_MINSIZEREL "CXX"
  IMPORTED_LINK_INTERFACE_LIBRARIES_MINSIZEREL "sys;math;simd;lexers;tasking"
  IMPORTED_LOCATION_MINSIZEREL "${_IMPORT_PREFIX}/lib/libembree3.a"
  )

list(APPEND _IMPORT_CHECK_TARGETS embree )
list(APPEND _IMPORT_CHECK_FILES_FOR_embree "${_IMPORT_PREFIX}/lib/libembree3.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
