#----------------------------------------------------------------
# Generated CMake target import file for configuration "MinSizeRel".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "lexers" for configuration "MinSizeRel"
set_property(TARGET lexers APPEND PROPERTY IMPORTED_CONFIGURATIONS MINSIZEREL)
set_target_properties(lexers PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_MINSIZEREL "CXX"
  IMPORTED_LINK_INTERFACE_LIBRARIES_MINSIZEREL "sys;math"
  IMPORTED_LOCATION_MINSIZEREL "${_IMPORT_PREFIX}/lib/liblexers.a"
  )

list(APPEND _IMPORT_CHECK_TARGETS lexers )
list(APPEND _IMPORT_CHECK_FILES_FOR_lexers "${_IMPORT_PREFIX}/lib/liblexers.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
