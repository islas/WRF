# Find pio
# Eventually replace with PIO's actual config if using that
# Once found this file will define:
#  PIO_FOUND - System has pnetcdf
#  PIO_INCLUDE_DIRS - The pnetcdf include directories
#  PIO_LIBRARIES - The libraries needed to use pnetcdf
#  PIO_DEFINITIONS - Compiler switches required for using pnetcdf

# exit early if we don't even need to be here
if ( PIO_FOUND )
  return()
endif()

# Use libpio.settings
find_file( 
          PIO_SETTINGS
          NAMES libpio.settings
          QUIET
          PATH_SUFFIXES lib lib64
          HINTS ENV PIO_ROOT ENV PIO ENV PIO_PATH
          )

if ( ${PIO_SETTINGS} MATCHES "-NOTFOUND$" )
  message( STATUS "No pio settings found" )
else()
  message( STATUS "Found PIO_SETTINGS : ${PIO_SETTINGS}" )
  file( READ ${PIO_SETTINGS} PIO_SETTINGS_LINES )
  set(
      PIO_SETTING_QUERIES
      "PIO Version"
      "Install Prefix"
      "PnetCDF Support"
      "NetCDF Integration"
      "Fortran"
      )
  set(
      PIO_SETTING_NAMES
      PIO_VERSION
      PIO_PREFIX
      PIO_HAS_PNETCDF
      PIO_HAS_NETCDF4
      PIO_HAS_FORTRAN
      )
  foreach( PIO_QUERY PIO_SETTING_NAME IN ZIP_LISTS PIO_SETTING_QUERIES PIO_SETTING_NAMES )
    string( REGEX MATCH "${PIO_QUERY}:[ \t]*([^ ]+)\n|$" SETTING_LINE ${PIO_SETTINGS_LINES} )

    if ( "${SETTING_LINE}" STREQUAL "" )
      # No setting found
      set( ${PIO_SETTING_NAME} NOTFOUND )
    else()
      string( REGEX REPLACE " |\t" ";" SETTING_LINE_LIST ${SETTING_LINE} )
      list( GET SETTING_LINE_LIST -1 ${PIO_SETTING_NAME} )
      string( REPLACE "\n" "" ${PIO_SETTING_NAME} ${${PIO_SETTING_NAME}} )
      if ( "${${PIO_SETTING_NAME}}" STREQUAL "yes" )
        set( ${PIO_SETTING_NAME} YES )
      elseif( "${${PIO_SETTING_NAME}}" STREQUAL "no" )
        set( ${PIO_SETTING_NAME} NO )
      endif()
    endif()
  endforeach()

  find_path(
            PIO_INCLUDE_DIR
            NAMES pio.h
            PATH_SUFFIXES include
            PATHS ${PIO_PREFIX}
            NO_DEFAULT_PATH
            )

  # Because we may need this for in-situ manual preprocessing do not use genex
  set( PIO_INCLUDE_DIRS ${PIO_INCLUDE_DIR} )

  # Find the actual name of the library
  find_library(
                PIO_C_LIBRARY
                pioc
                PATHS ${PIO_PREFIX}
                PATH_SUFFIXES lib lib64
                NO_DEFAULT_PATH
                )
  find_library(
                PIO_F_LIBRARY
                piof
                PATHS ${PIO_PREFIX}
                PATH_SUFFIXES lib lib64
                NO_DEFAULT_PATH
                )

  if ( NOT ${PIO_C_LIBRARY} MATCHES "-NOTFOUND$" )
    set( PIO_LIBRARIES ${PIO_LIBRARIES} ${PIO_C_LIBRARY} )
    set( PIO_LIBRARY ${PIO_C_LIBRARY} )
  endif()

  if ( NOT ${PIO_F_LIBRARY} MATCHES "-NOTFOUND$" )
    set( PIO_LIBRARIES ${PIO_LIBRARIES} ${PIO_F_LIBRARY} )
  endif()

  set( PIO_DEFINITIONS  )
endif()

include(FindPackageHandleStandardArgs)

# handle the QUIETLY and REQUIRED arguments and set PIO_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(
                                  PIO
                                  FOUND_VAR PIO_FOUND
                                  REQUIRED_VARS
                                    PIO_PREFIX
                                    PIO_INCLUDE_DIRS
                                    PIO_LIBRARIES
                                    PIO_VERSION
                                  VERSION_VAR PIO_VERSION
                                  HANDLE_VERSION_RANGE
                                )

if ( PIO_FOUND )
  # Do pnetcdf first as that might affect finding netcdf4
  if ( ${PIO_HAS_PNETCDF} )
    # PIO doesn't give quite a lot of info on where this netCDF might be
    find_package( PnetCDF REQUIRED )
  endif()

  if ( ${PIO_HAS_NETCDF4} )
    # PIO doesn't give quite a lot of info on where this netCDF might be
    find_package( netCDF REQUIRED )
  endif()

  if ( NOT TARGET PIO::PIO )
    set( PIO_LINK_INTERFACE_LANGUAGES )
    list( APPEND PIO_LINK_INTERFACE_LANGUAGES C )
    add_library( PIO::PIO_C UNKNOWN IMPORTED )
    set_target_properties(
                          PIO::PIO_C
                          PROPERTIES
                            IMPORTED_LOCATION                   "${PIO_C_LIBRARY}"
                            IMPORTED_LINK_INTERFACE_LANGUAGES   "${PIO_LINK_INTERFACE_LANGUAGES}"
                            INTERFACE_INCLUDE_DIRECTORIES       "${PIO_INCLUDE_DIRS}"
                          )
    if ( ${PIO_HAS_NETCDF4} )
      target_link_libraries( PIO::PIO_C INTERFACE netCDF::netcdf )
    endif()
    if ( ${PIO_HAS_PNETCDF} )
      target_link_libraries( PIO::PIO_C INTERFACE PnetCDF::pnetcdf )
    endif()
  endif()

  if( ${PIO_HAS_FORTRAN} AND NOT TARGET PIO::PIO_Fortran )
    list( APPEND PIO_LINK_INTERFACE_LANGUAGES Fortran )
    add_library( PIO::PIO_Fortran UNKNOWN IMPORTED )
    set_target_properties(
                          PIO::PIO_Fortran
                          PROPERTIES
                            IMPORTED_LOCATION                   "${PIO_F_LIBRARY}"
                            IMPORTED_LINK_INTERFACE_LANGUAGES   "${PIO_LINK_INTERFACE_LANGUAGES}"
                            INTERFACE_INCLUDE_DIRECTORIES       "${PIO_INCLUDE_DIRS}"
                          )
    target_link_libraries( PIO::PIO_Fortran INTERFACE PIO::PIO_C )
  endif()
endif()

mark_as_advanced( PIO_PREFIX PIO_C_LIBRARY PIO_F_LIBRARY )
