# Find netcdf
# Eventually replace with netCDF's actual config if using that
# Once found this file will define:
#  netCDF_FOUND - System has netcdf
#  netCDF_INCLUDE_DIRS - The netcdf include directories
#  netCDF_LIBRARIES - The libraries needed to use netcdf
#  netCDF_DEFINITIONS - Compiler switches required for using netcdf

# list( REMOVE_ITEM CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR} )
# find_package( netCDF )
# list( APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR} )

# exit early if we don't even need to be here
if ( netCDF_FOUND )
  return()
endif()


###############################################################################
# First try to find using netCDF native cmake build
# TODO : Enable this when netCDF native cmake build works well as an imported package
# find_package( netCDF CONFIG )
# if ( netCDF_FOUND )
#   message( STATUS "Found netCDF through native cmake build" )
#   return()
# endif()
###############################################################################


# else
# Use nc-config
find_program( 
                NETCDF_PROGRAM
                nc-config
                QUIET
                )

if ( ${NETCDF_PROGRAM} MATCHES "-NOTFOUND$" )
  message( STATUS "No nc-config found" )
else()
  message( STATUS "Found NETCDF_PROGRAM : ${NETCDF_PROGRAM}" )

  execute_process( COMMAND ${NETCDF_PROGRAM} --includedir   OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE netCDF_INCLUDE_DIR )
  execute_process( COMMAND ${NETCDF_PROGRAM} --libdir       OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE netCDF_LIBRARY_DIR )
  execute_process( COMMAND ${NETCDF_PROGRAM} --prefix       OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE netCDF_PREFIX )
  execute_process( COMMAND ${NETCDF_PROGRAM} --libs         OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE netCDF_CLIBS   )
  execute_process( COMMAND ${NETCDF_PROGRAM} --version      OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE netCDF_VERSION_RAW )

  # These do not pull all options available from nc-config out, but rather mirrors what is available from netCDFConfig.cmake.in
  set(
      netCDF_QUERY_YES_OPTIONS
      dap
      dap2
      dap4
      nc2
      nc4
      hdf5
      hdf4
      pnetcdf
      parallel

      # These are not part of the config but used in this to provide the properly linking
      szlib
      zstd
      )

  foreach( NC_QUERY ${netCDF_QUERY_YES_OPTIONS} )
    execute_process( COMMAND ${NETCDF_PROGRAM} --has-${NC_QUERY} OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE netCDF_${NC_QUERY}_LOWERCASE )
    string( TOUPPER ${NC_QUERY}                     NC_QUERY_UPPERCASE )
    string( TOUPPER ${netCDF_${NC_QUERY}_LOWERCASE} NC_ANSWER_UPPERCASE )
    # Convert to netCDF_HAS_* = YES/NO - Note this cannot be generator expression if you want to use it during configuration time
    set( netCDF_HAS_${NC_QUERY_UPPERCASE} ${NC_ANSWER_UPPERCASE} )
  endforeach()

  # check for large file support
  find_file( netCDF_INCLUDE_FILE netcdf.h ${netCDF_INCLUDE_DIR} )
  file( READ ${netCDF_INCLUDE_FILE} netCDF_INCLUDE_FILE_STR )
  string( FIND "${netCDF_INCLUDE_FILE_STR}" "NC_FORMAT_64BIT_DATA" netCDF_LARGE_FILE_SUPPORT_FOUND )
  if ( ${netCDF_LARGE_FILE_SUPPORT_FOUND} EQUAL -1 )
    set( netCDF_LARGE_FILE_SUPPORT "NO" )
  else()
    set( netCDF_LARGE_FILE_SUPPORT "YES" )
  endif()

  # Sanitize version
  string( REPLACE " " ";" netCDF_VERSION_LIST ${netCDF_VERSION_RAW} )
  list( GET netCDF_VERSION_LIST -1 netCDF_VERSION )

  set( netCDF_DEFINITIONS  )
  set( netCDF_LIBRARIES
      # All supported language variants will need this regardless - this may conflict with the RPATH in any
      # supplemental packages so be careful to use compatible langauge versions of netCDF
      $<$<OR:$<LINK_LANGUAGE:C>,$<LINK_LANGUAGE:Fortran>>:${netCDF_CLIBS}>
      )
  # Because we may need this for in-situ manual preprocessing do not use genex
  set( netCDF_INCLUDE_DIRS ${netCDF_INCLUDE_DIR} )
endif()

find_package( PkgConfig )

include(FindPackageHandleStandardArgs)

# handle the QUIETLY and REQUIRED arguments and set netCDF_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args( netCDF  DEFAULT_MSG
                                   netCDF_INCLUDE_DIRS
                                   netCDF_LIBRARY_DIR
                                   netCDF_CLIBS
                                   netCDF_VERSION
                                  )

mark_as_advanced( netCDF_CLIBS netCDF_PREFIX netCDF_LIBRARY_DIR )