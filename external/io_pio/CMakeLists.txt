# WRF CMake Build
# PIO requires using module_domain and thus cannot be compiled separately
target_sources(
                ${PROJECT_NAME}_Core
                PRIVATE
                  wrf_data_pio.F90
                  pio_routines.F90
                  wrf_io.F90
                  field_routines.F90
                  read_bdy_routines.F90
                  module_wrfsi_static.F90
                  )


target_link_libraries( 
                      ${PROJECT_NAME}_Core
                      PRIVATE
                        PIO::PIO_Fortran
                      )

# target_include_directories( ${FOLDER_COMPILE_TARGET}
#                             PUBLIC
#                               $<TARGET_PROPERTY:${FOLDER_COMPILE_TARGET},Fortran_MODULE_DIRECTORY>
#                               $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/external/ioapi_share>
#                               $<INSTALL_INTERFACE:include/external/ioapi_share>
#                             PRIVATE
#                               ${CMAKE_CURRENT_SOURCE_DIR}
#                             )

# install(
#         TARGETS ${FOLDER_COMPILE_TARGET}
#         EXPORT  ${EXPORT_NAME}Targets
#         RUNTIME DESTINATION bin/
#         ARCHIVE DESTINATION lib/
#         LIBRARY DESTINATION lib/
#         )
