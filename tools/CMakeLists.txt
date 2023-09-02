# WRF CMake Build

#!TODO ORGANIZE THIS FOLDER
set( FOLDER_COMPILE_TARGET registry )

add_executable(
               ${FOLDER_COMPILE_TARGET}
              )

set( GEN_COMMS gen_comms.stub )
if ( ${USE_MPI} )
  set( GEN_COMMS ${PROJECT_SOURCE_DIR}/external/RSL_LITE/gen_comms.c )
else()
  # Account for the weird makefile nonsense of copying things around
  set_source_files_properties( 
                              gen_comms.stub
                              TARGET_DIRECTORY ${FOLDER_COMPILE_TARGET}
                              PROPERTIES
                                LANGUAGE C
                              )
endif()

target_sources(
                ${FOLDER_COMPILE_TARGET}
                PRIVATE
                  registry.c
                  my_strtok.c
                  reg_parse.c
                  data.c
                  type.c
                  misc.c
                  gen_defs.c
                  gen_allocs.c
                  gen_mod_state_descr.c
                  gen_scalar_indices.c
                  gen_args.c
                  gen_config.c
                  sym.c
                  symtab_gen.c
                  gen_irr_diag.c
                  gen_model_data_ord.c
                  gen_interp.c
                  # gen_comms.c
                  ${GEN_COMMS}
                  gen_scalar_derefs.c
                  set_dim_strs.c
                  gen_wrf_io.c
                  gen_streams.c
                  )

# set_target_properties( 
#                       ${FOLDER_COMPILE_TARGET} 
#                         PROPERTIES
#                           Fortran_MODULE_DIRECTORY ${CMAKE_INSTALL_PREFIX}/${FOLDER_COMPILE_TARGET}
#                           Fortran_FORMAT           FREE                          
#                       )

target_include_directories( ${FOLDER_COMPILE_TARGET}
                            PRIVATE
                              ${CMAKE_CURRENT_SOURCE_DIR}
                            )

install(
        TARGETS ${FOLDER_COMPILE_TARGET}
        RUNTIME DESTINATION bin/
        ARCHIVE DESTINATION lib/
        LIBRARY DESTINATION lib/
        )

# Extra stuff for weird registry stuff
set( REGISTRY_FILE "NO_REGISTRY_FILE_SET" )
if ( ${WRF_CORE} STREQUAL "ARW" )
  
  if ( ${ENABLE_CHEM} )
    set( REGISTRY_FILE ${PROJECT_SOURCE_DIR}/Registry/Registry.EM_CHEM )

  # This check does nothing
  # elseif ( ${WRF_DFI_RADAR} )
  #   set( REGISTRY_FILE ${PROJECT_SOURCE_DIR}/Registry/Registry.EM )

  else()
    set( REGISTRY_FILE ${PROJECT_SOURCE_DIR}/Registry/Registry.EM )

  endif()

elseif ( ${WRF_CORE} STREQUAL "PLUS" )
    set( REGISTRY_FILE ${PROJECT_SOURCE_DIR}/Registry/Registry.tladj )

elseif ( ${WRF_CORE} STREQUAL "CONVERT" )
  set( REGISTRY_FILE ${PROJECT_SOURCE_DIR}/Registry/Registry.CONVERT )

elseif ( ${WRF_CORE} STREQUAL "DA" OR ${WRF_CORE} STREQUAL "DA_4D_VAR" )
  if ( ${WRF_CHEM} )
    set( REGISTRY_FILE ${PROJECT_SOURCE_DIR}/Registry/Registry.wrfchemvar )
  else()
    set( REGISTRY_FILE ${PROJECT_SOURCE_DIR}/Registry/Registry.wrfvar )
  endif()

endif()

get_directory_property( DIR_DEFS DIRECTORY ${CMAKE_SOURCE_DIR} COMPILE_DEFINITIONS )
wrf_expand_definitions( 
                        RESULT_VAR  REGISTRY_DEFS
                        DEFINITIONS ${DIR_DEFS}
                        )

# How this is not a bigger thing or not resolved is beyond me
# https://gitlab.kitware.com/cmake/cmake/-/issues/18005
# Also the suggestion does not work
add_custom_command(
                    OUTPUT
                      ${CMAKE_BINARY_DIR}/inc/nl_config.inc
                      ${CMAKE_BINARY_DIR}/frame/module_state_description.F
                    WORKING_DIRECTORY
                      ${CMAKE_BINARY_DIR}
                    # Replicate what exists in project directory for registry
                    COMMAND
                      ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/Registry
                    COMMAND
                      ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/inc
                    COMMAND
                      ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/frame
                    COMMAND
                      ${CMAKE_BINARY_DIR}/tools/registry ${REGISTRY_DEFS} -DNEW_BDYS ${REGISTRY_FILE} > ${CMAKE_BINARY_DIR}/registry.log 2>&1
                    #!TODO Just have the registry code actually check for failure or better yet rewrite the 
                    # registry code to not be so obfuscating
                    COMMAND 
                      ${CMAKE_COMMAND} -E compare_files ${CMAKE_BINARY_DIR}/inc/nl_config.inc ${CMAKE_BINARY_DIR}/inc/nl_config.inc
                    DEPENDS
                      ${FOLDER_COMPILE_TARGET}
                    )

if ( ${USE_MPI} )

  set( NUM_SPLIT_COMM_NESTING 4 )
  math( EXPR NUM_SPLIT_COMM_NESTING_LOOP "${NUM_SPLIT_COMM_NESTING} - 1" )

  set( COMM_NESTING_DECL ${CMAKE_BINARY_DIR}/inc/REGISTRY_COMM_NESTING_DM_subs_decl.inc )
  set( COMM_NESTING_SUBS      )
  foreach( n RANGE 0 ${NUM_SPLIT_COMM_NESTING_LOOP} )
    list( APPEND COMM_NESTING_SUBS ${CMAKE_BINARY_DIR}/inc/REGISTRY_COMM_NESTING_DM_subs_${n}_def.F )
  endforeach()

  add_custom_command(
                      OUTPUT
                        ${COMM_NESTING_SUBS}
                        ${COMM_NESTING_DECL}
                      WORKING_DIRECTORY
                        ${CMAKE_BINARY_DIR}
                      # Replicate what exists in project directory for registry
                      COMMAND
                      # split across a maximum of 4 files
                        ${PROJECT_SOURCE_DIR}/arch/subroutine_split.py 
                          module_comm_nesting_dm                                        # module
                          ${CMAKE_BINARY_DIR}/inc/REGISTRY_COMM_NESTING_DM_subs.inc     # source
                          ${CMAKE_BINARY_DIR}/inc/                                      # ouptut dir
                          ${PROJECT_SOURCE_DIR}/arch/submodule.template.F               # template submod
                          ${NUM_SPLIT_COMM_NESTING}                                     # split amount
                      COMMAND 
                        ${CMAKE_COMMAND} -E compare_files ${COMM_NESTING_DECL} ${COMM_NESTING_DECL}
                      DEPENDS
                        ${CMAKE_BINARY_DIR}/inc/nl_config.inc
                        ${FOLDER_COMPILE_TARGET}
                      )

  set( NUM_SPLIT_COMM_DM 6 )
  math( EXPR NUM_SPLIT_COMM_DM_LOOP "${NUM_SPLIT_COMM_DM} - 1" )
  

  set( COMM_DM_PERIOD_DECL ${CMAKE_BINARY_DIR}/inc/REGISTRY_COMM_DM_PERIOD_subs_decl.inc )
  set( COMM_DM_PERIOD_SUBS       )
  foreach( n RANGE 0 ${NUM_SPLIT_COMM_DM_LOOP} )
    list( APPEND COMM_DM_PERIOD_SUBS ${CMAKE_BINARY_DIR}/inc/REGISTRY_COMM_DM_PERIOD_subs_${n}_def.F )
  endforeach()

  add_custom_command(
                      OUTPUT
                        ${COMM_DM_PERIOD_SUBS}
                        ${COMM_DM_PERIOD_DECL}
                      WORKING_DIRECTORY
                        ${CMAKE_BINARY_DIR}
                      # Replicate what exists in project directory for registry
                      COMMAND
                      # split across a maximum of 4 files
                        ${PROJECT_SOURCE_DIR}/arch/subroutine_split.py 
                          module_comm_dm_4_                                             # module
                          ${CMAKE_BINARY_DIR}/inc/REGISTRY_COMM_DM_PERIOD_subs.inc      # source
                          ${CMAKE_BINARY_DIR}/inc/                                      # ouptut dir
                          ${PROJECT_SOURCE_DIR}/arch/submodule.template.F               # template submod
                          ${NUM_SPLIT_COMM_DM}                                          # split amount
                      COMMAND 
                        ${CMAKE_COMMAND} -E compare_files ${COMM_DM_PERIOD_DECL} ${COMM_DM_PERIOD_DECL}
                      DEPENDS
                        ${CMAKE_BINARY_DIR}/inc/nl_config.inc
                        ${FOLDER_COMPILE_TARGET}
                      )

  set( NUM_COMM_DM_SUBS 4 )
  math( EXPR NUM_COMM_DM_SUBS_LOOP "${NUM_COMM_DM_SUBS} - 1" )
  set( COMM_DM_DECL )
  set( COMM_DM_SUBS )

  foreach( i RANGE 0 ${NUM_COMM_DM_SUBS_LOOP} )

    set( COMM_DM_DECL_${i}   ${CMAKE_BINARY_DIR}/inc/REGISTRY_COMM_DM_${i}_subs_decl.inc )
    set( COMM_DM_SUBS_${i}   )

    foreach( n RANGE 0 ${NUM_SPLIT_COMM_DM_LOOP} )
      list( APPEND COMM_DM_SUBS_${i} ${CMAKE_BINARY_DIR}/inc/REGISTRY_COMM_DM_${i}_subs_${n}_def.F )
    endforeach()

    list( APPEND COMM_DM_DECL ${COMM_DM_SUBS_${i}} )
    list( APPEND COMM_DM_SUBS ${COMM_DM_SUBS_${i}} )

    add_custom_command(
                        OUTPUT
                          ${COMM_DM_DECL_${i}}
                          ${COMM_DM_SUBS_${i}}
                        WORKING_DIRECTORY
                          ${CMAKE_BINARY_DIR}
                        # Replicate what exists in project directory for registry
                        COMMAND
                        # split across a maximum of 4 files
                          ${PROJECT_SOURCE_DIR}/arch/subroutine_split.py 
                            module_comm_dm_${i}                                           # module
                            ${CMAKE_BINARY_DIR}/inc/REGISTRY_COMM_DM_${i}_subs.inc        # source
                            ${CMAKE_BINARY_DIR}/inc/                                      # ouptut dir
                            ${PROJECT_SOURCE_DIR}/arch/submodule.template.F               # template submod
                            ${NUM_SPLIT_COMM_DM}                                          # split amount
                        COMMAND 
                          ${CMAKE_COMMAND} -E compare_files ${COMM_DM_DECL_${i}} ${COMM_DM_DECL_${i}}
                        DEPENDS
                          ${CMAKE_BINARY_DIR}/inc/nl_config.inc
                          ${FOLDER_COMPILE_TARGET}
                        )
  endforeach()

  set( COMM_NESTING_SUBS   ${COMM_NESTING_SUBS}   PARENT_SCOPE )
  set( COMM_DM_PERIOD_SUBS ${COMM_DM_PERIOD_SUBS} PARENT_SCOPE )
  set( COMM_DM_SUBS        ${COMM_DM_SUBS}        PARENT_SCOPE )

endif()


add_custom_target( 
                  registry_code
                  DEPENDS
                    ${CMAKE_BINARY_DIR}/inc/nl_config.inc
                    ${CMAKE_BINARY_DIR}/frame/module_state_description.F
                    ${COMM_NESTING_DECL}
                    ${COMM_DM_PERIOD_DECL}
                    ${COMM_DM_DECL}
                  )

