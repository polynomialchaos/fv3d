# - Try to find METIS
# Once done this will define
#  METIS_FOUND - System has METIS
#  METIS_INCLUDE_DIRS - The METIS include directories
#  METIS_C_LIBRARIES - The libraries needed to use METIS

find_path(
    METIS_INCLUDE_DIR "metis.h"
    HINTS "${METIS_DIR}"
    PATH_SUFFIXES "include"
)

find_library(
    METIS_LIBRARY NAMES "libmetis.${suffix}"
    HINTS "${METIS_DIR}"
    PATH_SUFFIXES "lib" "lib/x86_64-linux-gnu"
)

include( FindPackageHandleStandardArgs )
# handle the QUIETLY and REQUIRED arguments and set METIS_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(
    METIS DEFAULT_MSG
    METIS_LIBRARY
    METIS_INCLUDE_DIR
)

mark_as_advanced( METIS_LIBRARY )
mark_as_advanced( METIS_INCLUDE_DIR )

set( METIS_C_LIBRARIES ${METIS_LIBRARY} )
set( METIS_C_INCLUDE_PATH ${METIS_INCLUDE_DIR} )
set( METIS_C_COMPILE_FLAGS "" )
set( METIS_C_LINK_FLAGS "" )