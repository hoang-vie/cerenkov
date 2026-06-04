#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 3.22)
project(CherenkovWater)

#----------------------------------------------------------------------------
# Find Geant4 package
# Set WITH_GEANT4_UIVIS to OFF for batch-only mode
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include)

#----------------------------------------------------------------------------
# ROOT
#
find_package(ROOT REQUIRED)
include(${ROOT_USE_FILE})

#----------------------------------------------------------------------------
# CRY cosmic-ray generator
#
include(ExternalProject)
ExternalProject_Add(
    CRY-1.7
    URL http://nuclear.llnl.gov/simulation/cry_v1.7.tar.gz
    URL_MD5 85f240bebe81fe0b257e92ef1d390a83
    CONFIGURE_COMMAND ""
    INSTALL_COMMAND  ""
    BUILD_IN_SOURCE  1
    TEST_AFTER_INSTALL 1
)

ExternalProject_Get_Property(CRY-1.7 source_dir)
set(CRY_INCLUDE_DIRS "${source_dir}/src")
set(CRY_LIBRARY_DIR  "${source_dir}/lib")
set(CRY_LIBRARIES    "libCRY.a")
include_directories(${CRY_INCLUDE_DIRS})
link_directories(${CRY_LIBRARY_DIR})

#----------------------------------------------------------------------------
# Sources and headers
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Main executable
#
add_executable(cherenkov cherenkov.cc ${sources} ${headers})
add_dependencies(cherenkov CRY-1.7)

target_link_libraries(cherenkov
    ${Geant4_LIBRARIES}
    ${CRY_LIBRARIES}
    ${ROOT_LIBRARIES}
)

#----------------------------------------------------------------------------
# Copy macros to build directory
#
file(GLOB RUN_MACROS macros/*.mac)
foreach(_script ${RUN_MACROS})
  get_filename_component(shortname ${_script} NAME)
  configure_file(
    macros/${shortname}
    ${PROJECT_BINARY_DIR}/${shortname}
    COPYONLY
  )
endforeach()


#----------------------------------------------------------------------------
add_custom_target(CherenkovWater DEPENDS cherenkov)

install(TARGETS cherenkov DESTINATION bin)