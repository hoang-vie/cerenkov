cmake_minimum_required(VERSION 3.22)
project(CherenkovGeneral LANGUAGES CXX)

option(WITH_GEANT4_UIVIS "Build with Geant4 UI and visualization drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()
include(${Geant4_USE_FILE})

# Build/source provenance.  These values are embedded into RunConfig so every
# ROOT output identifies the exact source revision used to build the executable.
# The fallback values keep source tarballs/builds without a .git directory valid.
find_package(Git QUIET)
set(CHERENKOV_GIT_COMMIT "unknown")
set(CHERENKOV_GIT_BRANCH "unknown")
set(CHERENKOV_GIT_DIRTY 0)
if(Git_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
  execute_process(
    COMMAND ${GIT_EXECUTABLE} rev-parse HEAD
    WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
    OUTPUT_VARIABLE _cherenkov_git_commit
    OUTPUT_STRIP_TRAILING_WHITESPACE
    RESULT_VARIABLE _cherenkov_git_commit_rc
    ERROR_QUIET
  )
  if(_cherenkov_git_commit_rc EQUAL 0 AND NOT "${_cherenkov_git_commit}" STREQUAL "")
    set(CHERENKOV_GIT_COMMIT "${_cherenkov_git_commit}")
  endif()

  execute_process(
    COMMAND ${GIT_EXECUTABLE} rev-parse --abbrev-ref HEAD
    WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
    OUTPUT_VARIABLE _cherenkov_git_branch
    OUTPUT_STRIP_TRAILING_WHITESPACE
    RESULT_VARIABLE _cherenkov_git_branch_rc
    ERROR_QUIET
  )
  if(_cherenkov_git_branch_rc EQUAL 0 AND NOT "${_cherenkov_git_branch}" STREQUAL "")
    set(CHERENKOV_GIT_BRANCH "${_cherenkov_git_branch}")
  endif()

  execute_process(
    COMMAND ${GIT_EXECUTABLE} status --porcelain --untracked-files=all
    WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
    OUTPUT_VARIABLE _cherenkov_git_status
    OUTPUT_STRIP_TRAILING_WHITESPACE
    RESULT_VARIABLE _cherenkov_git_status_rc
    ERROR_QUIET
  )
  if(_cherenkov_git_status_rc EQUAL 0 AND NOT "${_cherenkov_git_status}" STREQUAL "")
    set(CHERENKOV_GIT_DIRTY 1)
  endif()
endif()

set(CHERENKOV_PROJECT_REVISION "${CHERENKOV_GIT_COMMIT}")
if(NOT CHERENKOV_GIT_COMMIT STREQUAL "unknown")
  string(SUBSTRING "${CHERENKOV_GIT_COMMIT}" 0 12 CHERENKOV_PROJECT_REVISION)
endif()
if(CHERENKOV_GIT_DIRTY)
  string(APPEND CHERENKOV_PROJECT_REVISION "-dirty")
endif()

configure_file(
  "${PROJECT_SOURCE_DIR}/include/BuildInfo.hh.in"
  "${PROJECT_BINARY_DIR}/generated/BuildInfo.hh"
  @ONLY
)

include(ExternalProject)
ExternalProject_Add(
  CRY-1.7
  URL https://nuclear.llnl.gov/simulation/cry_v1.7.tar.gz
  URL_MD5 85f240bebe81fe0b257e92ef1d390a83
  CONFIGURE_COMMAND ""
  INSTALL_COMMAND ""
  BUILD_IN_SOURCE 1
  TEST_AFTER_INSTALL 1
)
ExternalProject_Get_Property(CRY-1.7 source_dir)
set(CRY_INCLUDE_DIRS "${source_dir}/src")
set(CRY_LIBRARY_DIR "${source_dir}/lib")

file(GLOB sources CONFIGURE_DEPENDS ${PROJECT_SOURCE_DIR}/src/*.cc)
list(FILTER sources EXCLUDE REGEX ".*/cherenkov\\.cc$")
file(GLOB headers CONFIGURE_DEPENDS ${PROJECT_SOURCE_DIR}/include/*.hh)

add_executable(cherenkov cherenkov.cc ${sources} ${headers})
add_dependencies(cherenkov CRY-1.7)
target_compile_features(cherenkov PRIVATE cxx_std_17)
target_include_directories(cherenkov PRIVATE
  ${PROJECT_SOURCE_DIR}/include
  ${PROJECT_BINARY_DIR}/generated
  ${CRY_INCLUDE_DIRS}
)
target_link_directories(cherenkov PRIVATE ${CRY_LIBRARY_DIR})
target_link_libraries(cherenkov PRIVATE ${Geant4_LIBRARIES} libCRY.a)
target_compile_definitions(cherenkov PRIVATE
  CRY_DATA_DIR="${source_dir}/data"
  OPTICAL_DATA_DIR="${PROJECT_SOURCE_DIR}/data"
)

file(GLOB RUN_MACROS CONFIGURE_DEPENDS
     "${PROJECT_SOURCE_DIR}/macros/*.mac")

foreach(_script IN LISTS RUN_MACROS)
  get_filename_component(_name "${_script}" NAME)
  set(_dst "${PROJECT_BINARY_DIR}/${_name}")
  file(REMOVE "${_dst}")
  file(CREATE_LINK "${_script}" "${_dst}" SYMBOLIC)
endforeach()

install(TARGETS cherenkov DESTINATION bin)
