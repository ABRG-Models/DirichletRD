# AUTOMOC requires 2.8.6. Per-target include directories would require 2.8.8
cmake_minimum_required(VERSION 2.8.8)
project(Ermentrout2009)

message(STATUS "Install prefix: ${CMAKE_INSTALL_PREFIX}")
message(STATUS "  (This can be changed with `cmake -DCMAKE_INSTALL_PREFIX=/some/place`")

# From CMAKE_SYSTEM work out which of __OSX__, __GLN__, __NIX__ or
# __WIN__ are required
message(STATUS "Operating system: " ${CMAKE_SYSTEM})
if(CMAKE_SYSTEM MATCHES Linux.*)
  set(MORPH_HOST_DEFINITION "-D__GLN__")
elseif(CMAKE_SYSTEM MATCHES BSD.*)
  set(MORPH_HOST_DEFINITION "-D__NIX__")
elseif(APPLE)
  set(MORPH_HOST_DEFINITION "-D__OSX__")
elseif(CMAKE_SYSTEM MATCHES Win.*)
  set(MORPH_HOST_DEFINITION "-D__WIN__")
else()
  message(ERROR "Operating system not supported: " ${CMAKE_SYSTEM})
endif()

# Add the host definition to CXXFLAGS along with other switches
set(CMAKE_CXX_FLAGS "${MORPH_HOST_DEFINITION} -Wall -std=c++11")

# Lib finding - we have to do this, so that the correct include paths
# are present for OpenGL, OpenCV, etc; they're all included by
# morphologica header files.
find_package(OpenCV REQUIRED)
find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)
find_package(X11 REQUIRED)
include_directories(${OpenCV_INCLUDE_DIRS} ${OPENGL_INCLUDE_DIR} ${GLUT_INCLUDE_DIR} ${X11_INC_SEARCH_PATH})

# Find libmorphologica using pkgconfig
find_package(PkgConfig)
if(PKG_CONFIG_FOUND)
  pkg_check_modules(MORPH REQUIRED libmorphologica)
  if (MORPH_FOUND)
    include_directories(${MORPH_INCLUDEDIR})
  else()
    message(FATAL_ERROR "You need libmorphologica from github.com/ABRG_Models/morphologica")
  endif(MORPH_FOUND)
else()
  message(WARNING "There's no pkg-config on this system to check for libmorphologica. "
    "You can set MORPH_LIB with `cmake -DMORPH_LIB=/pathto/lib` ...")
endif()

# The main process executable
add_executable(process process.cpp)
target_link_libraries(process ${MORPH_LDFLAGS})

# For debugging of variables:
set(DEBUG_VARIABLES ON)
if(DEBUG_VARIABLES)
  get_cmake_property(_variableNames VARIABLES)
  foreach (_variableName ${_variableNames})
    message(STATUS "${_variableName}=${${_variableName}}")
  endforeach()
endif(DEBUG_VARIABLES)
