# 3.1 minimum version for set(CMAKE_CXX_STANDARD 1z)
cmake_minimum_required(VERSION 3.1)
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

set(CMAKE_CXX_STANDARD 17)

# Add the host definition to CXXFLAGS along with other switches
if (APPLE)
  # OpenMP: According to
  # https://stackoverflow.com/questions/39979836/using-openmp-with-c11-on-mac-os#39983933
  # it should now be a case of:
  #
  # You should be able to get away without the libomp that's mentioned
  # in that post, as I don't use any of the runtime library OpenMP
  # features. However, the code did not compile without it, so we've moved the links to omp to the later apple if case.
  #
  set(CMAKE_CXX_FLAGS "${MORPH_HOST_DEFINITION} -Wall -g -Xpreprocessor -fopenmp")

else()
  # To use Intel compiler, call cmake as: cmake -DCMAKE_CXX_COMPILER=icpc ..
  if (CMAKE_CXX_COMPILER_ID MATCHES Intel)
    set(CMAKE_CXX_FLAGS "${MORPH_HOST_DEFINITION} -Wall -g -std=c++17 -xHOST -O3 -qopenmp -D__ICC__")
  else() # GCC or Clang
    set(CMAKE_CXX_FLAGS "${MORPH_HOST_DEFINITION} -Wall -Wfatal-errors -g -Wno-unused-result -Wno-unknown-pragmas -march=native -O3 -fopenmp")
  endif()
endif()

# Tell clang to be quiet about brace initialisers
if(CMAKE_CXX_COMPILER_ID MATCHES Clang)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-missing-braces")
endif()

# Additional GL compiler flags.
#
# Following `cmake --help-policy CMP0072`
set(OpenGL_GL_PREFERENCE "GLVND")
# maybe -DGL_GLEXT_PROTOTYPES
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DGL3_PROTOTYPES -DGL_GLEXT_PROTOTYPES")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DGL3_PROTOTYPES -DGL_GLEXT_PROTOTYPES")
# An option that's useful for Ubuntu 16.04 builds
option(USE_GLEW "Link libglew.so (try if the linker can't find glCreateVertexArrays)" OFF)


# Lib finding - we have to do this, so that the correct include paths
# are present for OpenGL, OpenCV, etc; they're all included by
# morphologica header files.
find_package(OpenCV REQUIRED)
find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)
find_package(X11 REQUIRED)


if(APPLE)
  find_package (HDF5 REQUIRED)
  set (HDF5LIBLINK ${HDF5_C_LIBRARY_hdf5})
  set (JSONLIBLINK /usr/local/lib/libjsoncpp.dylib)
  include_directories(/opt/local/include/libomp)
  link_libraries(-L/opt/local/lib/libomp omp)
else()
  set (LIB_TYPE SHARED) # STATIC or SHARED
  string(TOLOWER ${LIB_TYPE} SEARCH_TYPE)
  find_package (HDF5 NAMES hdf5 COMPONENTS C ${SEARCH_TYPE} REQUIRED)
  set (HDF5LIBLINK ${HDF5_C_${LIB_TYPE}_LIBRARY})
  set (JSONLIBLINK jsoncpp)
endif(APPLE)


# glfw3 library finding. However, it's possible to compile from
# source, and there are instructions on the glfw3 website on how to do
# this with cmake.
#
# glfw3 provides a cmake find_package script, though it doesn't set
# lib paths etc. Then you just link with "glfw3" if you happened to
# compile the static libglfw3.a or with "glfw" if you compiled the
# dynamic libglfw.so
find_package(glfw3 QUIET)
# You can find GLFW library with pkgconfig, but it's switched off for now
if(0)
  find_package(PkgConfig)
  if(PKG_CONFIG_FOUND)
    # Use pkg-config to check GLFW3 is present:
    pkg_check_modules(GLFW3 REQUIRED glfw3)
    if (GLFW3_FOUND)
      # Good....
      message(INFO, "pkg-config found glfw3. Excellent.")
    else(GLFW3_FOUND)
      # Bad....
      message(ERROR, "glfw3 was NOT found by pkg-config. Uh oh.")
    endif(GLFW3_FOUND)
  else()
    message(WARNING "There's no pkg-config on this system to check for libglfw3. You may need to try `sudo apt-get install pkg-config`")
  endif()
endif()

if(USE_GLEW)
  find_package(GLEW REQUIRED)
endif(USE_GLEW)

# pkgconfig is used to find JSON adn also to find lib paths for glfw3.
find_package(PkgConfig)

#include_directories(${OpenCV_INCLUDE_DIRS} ${OPENGL_INCLUDE_DIR} ${GLUT_INCLUDE_DIR} ${X11_INC_SEARCH_PATH} ${HDF5_INCLUDE_DIR})
include_directories(${OpenCV_INCLUDE_DIRS} ${OPENGL_INCLUDE_DIR} ${GLUT_INCLUDE_DIR} ${X11_INC_SEARCH_PATH} ${HDF5_INCLUDE_DIR} ${ARMADILLO_INCLUDE_DIR} ${GLFW3_INCLUDE_DIR} ${JSONCPP_INCLUDE_DIRS})

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


# JSON
if(${PKG_CONFIG_FOUND})
  message(STATUS "pkg-config was found. Search for JSON.")
  pkg_search_module(JSONCPP REQUIRED jsoncpp)
  if(${JSONCPP_FOUND})
    message(STATUS "pkg-config found jsoncpp. Excellent.")
  else(${JSONCPP_FOUND})
    message(WARNING ": jsoncpp was NOT found by pkg-config. Try default: /usr/local/lib/jsoncpp.a")
    set(JSONCPP_LDFLAGS "/usr/local/lib/libjsoncpp.a")
    set(JSONCPP_INCLUDE_DIRS "/usr/local/include")
  endif(${JSONCPP_FOUND})
else(${PKG_CONFIG_FOUND})
  message(STATUS "No pkg-config. Tryign default /usr/local/lib/libjsoncpp.a for jsoncpp.")
  set(JSONCPP_LDFLAGS "/usr/local/lib/libjsoncpp.a")
  set(JSONCPP_INCLUDE_DIRS "/usr/local/include")
endif(${PKG_CONFIG_FOUND})


# Find Armadillo using find_package
find_package(Armadillo REQUIRED)

# The main process executable
add_executable(process process.cpp)
target_link_libraries(process ${MORPH_LDFLAGS} ${HDF5LIBLINK} ${JSONLIBLINK} ${ARMADILLO_LIBRARIES})

#use online display
add_executable(processDisplay processDisplay.cpp)
target_compile_definitions(processDisplay PUBLIC FLT=double COMPILE_PLOTTING)
target_link_libraries(processDisplay ${MORPH_LDFLAGS} ${HDF5LIBLINK} ${JSONLIBLINK} ${ARMADILLO_LIBRARIES} ${OPENGL_LIBRARIES})

#run without display, e.g. remotely
add_executable(processNoDisplay processNoDisplay.cpp)
target_compile_definitions(processNoDisplay PUBLIC FLT=double)
target_link_libraries(processNoDisplay ${MORPH_LDFLAGS} ${HDF5LIBLINK} ${JSONLIBLINK} ${ARMADILLO_LIBRARIES})

#run just to create display from noDisplay run
add_executable(processField processField.cpp)
target_compile_definitions(processField PUBLIC FLT=double COMPILE_PLOTTING)
target_link_libraries(processField ${MORPH_LDFLAGS} ${HDF5LIBLINK} ${JSONLIBLINK} morphdisplay0 armadillo)


#use json for parameters
add_executable(processJson processJson.cpp)
target_compile_definitions(processJson PUBLIC FLT=double COMPILE_PLOTTING)
target_link_libraries(processJson ${MORPH_LDFLAGS} ${HDF5LIBLINK} ${JSONLIBLINK} ${ARMADILLO_LIBRARIES})

#to create random centres
add_executable(setCentres setCentres.cpp)
target_compile_definitions(setCentres PUBLIC FLT=double COMPILE_PLOTTING)
target_link_libraries(setCentres ${MORPH_LDFLAGS} ${HDF5LIBLINK})

#to create random centres
add_executable(testCorrelateVector testCorrelateVector.cpp)
target_compile_definitions(testCorrelateVector PUBLIC FLT=double COMPILE_PLOTTING)
target_link_libraries(testCorrelateVector ${MORPH_LDFLAGS} ${HDF5LIBLINK} armadillo)

#to create random centres
#add_executable(honda honda.cpp)
#target_compile_definitions(honda PUBLIC FLT=double)
#target_link_libraries(honda ${MORPH_LDFLAGS})

add_executable(honda honda.cpp)
target_compile_definitions(honda PUBLIC FLT=float)
target_link_libraries(honda ${MORPH_LDFLAGS})

#to resize centres
add_executable(transform110 transform110.cpp)
target_compile_definitions(transform110 PUBLIC FLT=double COMPILE_PLOTTING)
target_link_libraries(transform110 ${MORPH_LDFLAGS})


#to test Visual
add_executable(testVisual testVisual.cpp)
target_compile_definitions(testVisual PUBLIC FLT=double COMPILE_PLOTTING)
#target_link_libraries(testVisual ${MORPH_LDFLAGS} ${HDF5LIBLINK} ${JSONLIBLINK} ${X11_LIBRARIES} GL GLU glut morphologica armadillo glfw3)
target_link_libraries(testVisual ${MORPH_LDFLAGS} ${OPENGL_LIBRARIES} ${HDF5LIBLINK} ${JSONLIBLINK} ${X11_LIBRARIES}  morphologica armadillo glfw3)

# For debugging of variables:
set(DEBUG_VARIABLES ON)
if(DEBUG_VARIABLES)
  get_cmake_property(_variableNames VARIABLES)
  foreach (_variableName ${_variableNames})
    message(STATUS "${_variableName}=${${_variableName}}")
  endforeach()
endif(DEBUG_VARIABLES)
