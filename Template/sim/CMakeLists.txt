# Ermentrout (2009) model
add_executable(ermentrout ermentrout.cpp)
target_link_libraries(ermentrout ${MORPH_LDFLAGS} ${HDF5_C_${LIB_TYPE}_LIBRARY})
