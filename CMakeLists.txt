cmake_minimum_required(VERSION 3.18)
project(pysme_astro LANGUAGES C CXX Fortran)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Find Python and NumPy
find_package(Python REQUIRED COMPONENTS Interpreter Development.Module NumPy)

# Generate platform.h
set(PLATFORM "${CMAKE_SYSTEM_NAME}")
string(TOLOWER "${PLATFORM}" PLATFORM)
configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/platform.h.in
    ${CMAKE_CURRENT_BINARY_DIR}/platform.h
    @ONLY
)

# SMElib sources
set(SMELIB_SOURCES
    smelib/src/sme/sme_synth_faster.cpp
    smelib/src/eos/eos.f
    smelib/src/eos/eos_eqns.f
    smelib/src/eos/eos_math_special.f
    smelib/src/eos/kp_q_spln.f
    smelib/src/sme/hlinop.f
    smelib/src/sme/hlinprof.f
)

# Build libsme as shared library
add_library(sme SHARED ${SMELIB_SOURCES})
target_include_directories(sme PRIVATE
    ${CMAKE_CURRENT_BINARY_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/smelib/src/sme
)
target_compile_options(sme PRIVATE -O3)

if(APPLE)
    set_target_properties(sme PROPERTIES
        INSTALL_RPATH "@loader_path"
        BUILD_WITH_INSTALL_RPATH TRUE
        INSTALL_NAME_DIR "@rpath"
        MACOSX_RPATH TRUE
    )
elseif(UNIX)
    set_target_properties(sme PROPERTIES
        INSTALL_RPATH "$ORIGIN"
        BUILD_WITH_INSTALL_RPATH TRUE
    )
endif()

# Build Python extension module
Python_add_library(_smelib MODULE smelib/pymodule/_smelib.cpp)
target_include_directories(_smelib PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/smelib/src/sme
    ${Python_NumPy_INCLUDE_DIRS}
)
target_link_libraries(_smelib PRIVATE sme)

if(APPLE)
    set_target_properties(_smelib PROPERTIES
        INSTALL_RPATH "@loader_path/lib"
        BUILD_WITH_INSTALL_RPATH TRUE
    )
elseif(UNIX)
    set_target_properties(_smelib PROPERTIES
        INSTALL_RPATH "$ORIGIN/lib"
        BUILD_WITH_INSTALL_RPATH TRUE
    )
endif()

# Install targets
install(TARGETS sme
    LIBRARY DESTINATION ${SKBUILD_PLATLIB_DIR}/pysme/smelib/lib
    RUNTIME DESTINATION ${SKBUILD_PLATLIB_DIR}/pysme/smelib/lib
)

install(TARGETS _smelib
    LIBRARY DESTINATION ${SKBUILD_PLATLIB_DIR}/pysme/smelib
)

# Install smelib data files
install(DIRECTORY smelib/src/data/
    DESTINATION ${SKBUILD_PLATLIB_DIR}/pysme/smelib/data
)

# Install pysme data files
install(FILES
    src/pysme/config_default.json
    src/pysme/datafiles_atmospheres.json
    src/pysme/datafiles_nlte.json
    src/pysme/lineprof.dat.gz
    DESTINATION ${SKBUILD_PLATLIB_DIR}/pysme
)
