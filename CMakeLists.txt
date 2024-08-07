cmake_minimum_required(VERSION 3.12)

#--------------------------------------------------------------------
# List of available command-line options:
#   -DDEBUG_MODE:BOOL=<ON or OFF>
#   -DROOT_CONFIG_DIR:PATH=<root-config directory>
#   -DBAT_CONFIG_DIR:PATH=<bat-config directory>
#   -DCMAKE_BUILD_TYPE:STRING=<Debug or Release>
#   -DCMAKE_INSTALL_PREFIX:PATH=<installation directory>
#--------------------------------------------------------------------

if(DEBUG_MODE)
  set(CMAKE_VERBOSE_MAKEFILE ON)
endif(DEBUG_MODE)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

project(GammaDDbar)

##########  ROOT  ##########

if(NOT ROOT_CONFIG_DIR)
  find_program(ROOT_CONFIG root-config DOC "root-config")
  if(ROOT_CONFIG)
    string(REGEX REPLACE "/root-config$" "" TMP_PATH ${ROOT_CONFIG})
    set(ROOT_CONFIG_DIR ${TMP_PATH} CACHE PATH "Path to root-config")
    unset(TMP_PATH)
  endif()
else()
  set(ROOT_CONFIG_DIR ${ROOT_CONFIG_DIR} CACHE PATH "Path to root-config" FORCE)
  find_program(ROOT_CONFIG root-config PATHS ${ROOT_CONFIG_DIR} NO_DEFAULT_PATH DOC "root-config")
endif()
if(NOT ROOT_CONFIG)
  unset(ROOT_CONFIG CACHE)
  message(FATAL_ERROR "root-config could not be found. Please specify the directory of root-config with ${RedBold}-DROOT_CONFIG_DIR=<directory>${ColReset}.")
else()
  message(STATUS "Found root-config in ${ROOT_CONFIG_DIR}")
  execute_process(COMMAND ${ROOT_CONFIG} --version
                  OUTPUT_VARIABLE ROOT_VERSION
                  OUTPUT_STRIP_TRAILING_WHITESPACE)
  message(STATUS "ROOT version: ${ROOT_VERSION}")
  execute_process(COMMAND ${ROOT_CONFIG} --cflags
                  OUTPUT_VARIABLE ROOT_CFLAGS
                  OUTPUT_STRIP_TRAILING_WHITESPACE)
  message(STATUS "ROOT cflags: ${ROOT_CFLAGS}")
  execute_process(COMMAND ${ROOT_CONFIG} --libs
                  OUTPUT_VARIABLE ROOT_LIBS
                  OUTPUT_STRIP_TRAILING_WHITESPACE)
  message(STATUS "ROOT libs: ${ROOT_LIBS}")
  execute_process(COMMAND ${ROOT_CONFIG} --ldflags
                  OUTPUT_VARIABLE ROOT_LDFLAGS
                  OUTPUT_STRIP_TRAILING_WHITESPACE)
  message(STATUS "ROOT ldflags: ${ROOT_LDFLAGS}")
  execute_process(COMMAND ${ROOT_CONFIG} --libdir
                  OUTPUT_VARIABLE ROOT_LIB_PATH
                  OUTPUT_STRIP_TRAILING_WHITESPACE)
  string(REGEX REPLACE "/$" "" ROOT_SYS ${ROOT_CONFIG_DIR})
  string(REGEX REPLACE "/bin$" "" ROOT_SYS ${ROOT_SYS})
  unset(ROOT_CONFIG CACHE)
endif()

##########  BAT  ##########

if(NOT BAT_CONFIG_DIR)
  find_program(BAT_CONFIG bat-config DOC "bat-config")
  if(BAT_CONFIG)
    string(REGEX REPLACE "/bat-config$" "" TMP_PATH ${BAT_CONFIG})
    set(BAT_CONFIG_DIR ${TMP_PATH} CACHE PATH "Path to bat-config")
    unset(TMP_PATH)
  endif()
else()
  set(BAT_CONFIG_DIR ${BAT_CONFIG_DIR} CACHE PATH "Path to bat-config" FORCE)
  find_program(BAT_CONFIG bat-config PATHS ${BAT_CONFIG_DIR} NO_DEFAULT_PATH DOC "bat-config")
endif()
if(NOT BAT_CONFIG)
  unset(BAT_CONFIG CACHE)
  message(FATAL_ERROR "bat-config could not be found. Please specify the directory of bat-config with ${RedBold}-DBAT_CONFIG_DIR=<directory>${ColReset}.")
else()
  message(STATUS "Found bat-config in ${BAT_CONFIG_DIR}")
  execute_process(COMMAND ${BAT_CONFIG} --cflags
                  OUTPUT_VARIABLE BAT_CFLAGS
                  OUTPUT_STRIP_TRAILING_WHITESPACE)
  message(STATUS "BAT cflags: ${BAT_CFLAGS}")
  execute_process(COMMAND ${BAT_CONFIG} --libs
                  OUTPUT_VARIABLE BAT_LIBS
                  OUTPUT_STRIP_TRAILING_WHITESPACE)
  message(STATUS "BAT libs: ${BAT_LIBS}")
  execute_process(COMMAND ${BAT_CONFIG} --libdir
                  OUTPUT_VARIABLE BAT_LIB_PATH
                  OUTPUT_STRIP_TRAILING_WHITESPACE)
  message(STATUS "BAT libdir: ${BAT_LIB_PATH}")
  execute_process(COMMAND ${BAT_CONFIG} --incdir
                  OUTPUT_VARIABLE BAT_INCLUDE_PATH
                  OUTPUT_STRIP_TRAILING_WHITESPACE)
  message(STATUS "BAT incdir: ${BAT_INCLUDE_PATH}")
  unset(BAT_CONFIG CACHE)
endif()

set(CMAKE_CXX_FLAGS "-Wall ${BAT_CFLAGS} ${ROOT_CFLAGS}")
set(CMAKE_CXX_LINK_FLAGS "")
set(CMAKE_CXX_LINK_FLAGS "${CMAKE_CXX_LINK_FLAGS} -L${BAT_LIB_PATH}")

# Add the source files in the "Codes" directory
file(GLOB SOURCES "Codes/*.cpp")

link_directories(${BAT_LIB})

# Create the executable
add_executable(${PROJECT_NAME} ${SOURCES})

# Link against the libraries
target_link_libraries(${PROJECT_NAME} ${BAT_LIBS})
target_link_libraries(${PROJECT_NAME} ${ROOT_LIBS})


INSTALL(TARGETS GammaDDbar DESTINATION bin COMPONENT executable)