cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(ab_sim)

set(exec_name "main")
set(3rdparty "3rdparty")

option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)

find_package(Geant4 REQUIRED ui_all vis_all)

include(${Geant4_USE_FILE})

find_package(ROOT REQUIRED)

include_directories(${PROJECT_SOURCE_DIR}/include
    ${Geant4_INCLUDE_DIR}
    ${ROOT_INCLUDE_DIRS}
    ${ROOT_INCLUDE_DIR}
    ${3rdparty}
    )

file(GLOB sources 
    ${PROJECT_SOURCE_DIR}/src/*cxx
    ${PROJECT_SOURCE_DIR}/src/*cpp
    )
file(GLOB headers 
  ${PROJECT_SOURCE_DIR}/include/*hpp
  ${PROJECT_SOURCE_DIR}/include/*h
  ${PROJECT_SOURCE_DIR}/include/*hh
  )

add_executable(${exec_name} ${sources} ${headers})
target_link_libraries(${exec_name} ${Geant4_LIBRARIES} ${ROOT_LIBRARIES})
