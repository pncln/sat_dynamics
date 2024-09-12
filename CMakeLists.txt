cmake_minimum_required(VERSION 3.16)

project(sat_dynamics LANGUAGES CXX)

set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(QT NAMES Qt6 Qt5 REQUIRED COMPONENTS Core Widgets)
find_package(Qt${QT_VERSION_MAJOR} REQUIRED COMPONENTS Core Widgets)
find_package(Qt6 REQUIRED COMPONENTS
    3DCore
    3DRender
    3DInput
    3DLogic
    3DExtras
    3DAnimation
)

add_executable(sat_dynamics
  main.cpp
  SatelliteModel.h SatelliteModel.cpp
  AttitudeControlSystem.h AttitudeControlSystem.cpp
  Constants.h
  functions/keplerian_to_cartesian.h functions/keplerian_to_cartesian.cpp
  functions/cartesian_to_keplerian.h functions/cartesian_to_keplerian.cpp
  functions/cross_product.h functions/cross_product.cpp
  functions/magnitude.h functions/magnitude.cpp
  functions/dot_product.h functions/dot_product.cpp
  functions/subtract_vectors.h functions/subtract_vectors.cpp
  functions/scale_vector.h functions/scale_vector.cpp
  functions/multiply_vector.h functions/multiply_vector.cpp
  functions/add_vectors.h functions/add_vectors.cpp
)
target_link_libraries(sat_dynamics Qt${QT_VERSION_MAJOR}::Core
    Qt6::3DCore
    Qt6::3DRender
    Qt6::3DLogic
    Qt6::3DExtras
    Qt6::3DAnimation)

include(GNUInstallDirs)
install(TARGETS sat_dynamics
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
