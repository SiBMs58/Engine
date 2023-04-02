cmake_minimum_required(VERSION 3.5)

project(Engine)

############################################################
# Set compiler flags
############################################################
# If supported by your compiler, you can add the -Wall, -Wextra, –fstack-protector-all and -g3 flags here.
set(OWN_GXX_FLAGS "-std=c++14")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${OWN_GXX_FLAGS}")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${OWN_GXX_FLAGS}")
set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${OWN_GXX_FLAGS}")

############################################################
# List all sources
############################################################
set(engine_sources easy_image.cc engine.cc ini_configuration.cc l_parser.cc vector3d.cc
        2DLsystem/2DLine/Color.cpp 2DLsystem/2DLine/Point2D.cpp 2DLsystem/2DLine/Line2D.cpp 2DLsystem/Lsystem2D.cpp
        3DLsystem/3DFigure/Figure.cpp 3DLsystem/3DFigure/Face.cpp 3DLsystem/Lsystem.cpp 3DLsystem/Wireframe3D.cpp)

############################################################
# Create an executable
############################################################
set(exe_name "engine")
add_executable( ${exe_name} ${engine_sources} )
install( TARGETS ${exe_name} DESTINATION ${PROJECT_SOURCE_DIR}/ )
