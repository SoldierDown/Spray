set(ENABLE_SPRAY OFF CACHE BOOL "Enable building the Spray project tree.")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake/")
cmake_policy(SET CMP0054 NEW)

if(ENABLE_NON_FICKS_DIFFUSION)
add_subdirectory(spray_2d)
# add_subdirectory(spray_3d)
add_subdirectory(plugin_Spray)
endif(ENABLE_NON_FICKS_DIFFUSION)
