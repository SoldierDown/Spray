set(ENABLE_SPRAY OFF CACHE BOOL "Enable Spray.")
if(ENABLE_SPRAY)
  find_package(assimp)

  if(NOT assimp_FOUND)
    message(FATAL_ERROR "Cannot build embedded deformables plugin without the Assimp library being present.")
  endif()

  add_library(plugin_Spray
    ../../Nova_Examples/opengl/plugins/Assimp/AssimpRenderable_Model.cpp
    ../../Nova_Examples/opengl/plugins/Assimp/AssimpRenderable_Mesh.cpp
    main.cpp
    )

  if(WIN32)
    target_compile_definitions(plugin_Spray PUBLIC OPENGL_3D_WIN32=1)
  endif()
  
  if(UNIX)
    target_compile_definitions(plugin_Spray PUBLIC OPENGL_3D_LINUX=1)
  endif()

  target_include_directories(plugin_Spray
                             PUBLIC ${VISUALIZER_KERNEL_INCLUDES})
  
  
  target_link_libraries(plugin_Spray
    PRIVATE
    nova::dynamics
    nova::spgrid
    nova::tools
    ${Boost_LIBRARIES}
    ${assimp_LIBRARIES}
    plugin_Assimp_miniball
    )
  
  target_compile_features(plugin_Spray
    PRIVATE
    cxx_auto_type
    )
endif()
