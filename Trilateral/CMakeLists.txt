set(PROJECT_NAME Trilateral)

################################################################################
# Source groups
################################################################################
set(Header_Files
    "src/Application/Include/Callback.h"
    "src/Application/Include/imgui_stuff.h"
    "src/Application/Include/Mesh.h"
    "src/Application/Include/Mesh_imgui.h"
    "src/Application/Include/MeshFactory.h"
    "src/Application/Include/MeshFunctions.h"
    "src/Application/Include/Prototypes.h"
    "src/Application/Include/RandomRocket.h"
    "src/Application/Include/Shader.h"
    "src/Application/Include/Simulation.h"
    "src/Application/Include/TrilateralMap.h"
    "src/Application/Include/TrueInformation.h"
    "src/Application/Include/Sampling.h"
    "src/Application/Include/Laplace-Beltrami.h"
    "src/Application/Include/DominantSymmetry.h"
    "src/Application/Include/CoreTypeDefs.h"


)
source_group("Header Files" FILES ${Header_Files})

set(Source_Files
    "src/Application/main.cpp"
    "src/Application/Source/Mesh.cpp"
    "src/Application/Source/MeshFactory.cpp"
    "src/Application/Source/Vehicle.cpp"
    "src/Application/Source/Sampling.cpp"
    "src/Application/Source/Laplace-Beltrami.cpp"
    "src/Application/Source/DominantSymmetry.cpp"
    "src/Application/Source/CoreTypeDefs.cpp"
    "src/Application/Source/TrilateralMap.cpp"
    "src/External/Include/imgui/imgui.cpp"
    "src/External/Include/imgui/imgui_demo.cpp"
    "src/External/Include/imgui/imgui_draw.cpp"
    "src/External/Include/imgui/imgui_impl_glfw.cpp"
    "src/External/Include/imgui/imgui_impl_opengl3.cpp"
    "src/External/Include/imgui/imgui_stdlib.cpp"
    "src/External/Include/imgui/imgui_tables.cpp"
    "src/External/Include/imgui/imgui_widgets.cpp"
    "src/External/Include/imgui/implot.cpp"
    "src/External/Include/imgui/implot_demo.cpp"
    "src/External/Include/imgui/implot_items.cpp"

)
source_group("Source Files" FILES ${Source_Files})

set(ALL_FILES
    ${Header_Files}
    ${Source_Files}
)

################################################################################
# Target
################################################################################
add_executable(${PROJECT_NAME} ${ALL_FILES})

use_props(${PROJECT_NAME} "${CMAKE_CONFIGURATION_TYPES}" "${DEFAULT_CXX_PROPS}")
set(ROOT_NAMESPACE Trilateral)

set_target_properties(${PROJECT_NAME} PROPERTIES
    VS_GLOBAL_KEYWORD "Win32Proj"
)
if("${CMAKE_VS_PLATFORM_NAME}" STREQUAL "x64")
    set_target_properties(${PROJECT_NAME} PROPERTIES
        INTERPROCEDURAL_OPTIMIZATION_RELEASE "TRUE"
    )
elseif("${CMAKE_VS_PLATFORM_NAME}" STREQUAL "x86")
    set_target_properties(${PROJECT_NAME} PROPERTIES
        INTERPROCEDURAL_OPTIMIZATION_RELEASE "TRUE"
    )
endif()
################################################################################
# Include directories
################################################################################
if("${CMAKE_VS_PLATFORM_NAME}" STREQUAL "x64")
    target_include_directories(${PROJECT_NAME} PUBLIC
        "$<$<CONFIG:Debug>:"
            "$ENV{SolutionDIr}/Trilateral/src/Application/Include;"
            "${CMAKE_CURRENT_SOURCE_DIR}/src/External/Include"
        ">"
    )
endif()

################################################################################
# Compile definitions
################################################################################
if("${CMAKE_VS_PLATFORM_NAME}" STREQUAL "x64")
    target_compile_definitions(${PROJECT_NAME} PRIVATE
        "$<$<CONFIG:Debug>:"
            "_DEBUG;"
            "_CRT_SECURE_NO_WARNINGS"
        ">"
        "$<$<CONFIG:Release>:"
            "NDEBUG"
        ">"
        "_CONSOLE;"
        "UNICODE;"
        "_UNICODE"
    )
elseif("${CMAKE_VS_PLATFORM_NAME}" STREQUAL "x86")
    target_compile_definitions(${PROJECT_NAME} PRIVATE
        "$<$<CONFIG:Debug>:"
            "_DEBUG"
        ">"
        "$<$<CONFIG:Release>:"
            "NDEBUG"
        ">"
        "WIN32;"
        "_CONSOLE;"
        "UNICODE;"
        "_UNICODE"
    )
endif()

################################################################################
# Compile and link options
################################################################################
if(MSVC)
    if("${CMAKE_VS_PLATFORM_NAME}" STREQUAL "x64")
        target_compile_options(${PROJECT_NAME} PRIVATE
            $<$<CONFIG:Release>:
                /Oi;
                /Gy
            >
            /permissive-;
            /sdl;
            /W3;
            ${DEFAULT_CXX_DEBUG_INFORMATION_FORMAT};
            ${DEFAULT_CXX_EXCEPTION_HANDLING}
        )
    elseif("${CMAKE_VS_PLATFORM_NAME}" STREQUAL "x86")
        target_compile_options(${PROJECT_NAME} PRIVATE
            $<$<CONFIG:Release>:
                /Oi;
                /Gy
            >
            /permissive-;
            /sdl;
            /W3;
            ${DEFAULT_CXX_DEBUG_INFORMATION_FORMAT};
            ${DEFAULT_CXX_EXCEPTION_HANDLING}
        )
    endif()
    if("${CMAKE_VS_PLATFORM_NAME}" STREQUAL "x64")
        target_link_options(${PROJECT_NAME} PRIVATE
            $<$<CONFIG:Release>:
                /OPT:REF;
                /OPT:ICF
            >
            /DEBUG;
            /SUBSYSTEM:CONSOLE
        )
    elseif("${CMAKE_VS_PLATFORM_NAME}" STREQUAL "x86")
        target_link_options(${PROJECT_NAME} PRIVATE
            $<$<CONFIG:Release>:
                /OPT:REF;
                /OPT:ICF
            >
            /DEBUG;
            /SUBSYSTEM:CONSOLE
        )
    endif()
endif()

################################################################################
# Dependencies
################################################################################
if("${CMAKE_VS_PLATFORM_NAME}" STREQUAL "x64")
    set(ADDITIONAL_LIBRARY_DEPENDENCIES
        "$<$<CONFIG:Debug>:"
            "glew32;"
            "glfw3;"
            "opengl32"
        ">"
    )
endif()
target_link_libraries(${PROJECT_NAME} PRIVATE "${ADDITIONAL_LIBRARY_DEPENDENCIES}")

if("${CMAKE_VS_PLATFORM_NAME}" STREQUAL "x64")
    target_link_directories(${PROJECT_NAME} PRIVATE
        "$<$<CONFIG:Debug>:"
            "${CMAKE_SOURCE_DIR}//Trilateral/Dependencies"
        ">"
    )
endif()

