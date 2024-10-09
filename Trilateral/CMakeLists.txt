set(PROJECT_NAME Trilateral)

################################################################################
# Source groups
################################################################################
set(Header_Files
    "src/Application/Include/Callback.h"
    "src/Application/Include/imgui_stuff.h"
    "src/Application/Include/TrilateralMesh.h"
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
    "src/Application/Include/SymmetryAwareEmbeddingForShapeCorrespondence.h"
    "src/Application/Include/Skeleton.h"
    "src/Application/Include/MetricCalculations.h"
    "src/Application/Include/NLateralDescriptor.h"
    "src/Application/Include/ShapeDiameter.h"
    "src/Application/Include/Ray.h"
    "src/Application/Include/SkeletalNLateral.h"
    "src/Application/Include/FuzzyGeodesic.h"
    "src/Application/Include/Geodesic.h"
    "src/Application/Include/DvorakEstimatingApprox.h"
    "src/Application/Include/HeatKernelSignature.h"
    "src/Application/Include/ImguiMenuBar.h"
    "src/Application/Include/RaylibExtras.h"
    "src/Application/Include/ROI.h"
    "src/Application/Include/Histogram.h"
    "src/Application/Include/SpinImage.h"
    "src/Application/Include/HistogramFunctions.h"

    
    "src/External/Include/rlImGui/imgui_impl_raylib.h"
    "src/External/Include/rlImGui/rlImgui.h"
    "src/External/Include/rlImGui/rlImGuiColors.h"
    
    "src/External/Include/imgui/ImGuiFileDialog.h"
    "src/External/Include/imgui/ImGuiFileDialogConfig.h"


    


)
source_group("Header Files" FILES ${Header_Files})

set(Source_Files
    "src/Application/main.cpp"
    "src/Application/Source/TrilateralMesh.cpp"
    "src/Application/Source/MeshFactory.cpp"
    "src/Application/Source/Sampling.cpp"
    "src/Application/Source/Laplace-Beltrami.cpp"
    "src/Application/Source/DominantSymmetry.cpp"
    "src/Application/Source/CoreTypeDefs.cpp"
    "src/Application/Source/TrilateralMap.cpp"
    "src/External/Include/imgui/imgui.cpp"
    "src/External/Include/imgui/imgui_demo.cpp"
    "src/External/Include/imgui/imgui_draw.cpp"
    "src/External/Include/imgui/imgui_impl_opengl3.cpp"
    "src/External/Include/imgui/imgui_stdlib.cpp"
    "src/External/Include/imgui/imgui_tables.cpp"
    "src/External/Include/imgui/imgui_widgets.cpp"
    "src/External/Include/imgui/implot.cpp"
    "src/External/Include/imgui/implot_demo.cpp"
    "src/External/Include/imgui/implot_items.cpp"
    "src/Application/Source/SymmetryAwareEmbeddingForShapeCorrespondence.cpp"
    "src/Application/Source/Skeleton.cpp"
    "src/Application/Source/MetricCalculations.cpp"
    "src/Application/Source/NLateralDescriptor.cpp"
    "src/Application/Source/ShapeDiameter.cpp"
    "src/Application/Source/Ray.cpp"
    "src/Application/Source/SkeletalNLateral.cpp"
    "src/Application/Source/FuzzyGeodesic.cpp"
    "src/Application/Source/Geodesic.cpp"
    "src/Application/Source/DvorakEstimatingApprox.cpp"
    "src/Application/Source/HeatKernelSignature.cpp"
    "src/Application/Source/ImguiMenuBar.cpp"
    "src/Application/Source/RaylibExtras.cpp"
    "src/Application/Source/ROI.cpp"
    "src/Application/Source/Histogram.cpp"
    "src/Application/Source/SpinImage.cpp"
    "src/Application/Source/HistogramFunctions.cpp"


    "src/External/Include/imgui/ImGuiFileDialog.cpp"

    "src/External/Include/rlImGui/rlImgui.cpp"


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
target_include_directories(${PROJECT_NAME} PUBLIC  "$ENV{SolutionDIr}/Trilateral/src/Application/Include;"
            "${CMAKE_CURRENT_SOURCE_DIR}/src/External/Include")

target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src/External/Include/raylib/src)
target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src/External/Include/imgui)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/src/External/Include/raylib)
target_link_libraries(${PROJECT_NAME} PUBLIC raylib)

target_link_libraries(${PROJECT_NAME} PRIVATE "${ADDITIONAL_LIBRARY_DEPENDENCIES}"  "glew32"
"glfw3"
"opengl32")
target_include_directories(${PROJECT_NAME} PUBLIC  ${CMAKE_CURRENT_SOURCE_DIR}/src/External/Include/spectra)
target_include_directories(${PROJECT_NAME} PUBLIC  ${CMAKE_CURRENT_SOURCE_DIR}/src/External/Include/eigen)
target_link_directories(${PROJECT_NAME} PUBLIC   "${CMAKE_SOURCE_DIR}//Trilateral/Dependencies")


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


