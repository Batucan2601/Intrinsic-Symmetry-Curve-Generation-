#pragma once
#include <imgui/implot.h>
#include <imgui/implot_internal.h>
void imgui_new_frame()
{
    ImGui_ImplOpenGL3_NewFrame();
    ImGui_ImplGlfw_NewFrame();
    ImGui::NewFrame();
}
void imgui_setup()
{
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImPlot::CreateContext();
    ImGuiIO& io = ImGui::GetIO(); (void)io;
    ImGui::StyleColorsDark();
    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init((char*)glGetString(GL_NUM_SHADING_LANGUAGE_VERSIONS));
}
void imgui_render()
{
    ImGui::Render();
    ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
}
void imgui_close()
{
    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();
    ImPlot::DestroyContext();
}
void imgui_histogram( std::vector<float> &histogram  , int partition_size )
{
    ImGui::Begin("My Window");
    if (ImPlot::BeginPlot("My Plot")) {
        ImPlot::PlotBars("My Bar Plot", &histogram[0], partition_size);
        ImPlot::EndPlot();
    }
    ImGui::End();
}








