Official implementation of Intrinsic Symmetry Curve Generation In 3D Meshes
The executables are on Executable folder, GraphicsProgram launches the graphics interface that we mainly used for developing the program,
The Console.exe is just for generating outputs.

**How to use**
GraphicsProgram.exe should launch without problme. From there you can load prebuild datasets namely SCAPE TOSCA and Princeton using Mesh->Dataset
Then with proper parameters you can run Mesh->NLateral->Generate Descritors with midpoint.
After initial matchings you can generate best suited bisector curve with Mesh->Nlateral->Voronoi test-> Get closest curve with matches
For pruning there is also a prune curve in the same section.

For console program you can do the following
Console.exe -I inputmesh.off -OCors correspondence.txt -OAxis bisector.txt  -SampleNo 5 
-I is for inputmesh
-OCors for correspondences
-OAxis for bisector region points
-SampleNo is for how many Minimum geodesic  cycles after Average geodesic distance function. 


You can also build the code yourselves, 
cmake -B ./build -A x64 
builds for graphics support
cmake -B ./build -A x64 -DCONSOLE_MODE=ON builds for console mode

