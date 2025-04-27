# Intrinsic Symmetry Curve Generation in 3D Meshes

Official implementation of intrinsic symmetry curve generation in 3D meshes.

---

## Folder Structure

- **Executable/**  
  - `GraphicsProgram.exe`  
    Launches the graphics interface used during development.  
  - `Console.exe`  
    Command-line tool for generating outputs.

---

## How to Use

### Graphics Interface

1. **Launch**  
   Run `GraphicsProgram.exe`.  
2. **Load Dataset**  
   Mesh → Dataset → choose one of SCAPE, TOSCA or Princeton.  
3. **Generate Descriptors**  
   Mesh → NLateral → Generate Descriptors (with midpoint).  
4. **Find Bisector Curve**  
   Mesh → NLateral → Voronoi test → Get closest curve with matches.  
5. **Prune Curve**  
   Mesh → NLateral → Prune curve.

### Console Program

```bash
Console.exe \
  -I inputmesh.off \
  -OCors correspondence.txt \
  -OAxis bisector.txt \
  -SampleNo 5

## Input Data

Sample input data can be found under the `Trilateral/Mesh/SCB/Data/` directory.  

Place your `.off` mesh files and any accompanying data in this folder before running either the graphics or console programs.  

