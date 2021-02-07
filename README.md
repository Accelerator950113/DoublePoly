# DoublePoly
code for **DoublePoly: A Network Visualization Method Based on Planar Polyline Drawing**  
## Environment
Julia 1.5.3  
Metapost 2.00 (Tex Live 2018)
## Get Started
```
cd code
julia -O3 Main.jl $networkName $drawingParameters
```
networkName is a string *e.g. Karate*  
drawingParameters = "n0,m0,NodeSelectMode,ConnetivityPreservationMode"  
n0,m0 are numerics *e.g. 20*  
NodeSelectMode = \[ n | ne \]  
ConnetivityPreservationMode = \[ KeepOne | KeepAll \]  
drawingParameters is optional  
### Example
```
cd code
julia -O3 Main.jl Karate 30,30,n,KeepAll
```
