# SolPosSharp

This library is a port of the NREL Solar Position and Intensity C code. Details of that code including
documentation can be found here: https://www.nrel.gov/grid/solar-resource/solpos.html.

This port preserves that code significantly except variable names are changed to be more idomatic to C#.

There is also an overload for the main `CalcSolarPosition` function that simply takes time and lat/lon
coordinates, to make it somewhat easier to call.
