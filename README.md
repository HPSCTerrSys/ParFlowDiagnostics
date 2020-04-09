# ana_parflow-diagnostics_PythonHeAT

Parallelizing parflow diagnostics by applying HeAT.

****I/O Interface:****

**Loading:**
Use Dictionaries to load netCDF Data into HeAT tensors. Using "File: [Variables]" loads the variables of the according files. The Method returns the HeAT Tensors in a Dictionary of the same structure: "File: Variable: Tensor".

**Saving:**
Pass the Saving Method a Dictionary of the structure "<Variable>: Tensor" and it will store the Tensor in the specified outFile. A File construct consisting of the time, levels, lat, lon dimensions will be generated. 'time' is generated as an unlimited dimension. Also, the saving Method allows slices, such that the tensor will be stored in only a part of the file, or appends the file. These slices also allow negative slicing i.e. slicing from the end. If an unlimited dimension is negatively sliced, 'the end' is defined as the currently stored data. If there are less slices than dimensions provided, these slices are used for the first dimensions (Same as numpy).

**Data structure:**
The central Data structure is the "Diagnostic" object. It loads and stores the parFlow Output and provides the methods to calculate the Diagnostic Variables. Loading of the data happens lazily, i.e. it gets loaded when it's needed (not yet implemented).

The I/O methods can be accessed without creating a "Diagnostic" object. However, "Diagnostic" objects provide wrapper functions for I/O that automatically pass the correct arguments.
