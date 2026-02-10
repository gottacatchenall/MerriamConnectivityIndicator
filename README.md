# MerriamConnectivityIndiactor.jl

An extension of Circuitscape and Omniscape for computing landscape connectivity using a sliding circular window, where current is injected at 'spokes' on the boundary of the sliding window, injects current into these spokes, and measures the effective resistance of current flow to the center of the sliding window, which is set to ground.


## File Structure

- `Project.toml`: Julia dependency management file
- src/
    - `MerriamConnectivityIndicator.jl`: Root module for the package
    - `types.jl`: Type definitions for model configuration (`MCIConfig`) and spokes (`SpokePoint`)
    - `spokes.jl`: Method for generating spoke point Cartesian indices for a given sliding window location, radius, and number of spokes.
    - `window.jl`: Methods for getting the relevant inputs for a given sliding window (conductance, source, and ground matrices)
    - `solver.jl`: Methods for converting sliding window inputs into inputs to Circuitscape, and using Circuitscape internals to solve for effective resistance.  
    - `mci.jl`: Main entrypoint for users; methods for taking a raster and running a sliding window over it. 
- notebooks/
    - `comparing_mci_to_circuitscape.ipynb`: Tests to compare MCI implementation to it computed directly with Circuitscape functions  
    - `example.ipynb`: An example MCI analysis run on a demo landscape from Omniscape

