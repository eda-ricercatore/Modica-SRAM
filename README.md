Modica-SRAM
===========

Design of a 32-kbit synchronous SRAM with 32-bit words, using 180 nm process technology.  

Developed *MATLAB* scripts to evaluate architectural trade-offs between performance (using logical effort analysis) and area usage; see the source code for the *HSPICE* decks and *MATLAB* scripts that are used during architectural trade-off evaluation, and characterization of inverters for different supply voltages (VDD) and temperatures.

Essentially, the aforementioned *MATLAB* scripts serve as a primitive memory compiler for the design space exploration of the SRAM's floorplan.

It also includes *HSPICE* decks for the characterization of the 6-transistor SRAM cell for different transistor ratios, and the SRAM read and write circuitry.  

Co-designed and co-developed the SRAM using schematic entry in *Cadence Virtuoso*.  

Performed functional and timing verification by simulating extracted *SPICE* netlist in *NanoSim*.


See [sram.pdf](https://github.com/eda-ricercatore/Modica-SRAM/blob/master/report/sram.pdf) for a report of this SRAM project.
