# Sample assembler
This program takes a .STEP file with one or multiple components and then fills the shapes accordingly with atoms.
The result is saved as a .xyz file.
## Prerequisites
- [ ] Download [FreeCAD](www.freecad.org/downloads.php?lang=en)
- [ ] Install python version according to the FreeCAD version

## Setup
1. Clone this repository
2. In `FreeCAD_manager.py`, change the absolute path to your FreeCAD folder

## How to use
1. In your CAD program of choice, create one or multiple components and rename them as you wish.
2. Export the components as a single .STEP file.
3. Open your main Python script and update the file path to point to your exported .STEP file. 
4. Retrieve your individual components by their CAD names using the import function (e.g., `part_A = my_parts["ComponentA"]`). 
5. Initialize your desired atomic cell structures (e.g., `FCC_Cell, HCP_Cell`) with your chosen lattice parameters. 
6. **Crucial:** Assign specific elements to each cell using the `.set_atoms()` method to ensure the exporter has element data to write. 
7. Map your CAD components to their corresponding atomic cells by grouping them into a list of tuples. 
8. Run the script. The program will automatically filter the lattice grids to your CAD boundaries and generate your final .xyz file.