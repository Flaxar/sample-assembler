import numpy as np
try:
    import FreeCAD
    import Part
except ImportError as e:
    print(f"IMPORT: Exception {e}")

class Atomic_cell:
    def __init__(self, type):
        pass



    def fill_volume_with_lattice(self, target_shape, atomic_cell, step_x, step_y, step_z):
        """
        Fills a target shape with a repeated atomic cell.
        """
        bbox = target_shape.BoundBox
        cells = []

        print(f"Target Bounding Box: X[{bbox.XMin:.2f}, {bbox.XMax:.2f}], "
              f"Y[{bbox.YMin:.2f}, {bbox.YMax:.2f}], Z[{bbox.ZMin:.2f}, {bbox.ZMax:.2f}]")

        # Generate the grid. We use a small buffer to ensure the very edges aren't missed.
        x = bbox.XMin
        while x <= (bbox.XMax + step_x):
            y = bbox.YMin
            while y <= (bbox.YMax + step_y):
                z = bbox.ZMin
                while z <= (bbox.ZMax + step_z):
                    # Copy the base cell and move it to the current coordinate
                    cell_copy = atomic_cell.copy()
                    cell_copy.translate(FreeCAD.Vector(x, y, z))
                    cells.append(cell_copy)
                    z += step_z
                y += step_y
            x += step_x

        print(f"Generated {len(cells)} atomic cells. Performing boolean intersection...")

        # Combine all cells into one compound, then intersect with the target shape
        lattice_compound = Part.makeCompound(cells)
        filled_shape = target_shape.common(lattice_compound)

        return filled_shape