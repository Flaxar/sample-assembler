import FreeCAD_manager as cm
import Atomic_cells as atc
import math
try:
    import FreeCAD
    import Part
    import Import
    print("IMPORT: FreeCAD imported successfully!")
except ImportError as e:
    print(f"IMPORT: Failed to import FreeCAD: {e}")


def fill_volume_with_lattice(target_shape, atomic_cell, step_x, step_y, step_z):
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



if __name__ == '__main__':
    target_shape = cm.import_part_shape(r"C:\Programming\VUT\Sample assembler\test_object.step")

    radius = 5.0
    height = 10.0
    polygon_points = []
    for i in range(7):  # 7 points to close the loop back at the start
        angle = math.radians(60 * i)
        polygon_points.append(FreeCAD.Vector(radius * math.cos(angle), radius * math.sin(angle), 0))

    hex_wire = Part.makePolygon(polygon_points)
    hex_face = Part.Face(hex_wire)
    atomic_cell = hex_face.extrude(FreeCAD.Vector(0, 0, height))

    # 3. Define the step sizes for the grid
    # For hexagons to touch point-to-point in a simple grid (not staggered):
    dx = radius * 2.0
    dy = radius * math.sqrt(3)
    dz = height

    try:
        filled_result = fill_volume_with_lattice(target_shape, atomic_cell, dx, dy, dz)
        print(f"Successfully created filled shape with volume: {filled_result.Volume:.2f}")

        # 5. Export the final result directly to a new STEP file
        output_path = r"C:\path\to\your\FilledModel.step"
        filled_result.exportStep(output_path)
        print(f"Saved filled model to: {output_path}")

    except Exception as e:
        print(f"An error occurred during lattice generation: {e}")

