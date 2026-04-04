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



def export_lattice_to_txt(target_shape, cell: atc.AtomicCell, output_path: str):
    """
    Takes unit cell data, tiles it across the target shape's bounding box,
    filters out points outside the shape, and exports to a text file.
    """
    base_atoms = cell.atoms
    step_x = cell.step_x
    step_y = cell.step_y
    step_z = cell.step_z

    bbox = target_shape.BoundBox

    # Using a set automatically prevents duplicate coordinates from being recorded
    # if cell boundaries perfectly overlap during tiling.
    valid_points = set()

    print(f"Scanning Bounding Box: X[{bbox.XMin:.2f}, {bbox.XMax:.2f}], "
          f"Y[{bbox.YMin:.2f}, {bbox.YMax:.2f}], Z[{bbox.ZMin:.2f}, {bbox.ZMax:.2f}]")

    x = bbox.XMin - (bbox.XMin % step_x)
    while x <= bbox.XMax + step_x:
        y = bbox.YMin - (bbox.YMin % step_y)
        while y <= bbox.YMax + step_y:
            z = bbox.ZMin - (bbox.ZMin % step_z)
            while z <= bbox.ZMax + step_z:

                # Iterate through all atoms defined in our unit cell
                for (ax, ay, az) in base_atoms:
                    # Translate the base atom to the current grid location
                    pt = FreeCAD.Vector(x + ax, y + ay, z + az)

                    # Check if the translated point is inside the target shape
                    if target_shape.isInside(pt, 1e-6, True):
                        # Round to 9 decimal places to ensure clean deduplication in the set
                        valid_points.add((round(pt.x, 9), round(pt.y, 9), round(pt.z, 9)))

                z += step_z
            y += step_y
        x += step_x


    atom_count = len(valid_points)
    print(f"Found {atom_count} atoms inside the shape. Writing to file...")

    # Write the deduplicated points to the text file
    with open(output_path, 'w') as file:
        file.write(f"{atom_count}\n\n")
        for point in valid_points:
            file.write(f"H {point[0]:.9f}   {point[1]:.9f}   {point[2]:.9f}\n")

    print(f"Successfully exported coordinates to: {output_path}")

if __name__ == '__main__':
    target_shape = cm.import_part_shape(r"E:\Programming\VUT\sample-assembler\test.step")

    # bcc_cell = atc.AtomicCell()
    # bcc_cell.generate_BCC(side_length = 3)

    bcc_cell = atc.AtomicCell()
    bcc_cell.generate_BCC(side_length=3)

    try:
        export_lattice_to_txt(target_shape, bcc_cell, "bcc_test.xyz")

    except Exception as e:
        print(f"An error occurred during lattice generation: {e}")

