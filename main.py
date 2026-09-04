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


def export_multiple_lattices_to_txt(part_data, output_path: str):
    """
    Takes a list of (shape, cell) pairs, tiles the respective cells across
    each shape's bounding box, and exports all valid points to a single text file.

    Using a dictionary keyed by (X, Y, Z) ensures that if two parts share
    a boundary face, an atom on that face is only recorded once.
    """

    master_points = {}

    for target_shape, cell in part_data:
        base_atoms = cell.atom_pos
        step_x = cell.step_x
        step_y = cell.step_y
        step_z = cell.step_z

        bbox = target_shape.BoundBox

        print(f"Scanning Bounding Box: X[{bbox.XMin:.2f}, {bbox.XMax:.2f}], "
              f"Y[{bbox.YMin:.2f}, {bbox.YMax:.2f}], Z[{bbox.ZMin:.2f}, {bbox.ZMax:.2f}]")

        x = bbox.XMin - (bbox.XMin % step_x)
        while x <= bbox.XMax + step_x:
            y = bbox.YMin - (bbox.YMin % step_y)
            while y <= bbox.YMax + step_y:
                z = bbox.ZMin - (bbox.ZMin % step_z)
                while z <= bbox.ZMax + step_z:

                    # Iterate through all atoms defined in our unit cell
                    for index, (ax, ay, az) in enumerate(base_atoms):
                        pt = FreeCAD.Vector(x + ax, y + ay, z + az)

                        # Check if the translated point is inside the current shape
                        if target_shape.isInside(pt, 1e-6, True):
                            rx = round(pt.x, 9)
                            ry = round(pt.y, 9)
                            rz = round(pt.z, 9)
                            # Add to dictionary. If the coordinate already exists
                            # from an adjacent shape, this overwrites/merges it.
                            element_data = cell.atom_elem[index]
                            master_points[(rx, ry, rz)] = element_data

                    z += step_z
                y += step_y
            x += step_x

    atom_count = len(master_points)
    print(f"Found {atom_count} total atoms across all shapes. Writing to file...")

    # Write the unique points to the text file
    with open(output_path, 'w') as file:
        file.write(f"{atom_count}\n\n")

        for coords, element_data in master_points.items():
            # coords is (X, Y, Z)
            # element_data maintains your original cell.atom_elem structure
            file.write(f"{element_data[1]} {coords[0]:.9f}   {coords[1]:.9f}   {coords[2]:.9f}\n")

    print(f"Successfully exported coordinates to: {output_path}")


if __name__ == '__main__':
    """
    Initialize the atomic cells you want to use. 
    You can also define your own atomic cells in Atomic_cells.py
    """
    bcc_cell = atc.BCC_Cell(6)
    bcc_cell.set_atoms(atc.Atoms.hydrogen, atc.Atoms.helium)

    hcp_cell = atc.HCP_Cell(6)
    hcp_cell.set_atoms(atc.Atoms.helium, atc.Atoms.hydrogen)

    fcc_cell = atc.FCC_Cell(6)
    fcc_cell.set_atoms(atc.Atoms.carbon, atc.Atoms.hydrogen)

    # Change the path to your .step file.
    # multi-body-test.step is included in the repository for a test run.
    bodies = cm.import_multi_body_step(r"E:\Programming\VUT\sample-assembler\multi-body-test.step")
    compA = bodies["ComponentA"] # Use the same name as you did in the CAD
    compB = bodies["ComponentB"]
    compC = bodies["ComponentC"]

    # Here, you define which shape will use which atomic cell.
    assembly_data = [
        (compA, bcc_cell),
        (compB, hcp_cell),
        (compC, fcc_cell)
    ]

    try:
        export_multiple_lattices_to_txt(assembly_data, "Assembly_Atoms.xyz")
    except Exception as e:
        print(f"An error occurred during lattice generation: {e}")

