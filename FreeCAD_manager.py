import sys

## ---------------- FreeCAD IMPORT --------------------

# Define your FreeCAD paths (using 'r' before the string handles the Windows backslashes)
# NOTEBOOK
# freecad_bin_path = r"C:\Program Files\FreeCAD 1.1\bin"
# freecad_lib_path = r"C:\Program Files\FreeCAD 1.1\bin\Lib"

# PC
freecad_bin_path = r"C:\Program Files\FreeCAD 1.1\bin"
freecad_lib_path = r"C:\Program Files\FreeCAD 1.1\bin\Lib"

# Append the paths to Python's system path if they aren't already there
if freecad_bin_path not in sys.path:
    sys.path.append(freecad_bin_path)

if freecad_lib_path not in sys.path:
    sys.path.append(freecad_lib_path)

# Now you can import FreeCAD and its components
try:
    import FreeCAD
    import Part
    import Import

    print("IMPORT: FreeCAD imported successfully!")
except ImportError as e:
    print(f"IMPORT: Failed to import FreeCAD: {e}")

## ----------------------------------------------------

def import_part_shape(file_path):
    try:
        new_shape = Part.Shape()
        new_shape.read(file_path)
        print(f"IMPORT: New Part.Shape sucessfully imported: {new_shape}")
    except Exception as e:
        print(f"IMPORT: Failed to import Part.Shape: {e}")
        return
    return new_shape

def import_multi_body_step(filepath):
    """Imports a STEP file and returns a list of individual solid shapes."""
    try:
        doc = FreeCAD.newDocument("TempDoc")
        Import.insert(filepath, doc.Name)
        shapes_dict = {}
        # Loop through the imported objects and grab their Label (Name) and Shape
        for obj in doc.Objects:
            # Fusion sometimes imports compound groups, we only want the actual solid shapes
            if hasattr(obj, "Shape") and not obj.Shape.isNull() and obj.Shape.Volume > 0:
                if "(Unsaved)" not in obj.Label:
                    shapes_dict[obj.Label] = obj.Shape
                    print(f"Found part named: {obj.Label}")
        return shapes_dict
    except Exception as e:
        print(f"IMPORT: Failed to import STEP file: {e}")
        return None
    # If the STEP file has multiple parts, FreeCAD imports it as a Compound.
    # We can extract the individual solids using .Solids






