import sys

## ---------------- FreeCAD IMPORT --------------------

# Define your FreeCAD paths (using 'r' before the string handles the Windows backslashes)
freecad_bin_path = r"C:\Program Files\FreeCAD 1.1\bin"
freecad_lib_path = r"C:\Program Files\FreeCAD 1.1\bin\Lib"

# Append the paths to Python's system path if they aren't already there
if freecad_bin_path not in sys.path:
    sys.path.append(freecad_bin_path)

if freecad_lib_path not in sys.path:
    sys.path.append(freecad_lib_path)

try:
    import FreeCAD
    import Part
    import Import

    print("IMPORT: FreeCAD imported successfully!")
except ImportError as e:
    print(f"IMPORT: Failed to import FreeCAD: {e}")

## ----------------------------------------------------

def __import_part_shape(file_path):
    """
    Deprecated. Use import_multi_body_step() instead.
    This function returns the object as a variable directly.
    However, the main function works with a dictionary.
    """
    try:
        new_shape = Part.Shape()
        new_shape.read(file_path)
        print(f"IMPORT: New Part.Shape sucessfully imported: {new_shape}")
    except Exception as e:
        print(f"IMPORT: Failed to import Part.Shape: {e}")
        return
    return new_shape

def import_multi_body_step(filepath):
    """
    Imports a STEP file and returns a dictionary of individual solid shapes, like this:
    {
        "Shape 1 name": shape1_object,
        "Shape 2 name": shape2_object,
        ...
    }
    The shape names are defined in the CAD software before exporting into a STEP file.
    """
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






