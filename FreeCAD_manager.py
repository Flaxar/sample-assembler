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





