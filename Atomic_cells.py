import numpy as np
from enum import Enum
import math
try:
    import FreeCAD
    import Part
except ImportError as e:
    print(f"IMPORT: Exception {e}")

# Define all the atoms you want to use in your atomic cells.
class Atoms:
    hydrogen = (1, "H")
    helium = (2, "He")
    carbon = (6, "C")

# Parent class
class AtomicCell:
    def __init__(self):
        self.atom_elem = []
        self.atom_pos = []
        self.step_x = 0
        self.step_y = 0
        self.step_z = 0


# There are 3 basic atomic cells defined. However, you can define your own, using the parent class.
class BCC_Cell(AtomicCell):
    def __init__(self, side_length):
        super().__init__()
        self.atom_elem = []
        self.atom_pos = [
            (0, 0, 0),
            (side_length / 2, side_length / 2, side_length / 2),
        ]

        self.step_x = side_length
        self.step_y = side_length
        self.step_z = side_length
        return

    def set_atoms(self, corner,  central):
        self.atom_elem = [
            corner,
            central
        ]


class FCC_Cell(AtomicCell):
    def __init__(self, side_length):
        super().__init__()
        self.atom_pos = [
            (0, 0, 0),
            (0, side_length / 2, side_length / 2),
            (side_length / 2, 0, side_length / 2),
            (side_length / 2, side_length / 2, 0),
        ]

        self.step_x = side_length
        self.step_y = side_length
        self.step_z = side_length
        return

    def set_atoms(self, corner, sides):
        self.atom_elem = [
            corner,
            sides,
            sides,
            sides
        ]


class HCP_Cell(AtomicCell):
    def __init__(self, side_length, c_length=None):
        super().__init__()
        # 'side_length' is the basal lattice parameter 'a'
        a = side_length

        # If no 'c' height is provided, calculate the ideal HCP c/a ratio (~1.633)
        c = c_length if c_length is not None else a * math.sqrt(8 / 3)

        # The orthogonal unit cell requires 4 atoms to tile properly in 3D cartesian space.
        self.atom_pos = [
            (0, 0, 0),  # Layer A
            (a / 2, a * math.sqrt(3) / 2, 0),  # Layer A
            (a / 2, a * math.sqrt(3) / 6, c / 2),  # Layer B (nested in hollows of A)
            (0, a * math.sqrt(3) * 2 / 3, c / 2)  # Layer B
        ]

        # The steps required to seamlessly tile the 4-atom block
        self.step_x = a
        self.step_y = a * math.sqrt(3)
        self.step_z = c
        return

    def set_atoms(self, layer_a_elem, layer_b_elem):
        """
        Assign elements to the atoms.
        """
        self.atom_elem = [
            layer_a_elem,
            layer_a_elem,
            layer_b_elem,
            layer_b_elem
        ]
