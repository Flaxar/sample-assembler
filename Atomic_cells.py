import numpy as np
from enum import Enum
import math
try:
    import FreeCAD
    import Part
except ImportError as e:
    print(f"IMPORT: Exception {e}")

class Atoms:
    hydrogen = (1, "H")
    helium = (2, "He")


class AtomicCell:
    def __init__(self):
        # self.atomic_cell_type = atomic_cell_type
        self.atom_elem = []
        self.atom_pos = []
        self.step_x = 0
        self.step_y = 0
        self.step_z = 0



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






# def generate_square_atomic_cell(side_len):
#     return Part.makeBox(side_len, side_len, side_len)
#
#
#
# def generate_hexagonal_atomic_cell(radius, height):
#     polygon_points = []
#     for i in range(7):  # 7 points to close the loop back at the start
#         angle = math.radians(60 * i)
#         polygon_points.append(FreeCAD.Vector(radius * math.cos(angle), radius * math.sin(angle), 0))
#
#     hex_wire = Part.makePolygon(polygon_points)
#     hex_face = Part.Face(hex_wire)
#     return hex_face.extrude(FreeCAD.Vector(0, 0, height))
