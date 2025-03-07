from disk_packing.solver import circleTangentToCCL
from disk_packing.circle_tools import tangency_correct_part_of_element, correct_side_edge


def disk_clc(elm1_tup, elm2_tup, elm3_tup):
    """
    Function that, given a disk, an edge and a third element, a disk, return the possible disks that are tangent to all three elements.

    Args:
        elm1_tup (tuple): Tuple of the element (disk) and the corresponding disk
        elm2_tup (tuple): Tuple of the element (disk) and the corresponding disk
        elm3_tup (tuple): Tuple of the element (edge) and the corresponding disk
        mid_elm (Disk): Disk that is the middle element between the two first elements - if given

    Returns:
        List of solution disks
    """

    d1 = elm1_tup[1]
    e = elm2_tup[1]
    d2 = elm3_tup[1]

    solutions = circleTangentToCCL(d1, d2, e, -1, -1)

    # Filter the solutions that are tangent to the correct part
    solutions = tangency_correct_part_of_element(solutions, elm1_tup)
    solutions = tangency_correct_part_of_element(solutions, elm2_tup)
    solutions = tangency_correct_part_of_element(solutions, elm3_tup)

    solutions = correct_side_edge(solutions, e)

    return solutions
