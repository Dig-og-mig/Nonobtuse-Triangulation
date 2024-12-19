from disk_packing.solver import circleTangentToCCL
from disk_packing.circle_tools import tangency_correct_part_of_element, correct_side_edge


def disk_lcc(elm1_tup, elm2_tup, elm3_tup):
    """
    Function that, given an edge and two disks, return the possible disks that are tangent to all three elements.

    Args:
        elm1_tup (tuple): Tuple of the element (edge) and the corresponding disk
        elm2_tup (tuple): Tuple of the element (disk) and the corresponding disk
        elm3_tup (tuple): Tuple of the element (disk) and the corresponding disk
        mid_elm (Disk): Disk that is the middle element between the two first elements - if given

    Returns:
        List of solution disks
    """

    l = elm1_tup[1]
    c1 = elm2_tup[1]
    c2 = elm3_tup[1]

    solutions = circleTangentToCCL(c1, c2, l, -1, -1)

    # Filter the solutions that are tangent to the correct part
    solutions = tangency_correct_part_of_element(solutions, elm1_tup)
    solutions = tangency_correct_part_of_element(solutions, elm2_tup)
    solutions = tangency_correct_part_of_element(solutions, elm3_tup)

    solutions = correct_side_edge(solutions, l)

    return solutions
