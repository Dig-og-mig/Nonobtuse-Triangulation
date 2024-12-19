from edge import Edge
from solver import circleTangentToCCL
from circle_tools import not_intersect, tangency_correct_part_of_element, correct_side_edge


def disk_ccl(elm1_tup, elm2_tup, elm3_tup, mid_elm = None):
    """
    Function that, given two disks and an edge, return the possible disks that are tangent to to all three elements.

    Args:
        elm1_tup (tuple): Tuple of the element (disk) and the corresponding disk
        elm2_tup (tuple): Tuple of the element (disk) and the corresponding disk
        elm3_tup (tuple): Tuple of the element (edge) and the corresponding disk
        mid_elm (Disk): Disk that is the middle element between the two first elements - if given

    Returns:
        List of solution disks
    """

    c1 = elm1_tup[1]
    c2 = elm2_tup[1]
    l = elm3_tup[1]

    solutions = circleTangentToCCL(c1, c2, l, -1, -1)

    # If the middle element is an edge, we need to check if the solutions intersect the edge
    if isinstance(mid_elm, Edge):
        solutions = not_intersect(solutions, mid_elm)

    # Filter the solutions that are tangent to the correct part
    solutions = tangency_correct_part_of_element(solutions, elm1_tup)
    solutions = tangency_correct_part_of_element(solutions, elm2_tup)
    solutions = tangency_correct_part_of_element(solutions, elm3_tup)

    # Filter solutions - for edges - that are on the correct side of the edge
    solutions = correct_side_edge(solutions, l)

    return solutions
