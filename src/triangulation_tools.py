import math
from circle_tools import dist_point_to_point
from disk import Disk
from edge import Edge
from circle_tools import tangent_point_dd, tangent_point_ld, tangent_point_between_two_elements


def find_perpendicular_line_l1_p1(l1, p):
    """
    Function to find the perpendicular line to l1 passing through point p.

    Args:
        l1 (tuple): The line
        p (tuple): The point

    Returns:
        tuple: The perpendicular line (y=ax+b)
    """

    # Unpack
    (x1, y1), (x2, y2) = l1
    x, y = p

    # Calculate slope and intercept for L1
    if (math.isclose(x2, x1, abs_tol=1e-3)):
        a_l1 = math.inf
    else:
        a_l1 = (y2 - y1) / (x2 - x1)

    # Slope of line perpendicular to L1
    a = -1 / a_l1

    # Intercept of the perpendicular line passing through point p
    b = y - a * x

    return a, b


def find_perpendicular_point_l1_p1_to_l2(l1, p, l2):
    """
    Function to find the intersection point between the perpendicular line to l1 passing through p and l2.

    Args:
        l1 (tuple): The first line
        p (tuple): The point
        l2 (tuple): The second line

    Returns:
        tuple: The intersection point
    """

    # Unpack
    (x1, y1), (x2, y2) = l1
    (x3, y3), (x4, y4) = l2
    x, y = p

    # Calculate slope and intercept for L1
    if (math.isclose(x2, x1, abs_tol=1e-9)):
        a_l1 = math.inf
    else:
        a_l1 = (y2 - y1) / (x2 - x1)

    # Calculate slope and intercept for L2
    if (math.isclose(x4, x3, abs_tol=1e-9)):
        a_l2 = math.inf
    else:
        a_l2 = (y4 - y3) / (x4 - x3)
    b2 = y3 - a_l2 * x3

    # Slope of perpendicular to L1
    a = -1 / a_l1

    # Intercept of the perpendicular line passing through point p
    b = y - a * x

    # Find intersection between the perpendicular line and L2
    if a_l2 == math.inf:
        x_intersect = x3
        y_intersect = a * x_intersect + b
    else:
        x_intersect = (b2 - b) / (a - a_l2)
        y_intersect = a_l2 * x_intersect + b2

    return x_intersect, y_intersect


def find_inscribed_circle_to_Rp(boundary):
    """
    Function to find the inscribed circle to the sub polygon defined by the boundary.

    Args:
        boundary (list): The list of disks and edges defining the boundary of a sub polygon

    Returns:
        Disk: The inscribed circle
    """
    # Find the tangent points
    (x1, y1), (x2, y2), (x3, y3) = find_tangent_points_of_Rp(boundary)[0:3]

    A1 = 2*y2 - 2*y1
    B1 = 2*x2 - 2*x1
    C1 = x1**2 - x2**2 + y1**2 - y2**2

    A2 = 2*y3 - 2*y1
    B2 = 2*x3 - 2*x1
    C2 = x1**2 - x3**2 + y1**2 - y3**2

    # A1 != A2 and B1 != B2 at the same time (otherwise the disks are all tangent at the *same* point)
    if (math.isclose(A1, A2, abs_tol=1e-9) and math.isclose(B1, B2, abs_tol=1e-9)):
        raise ValueError(
            "Two or more tangent points are the same - something went wrong")

    if math.isclose(A1*B2, A2*B1, abs_tol=1e-9):
        raise ValueError("The disks are collinear - something went wrong")
    else:
        x = (A2*C1 - A1*C2) / (A1*B2 - A2*B1)
        y = (B1*C2 - B2*C1) / (A1*B2 - A2*B1)

        r = dist_point_to_point((x, y), (x1, y1))

        return Disk(center=(x, y), radius=r)


def find_tangent_points_of_Rp(boundary):
    """
    Function to find the tangent points of the sub polygon defined by the boundary.

    Args:
        boundary (list): The list of disks and edges defining the boundary of a sub polygon

    Returns:
        3-4 tuples: The tangent points
    """
    if len(boundary) == 3:
        # if not (isinstance(boundary[0], Disk) and isinstance(boundary[1], Disk) and isinstance(boundary[2], Disk)):
        #     raise ValueError("The boundary should be a list of disks when the boundary length is 3, otherwise it is handled differently")
        # The 3 tangent points p1, p2, p3
        (x1, y1) = tangent_point_between_two_elements(boundary[0], boundary[1])
        (x2, y2) = tangent_point_between_two_elements(boundary[1], boundary[2])
        (x3, y3) = tangent_point_between_two_elements(boundary[2], boundary[0])
        return (x1, y1), (x2, y2), (x3, y3)
    elif len(boundary) == 4:
        # The 4 tangent points p1, p2, p3, p4
        if isinstance(boundary[0], Disk) and isinstance(boundary[1], Disk) and isinstance(boundary[2], Disk) and isinstance(boundary[3], Disk):
            (x1, y1) = tangent_point_dd(boundary[0], boundary[1])
            (x2, y2) = tangent_point_dd(boundary[1], boundary[2])
            (x3, y3) = tangent_point_dd(boundary[2], boundary[3])
            (x4, y4) = tangent_point_dd(boundary[3], boundary[0])
        elif isinstance(boundary[0], Disk) and isinstance(boundary[1], Disk) and isinstance(boundary[2], Disk) and isinstance(boundary[3], Edge):
            (x1, y1) = tangent_point_dd(boundary[0], boundary[1])
            (x2, y2) = tangent_point_dd(boundary[1], boundary[2])
            (x3, y3) = tangent_point_ld(boundary[3], boundary[2])
            (x4, y4) = tangent_point_ld(boundary[3], boundary[0])
        elif isinstance(boundary[0], Disk) and isinstance(boundary[1], Disk) and isinstance(boundary[2], Edge) and isinstance(boundary[3], Disk):
            (x1, y1) = tangent_point_dd(boundary[0], boundary[1])
            (x2, y2) = tangent_point_ld(boundary[2], boundary[1])
            (x3, y3) = tangent_point_ld(boundary[2], boundary[3])
            (x4, y4) = tangent_point_dd(boundary[3], boundary[0])
        elif isinstance(boundary[0], Disk) and isinstance(boundary[1], Edge) and isinstance(boundary[2], Disk) and isinstance(boundary[3], Disk):
            (x1, y1) = tangent_point_ld(boundary[1], boundary[0])
            (x2, y2) = tangent_point_ld(boundary[1], boundary[2])
            (x3, y3) = tangent_point_dd(boundary[2], boundary[3])
            (x4, y4) = tangent_point_dd(boundary[3], boundary[0])
        elif isinstance(boundary[0], Edge) and isinstance(boundary[1], Disk) and isinstance(boundary[2], Disk) and isinstance(boundary[3], Disk):
            (x1, y1) = tangent_point_ld(boundary[0], boundary[1])
            (x2, y2) = tangent_point_dd(boundary[1], boundary[2])
            (x3, y3) = tangent_point_dd(boundary[2], boundary[3])
            (x4, y4) = tangent_point_ld(boundary[0], boundary[3])
        elif isinstance(boundary[0], Edge) and isinstance(boundary[1], Disk) and isinstance(boundary[2], Edge) and isinstance(boundary[3], Disk):
            (x1, y1) = tangent_point_ld(boundary[0], boundary[1])
            (x2, y2) = tangent_point_ld(boundary[2], boundary[1])
            (x3, y3) = tangent_point_ld(boundary[2], boundary[3])
            (x4, y4) = tangent_point_ld(boundary[0], boundary[3])
        elif isinstance(boundary[0], Disk) and isinstance(boundary[1], Edge) and isinstance(boundary[2], Disk) and isinstance(boundary[3], Edge):
            (x1, y1) = tangent_point_ld(boundary[1], boundary[0])
            (x2, y2) = tangent_point_ld(boundary[1], boundary[2])
            (x3, y3) = tangent_point_ld(boundary[3], boundary[2])
            (x4, y4) = tangent_point_ld(boundary[3], boundary[0])
        else:
            raise ValueError(
                "3 or more edges in a 4 remainder region - something went wrong")
        return (x1, y1), (x2, y2), (x3, y3), (x4, y4)

    else:
        raise ValueError("The number of disks should be 3 or 4")
