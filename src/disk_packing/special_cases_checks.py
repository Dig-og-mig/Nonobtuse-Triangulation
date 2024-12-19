from disk_packing.circle_tools import tangent_point_dd, tangent_point_ld, angle_between_points, cross_product
from triangulation.triangulation_tools import find_tangent_points_of_Rp
from models.disk import Disk
import math


def arcs_less_than_180(boundary: list[Disk]):
    """
    Function to check if the angles between the tangent points are less than 180 degrees (if the element is a disk).

    Args:
        boundary (list): The list of disks and edges defining the boundary of a sub polygon

    Returns:
        bool: Whether the angles are less than 180 degrees
    """

    for i in range(len(boundary)):
        if isinstance(boundary[i], Disk):
            # Find tangent points
            if isinstance(boundary[(i - 1) % len(boundary)], Disk):
                x1, y1 = tangent_point_dd(
                    boundary[(i - 1) % len(boundary)], boundary[i])
            else:
                x1, y1 = tangent_point_ld(
                    boundary[(i - 1) % len(boundary)], boundary[i])

            if isinstance(boundary[(i + 1) % len(boundary)], Disk):
                x2, y2 = tangent_point_dd(
                    boundary[i], boundary[(i + 1) % len(boundary)])
            else:
                x2, y2 = tangent_point_ld(
                    boundary[(i + 1) % len(boundary)], boundary[i])

            # Check if the angle is less than 180 degrees
            angle = angle_between_points(
                boundary[i].center, (x1, y1), (x2, y2), True)
            if angle > 180 and not math.isclose(angle, 180, abs_tol=1e-3):
                return False, boundary[i]
    return True, None


def center_is_inside_convex_hull(c: tuple[float, float], boundary: list[Disk]):
    """
    Function to check if the center of a disk is inside the convex hull defined by the boundary.

    Args:
        c (tuple): The center of the disk
        boundary (list): The list of disks and edges defining the boundary of a sub polygon

    Returns:
        bool: Whether the center is inside the convex hull
    """
    R_boundary = find_tangent_points_of_Rp(boundary)
    for i in range(len(boundary)):
        (x1, y1), (x2, y2) = R_boundary[i], R_boundary[(
            i + 1) % len(R_boundary)]
        if cross_product((x1, y1), (x2, y2), c) < 0:
            return False, boundary[(i+1) % len(boundary)]
    return True, None
