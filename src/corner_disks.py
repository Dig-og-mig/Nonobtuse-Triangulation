import math

from edge import Edge
from disk import Disk
from corner_disks_distance import dist_point_to_line_inside_corner
from circle_tools import angle_between_points, projection_scalar_for_line, dist_point_to_line


def disk_at_convex_corner(e1: Edge, e2: Edge, edges: list[Edge]):
    """
    Function to find a disk at a convex corner of the polygon

    Args:
        e1 (Edge): First edge at the corner
        e2 (Edge): Second edge at the corner
        - It is assumed that e1.end = e2.start
        edges (list[Edge]): List of the edges in the polygon (boundary and holes)

    Returns:
        Disk: Object representing the disk at the corner
    """
    smallestDistance = math.inf
    corner_point = e1.end  # = e2.start

    # Find smallest distance to any edge inside the corner
    for edge in edges:
        dist = dist_point_to_line_inside_corner(corner_point, edge, e1, e2)

        if dist < smallestDistance:
            smallestDistance = dist

    # Radius of the circle with center at the corner point
        # Any radius strictly smaller than d = smallestDistance/4 would work - we have chosen d/2
    radius = smallestDistance / 8

    # Find angle for solution center
    angle = angle_between_points(e1.end, e2.end, e1.start) / 2
    angle_sol = angle + \
        angle_between_points(
            corner_point, (corner_point[0] + 1, corner_point[1]), e2.end)

    # Solution center point placed with distance = radius from corner point and rotated to the middle of the two edges
    solution_center = (corner_point[0] + radius*math.cos(math.radians(
        angle_sol)), corner_point[1] + radius*math.sin(math.radians(angle_sol)))

    # Find projection to check if the solution is on the line segment
    proj_onto_line = projection_scalar_for_line(solution_center, e1)
    if proj_onto_line > 1.0 or proj_onto_line < 0.0:
        raise ValueError(
            "The projection is not on the line segment - something went wrong for a corner disk (convex case)")

    # Find the distance from the center perpendicular to one of the edges
    solution_radius = dist_point_to_line(solution_center, e1)

    if e1.hole:
        return Disk(center=solution_center, radius=solution_radius, hole=True)
    else:
        return Disk(center=solution_center, radius=solution_radius)


def disk_at_concave_corner(e1: Edge, e2: Edge, edges: list[Edge]):
    """
    Function to find two disks at a concave corner of the polygon

    Args:
        e1 (Edge): First edge at the corner
        e2 (Edge): Second edge at the corner
        - It is assumed that e1.end = e2.start
        edges (list[Edge]): List of the edges in the polygon (boundary and holes)

    Returns:
        (Disk, Disk): Tuple of Disk's representing the disk at the corner
    """

    # Corner point
    corner_point = e1.end  # = e2.start

    # Find edge in the middle of e1 and e2 (in the concave corner) (with same length as e2)
    # the (concave) angle between e1 and e2
    angle = angle_between_points(corner_point, e2.end, e1.start) / 2
    vec2_x, vec2_y = (e2.end[0] - e2.start[0],
                      e2.end[1] - e2.start[1])  # place e2 at origin
    mid_end_x = corner_point[0] + (vec2_x * math.cos(math.radians(
        angle)) - vec2_y * math.sin(math.radians(angle)))  # placed at corner_point
    mid_end_y = corner_point[1] + (vec2_x * math.sin(math.radians(
        angle)) + vec2_y * math.cos(math.radians(angle)))  # placed at corner_point

    mid_edge = Edge(start=corner_point, end=(mid_end_x, mid_end_y))

    # Find smalles distance to any edge inside the concave corner
    smallestDistance = math.inf
    for edge in edges:
        # Find the smallest distance from the corner to the given edge (inside the corner)
        dist = dist_point_to_line_inside_corner(corner_point, edge, e1, e2)

        # If the distance is smaller - update the smallest distance
        if dist < smallestDistance:
            smallestDistance = dist

    # Radius of the circle with center at the corner point
    # Any radius strictly smaller than d = smallestDistance/4 would work - we have chosen d/2
    radius = smallestDistance / 8

    #### First disk ####

    # Average angle between e1, e2:
    angle = angle_between_points(corner_point, mid_edge.end, e1.start) / 2
    # Total angle from the x-axis for center point
    angle_sol = angle + \
        angle_between_points(
            corner_point, (corner_point[0] + 1, corner_point[1]), mid_edge.end)
    # Solution center point placed with distance = radius from corner point and rotated to the middle of mid_edge and e1
    solution_center1 = (corner_point[0] + radius*math.cos(math.radians(
        angle_sol)), corner_point[1] + radius*math.sin(math.radians(angle_sol)))

    # Find projection to check if the solution is on edge e1
    proj_onto_line = projection_scalar_for_line(solution_center1, e1)
    if proj_onto_line > 1.0 or proj_onto_line < 0.0:
        raise ValueError(
            "The projection is not on the line segment - something went wrong for a corner disk (concave case - 1)")

    # Find the distance from the center perpendicular to one of the edges
    solution_radius1 = dist_point_to_line(solution_center1, e1)

    #### Second disk ####

    # Average angle between e1, e2:
    angle = angle_between_points(corner_point, e2.end, mid_edge.end) / 2
    # Total angle from the x-axis for center point
    angle_sol = angle + \
        angle_between_points(
            corner_point, (corner_point[0] + 1, corner_point[1]), e2.end)
    # Solution center point placed with distance = radius from corner point and rotated to the middle of e2 and mid_edge
    solution_center2 = (corner_point[0] + radius*math.cos(math.radians(
        angle_sol)), corner_point[1] + radius*math.sin(math.radians(angle_sol)))

    # Find projection to check if the solution is on edge e2
    proj_onto_line = projection_scalar_for_line(solution_center2, e2)
    if proj_onto_line > 1.0 or proj_onto_line < 0.0:
        raise ValueError(
            "The projection is not on the line segment - something went wrong for a corner disk (concave case - 2)")

    # Find the distance from the center perpendicular to one of the edges
    solution_radius2 = dist_point_to_line(solution_center2, e2)

    # Check if the radii are the same
    if not math.isclose(solution_radius1, solution_radius2, abs_tol=1e-3):
        raise ValueError(
            "The radii are not the same for the two disks - something went wrong")

    if e1.hole:
        return Disk(center=solution_center1, radius=solution_radius1, hole=True), Disk(center=solution_center2, radius=solution_radius2, hole=True)
    else:
        return Disk(center=solution_center1, radius=solution_radius1), Disk(center=solution_center2, radius=solution_radius2)
