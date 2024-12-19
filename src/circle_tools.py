from edge import Edge
from disk import Disk
import math
import mpmath as mp


def intersection_point_proj_dd(d1: Disk, d2: Disk):
    """
    Function to find the point of intersection between disk d1 and the line from the center of c1 to the center of c2.

    Args:
        d1 (Disk): The first disk
        d2 (Disk): The second disk

    Returns:
        tuple[float, float]: The point of tangency between the two disks
    """

    # Unpack
    (x1, y1), r1 = d1.center, d1.radius
    (x2, y2), r2 = d2.center, d2.radius

    distance = mp.sqrt((x2 - x1)**2 + (y2 - y1)**2)

    # Compute xt and yt
    xt = x1 + r1 / distance * (x2 - x1)
    yt = y1 + r1 / distance * (y2 - y1)

    return (xt, yt)


def intersection_point_proj_dl(d: Disk, e: Edge):
    """
    Function to find the point of intersection between disk d and the projection line from the center of d to the edge e.

    Args:
        d (Disk): The disk
        e (Edge): The line

    Returns:
        tuple[float, float]: The point of tangency between the disk and the line"""

    # Unpack
    (x, y), r = d.center, d.radius

    # Project the circle center (x, y) onto the line
    xt, yt = projection_point_on_line((x, y), e)

    # Compute the direction vector from circle's center to projected point
    dx = xt - x
    dy = yt - y

    # Find the length of the direction vector
    distance = mp.sqrt(dx**2 + dy**2)

    # Compute xt, yt
    if dx == 0:
        xt = x
    else:
        xt = x + r / distance * dx
    if dy == 0:
        yt = y
    else:
        yt = y + r / distance * dy

    return (xt, yt)


def tangent_point_dd(d1: Disk, d2: Disk):
    """
    Function to find the point of tangency between disk d1 and disk d2.
     - It is assumed that the disks are tangent to each other.

    Args:
        d1 (Disk): The first disk
        d2 (Disk): The second disk

    Returns:
        tuple[float, float]: The point of tangency between the two disks
    """

    xt, yt = intersection_point_proj_dd(d1, d2)

    return (xt, yt)


def tangent_point_ld(e: Edge, d: Disk):
    """
    Function to find the point of tangency between edge e and disk d.
        - It is assumed that the edge and the disk are tangent to each other.

    Args:
        e (Edge): The edge
        d (Disk): The disk

    Returns:
        tuple[float, float]: The point of tangency between the edge and the disk
    """

    xt, yt = intersection_point_proj_dl(d, e)

    return (xt, yt)


def tangent_point_between_two_elements(elm1, elm2):
    """
    Function to find the point of tangency between two elements.

    Args:
        elm1: The first element
        elm2: The second element

    Returns:
        tuple[float, float]: The point of tangency between the two elements
    """

    if isinstance(elm1, Disk) and isinstance(elm2, Disk):
        return tangent_point_dd(elm1, elm2)
    elif isinstance(elm1, Disk) and isinstance(elm2, Edge):
        return tangent_point_ld(elm2, elm1)
    elif isinstance(elm1, Edge) and isinstance(elm2, Disk):
        return tangent_point_ld(elm1, elm2)
    else:
        raise ValueError(
            "The two elements are not a valid combination of elements.")


def projection_point_on_line(p: tuple[float, float], e: Edge):
    """
    Function to find the projection of a point p onto a line e.

    Args:s
        p (tuple[float, float]): The point to project
        e (Edge): The line to project onto

    Returns:
        tuple[float, float]: The projected point
    """

    # Unpack
    (x1, y1), (x2, y2) = (e.start, e.end)

    # Projection of p onto the line
    t = projection_scalar_for_line(p, e)

    # Compute xt and yt
    xt = (1 - t) * x1 + t * x2
    yt = (1 - t) * y1 + t * y2

    x1 - t*x1 + t*x2

    return (xt, yt)


def projection_scalar_for_line(p: tuple[float, float], e: Edge):
    """
    Function to find the projection of a point p onto a line e.

    Args:
        p (tuple[float, float]): The point to project
        e (Edge): The line to project onto

    Returns:
        float: The projection of p onto the line
    """
    # Unpack
    (x, y) = p
    (x1, y1), (x2, y2) = (e.start, e.end)

    return ((x - x1) * (x2 - x1) + (y - y1) * (y2 - y1)) / ((x2 - x1)**2 + (y2 - y1)**2)


def dist_point_to_point(p1: tuple[float, float], p2: tuple[float, float]):
    """
    Function to find the distance between two points p1 and p2.

    Args:
        p1 (tuple[float, float]): The first point
        p2 (tuple[float, float]): The second point

    Returns:
        float: The distance between the two points
    """

    return mp.sqrt((p1[0] - p2[0]) ** 2 + (p1[1] - p2[1]) ** 2)


def dist_point_to_line(p: tuple[float, float], e: Edge):
    """
    Function to find the distance between a point p and an edge e.
        - It will find the distance to the nearest point on the edge

    Args:
        p (tuple[float, float]): The point
        e (Edge): The line

    Returns:
        float: The distance between the point and the line
    """
    # Unpack
    (x1, y1), (x2, y2) = (e.start, e.end)

    # Find the projection of the point onto the line
    t = projection_scalar_for_line(p, e)

    if 1.0 <= t or t <= 0.0:
        return min(dist_point_to_point(p, (x1, y1)), dist_point_to_point(p, (x2, y2)))

    # Calculate the distance using the formula
    xp = (1 - t) * x1 + t * x2 #= x1 - t*x1 + t*x2 = x1 + t*(x2 - x1)
    yp = (1 - t) * y1 + t * y2 #= y1 - t*y1 + t*y2 = y1 + t*(y2 - y1)

    return dist_point_to_point(p, (xp, yp))


def cross_product(p1: tuple[float, float],
                  p2: tuple[float, float],
                  p3: tuple[float, float]):
    """
    Function to compute the cross product of two vectors p1 -> p2 and p1 -> p3.

    Args:
        p1 (tuple[float, float]): The first point
        p2 (tuple[float, float]): The second point
        p3 (tuple[float, float]): The third point

    Returns:
        float: The cross product of the two vectors
    """
    # Unpack
    (x1, y1) = p1
    (x2, y2) = p2
    (x3, y3) = p3

    return (x2 - x1) * (y3 - y1) - (y2 - y1) * (x3 - x1)


def dot_product(p1: tuple[float, float],
                p2: tuple[float, float],
                p3: tuple[float, float]):
    """
    Function to compute the dot product of two vectors p1 -> p2 and p1 -> p3.

    Args:
        p1 (tuple[float, float]): The first point
        p2 (tuple[float, float]): The second point
        p3 (tuple[float, float]): The third point

    Returns:
        float: The dot product of the two vectors
    """
    # Unpack
    (x1, y1) = p1
    (x2, y2) = p2
    (x3, y3) = p3

    return (x2 - x1) * (x3 - x1) + (y2 - y1) * (y3 - y1)

def angle_between_edges(e1: Edge, e2: Edge):
    if e1.start == e2.start:
        return angle_between_points(e1.start, e1.end, e2.end, True)
    elif e1.start == e2.end:
        return angle_between_points(e1.start, e1.end, e2.start, True)
    elif e1.end == e2.start:
        return angle_between_points(e1.end, e1.start, e2.end, True)
    elif e1.end == e2.end:
        return angle_between_points(e1.end, e1.start, e2.start, True)
    else:
        raise ValueError("The edges are not connected to each other.")


def angle_between_points(p0: tuple[float, float],
                         p1: tuple[float, float],
                         p2: tuple[float, float],
                         clockwise: bool = False):
    """
    Function to compute the angle (in degrees) between two vectors p0 -> p1 and p0 -> p2 in the counter-clockwise order (as default).

    Args:
        p0 (tuple[float, float]): The first point
        p1 (tuple[float, float]): The second point
        p2 (tuple[float, float]): The third point
        clockwise (bool): Whether to compute the angle in the clockwise direction

    Returns:
        float: The angle between the two vectors
    """

    if math.isclose(p1[0], p2[0], abs_tol=1e-9) and math.isclose(p1[1], p2[1], abs_tol=1e-9):
        return 0.0

    # Dot product of p0->p1 and p0->p2
    dot_prod = dot_product(p0, p1, p2)

    # Lengths of p0->p1 and p0->p2
    length_p1 = dist_point_to_point(p0, p1)
    length_p2 = dist_point_to_point(p0, p2)

    length_mul = length_p1 * length_p2
    # Compute angle using acos (this will return the angle between 0 and 180 degrees)
    if dot_prod / length_mul > 1.0:
        length_mul = dot_prod
    elif (dot_prod / length_mul) < -1.0:
        length_mul = -dot_prod

    # Compute angle in radians and degrees
    angle_radians = mp.acos(dot_prod / (length_mul))
    angle_degrees = mp.degrees(angle_radians)

    # Compute cross_product to determine the direction (counter-clockwise or clockwise)
    cross_prod = cross_product(p0, p1, p2)

    # If cross_product > 0 then it is counter-clockwise and vice versa
    # So e.g. if we wish to find the counter-clockwise angle (even if it is concave) then we will have to modify the angle if cross_product < 0
    if (not clockwise and cross_prod < 0) or (clockwise and cross_prod > 0):
        angle_degrees = 360 - angle_degrees

    return angle_degrees


def tangency_correct_part_of_element(solutions: list[Disk], elm_tup):
    """
    Function to find the correct solutions for a disk or edge to be tangent to the element in the middle of the tuple.

    Args:
        solutions (list[Disk]): The list of solutions
        elm_tup (tuple): The tuple of the element and its neighbors

    Returns:
        list[Disk]: The list of valid solutions
    """

    # Unpack
    elm = elm_tup[1]
    neighbor1 = elm_tup[0]
    neighbor2 = elm_tup[2]
    if elm == neighbor1 or elm == neighbor2:# or neighbor1 == neighbor2:
        return solutions

    # In the (per default) counter-clockwise order of the elements is:
    # neighbor1 -> elm -> neighbor2

    if isinstance(elm, Disk):
        # Compute tangent point for neighbor1 and the element
        if isinstance(neighbor1, Disk):
            p1 = tangent_point_dd(elm, neighbor1)
        else:
            p1 = tangent_point_ld(neighbor1, elm)

        # Compute tangent point for neighbor2 and the element
        if isinstance(neighbor2, Disk):
            p2 = tangent_point_dd(elm, neighbor2)
        else:
            p2 = tangent_point_ld(neighbor2, elm)

        valid_solutions = []
        # The solutions are only valid if they tangent the element on the part between tangent points for neighbor1 and neighbor2
        for sol in solutions:
            # Compute tangent point for the element and the solution
            ps = tangent_point_dd(sol, elm)

            # Compute the angles from element.center -> p1 to element.center -> ps and element.center -> p1 to element.center -> p2
            angle_sol = angle_between_points(
                elm.center, p1, ps, True)  # Angle from p1 to sol
            angle_p2 = angle_between_points(
                elm.center, p1, p2, True)  # Angle from p1 to p2
            if angle_sol < angle_p2:
                valid_solutions.append(sol)
            if math.isclose(angle_p2, 0, abs_tol=1e-9):
                valid_solutions.append(sol)
        return valid_solutions
    else:  # Edge
        # Compute tangent point for neighbor1 and the element
        if isinstance(neighbor1, Disk):
            p1 = tangent_point_ld(elm, neighbor1)
        else:
            raise ValueError(
                "Tangency_correct_side_3rd_element: two consecutive elements cannot be edges")
        # Compute tangent point for neighbor2 and the element
        if isinstance(neighbor2, Disk):
            p2 = tangent_point_ld(elm, neighbor2)
        else:
            raise ValueError(
                "Tangency_correct_side_3rd_element: two consecutive elements cannot be edges")

        valid_solutions = []
        # The solutions are only valid if they tangent the element on the part between tangent points for neighbor1 and neighbor2
        for sol in solutions:
            # Compute tangent point for the element and the solution
            ps = tangent_point_ld(elm, sol)

            # Compute the following distances
            dist_p1_ps = dist_point_to_point(p1, ps)  # Distance from p1 to sol
            dist_p2_ps = dist_point_to_point(p2, ps)  # Distance from p2 to sol
            dist_p1_p2 = dist_point_to_point(p1, p2)  # Distance from p1 to p2

            # If math.isclose happen it is because the disk is tangent on the wrong side of the edge
            if max(dist_p1_ps, dist_p2_ps) < dist_p1_p2 and not math.isclose(min(dist_p1_ps, dist_p2_ps), 0, abs_tol=1e-9):
                # print("dist: ", min(dist_p1_ps, dist_p2_ps))
                valid_solutions.append(sol)
        return valid_solutions


def not_intersect(solutions: list[Disk], e: Edge):
    """
    Function to find the solution disks that do not intersect the edge e.

    Args:
        solutions (list[Disk]): The list of solutions
        e (Edge): The edge

    Returns:
        list[Disk]: The list of valid solutions
    """
    valid_solutions = []
    for sol in solutions:
        dist_center_to_line = dist_point_to_line(sol.center, e)
        if dist_center_to_line > sol.radius:
            valid_solutions.append(sol)
    return valid_solutions

def correct_side_edge(solutions: list[Disk], e: Edge):
    """
    Function to find the solution disks that are on the correct side of the edge e.

    Args:
        solutions (list[Disk]): The list of solutions
        e (Edge): The edge

    Returns:
        list[Disk]: The list of valid solutions
    """
    valid_solutions = []
    for sol in solutions:
        if cross_product(e.start, e.end, sol.center) > 0:
            valid_solutions.append(sol)
    return valid_solutions


#### OLD FUNCTIONS ####

# def find_solution_between_line_segments(e1: Edge, e2: Edge, solutions: list[Disk]):
#     (x1, y1), (x2, y2) = (e1.start, e1.end)
#     (x4, y4) = e2.end
#     (x1, y1) = (x1 - x2, y1 - y2)
#     (x4, y4) = (x4 - x2, y4 - y2)
#     dir_edges = x1 * y4 - y1 * x4
#     feasible_solutions = []

#     for sol in solutions:
#         xs, ys = sol.center
#         xs, ys = xs - x2, ys - y2
#         dir_edge1_sol = x1 * ys - y1 * xs  # same for both?
#         dir_sol_edge2 = xs * y4 - ys * x4

#         if ((dir_edges < 0 and dir_edge1_sol < 0 and dir_sol_edge2 < 0) or
#                 (dir_edges > 0 and dir_edge1_sol > 0 and dir_sol_edge2 > 0)):
#             feasible_solutions.append(sol)

#     if len(feasible_solutions) > 0:
#         return feasible_solutions[0]
#     else:
#         return (-1, -1, 1)

# def not_tangent_at_same_point(solutions: list[Disk], d: Disk, e: Edge):
#     p = tangent_point_ld(e, d)
#     valid_solutions = []
#     for sol in solutions:
#         ps = tangent_point_ld(e, sol)

#         if not (math.isclose(p[0], ps[0], abs_tol=1e-3) and math.isclose(p[1], ps[1], abs_tol=1e-3)):
#             valid_solutions.append(sol)
#     return valid_solutions


########## TANGENCY FUNCTIONS ###########

# def tangency_correct_side_ddl(solutions: list[Disk], d1: Disk, d2: Disk, e: Edge):
#     (x1, y1) = d1.center
#     (x3, y3), (x4, y4) = (e.start, e.end)

#     valid_solutions = []
#     for sol in solutions:
#         (x, y) = sol.center

#         proj_onto_line2 = projection_onto_line((x1, y1), e)
#         proj_x = (1 - proj_onto_line2) * x3 + proj_onto_line2 * x4
#         proj_y = (1 - proj_onto_line2) * y3 + proj_onto_line2 * y4

#         solution_tangency_point = tangent_point_dd(d1, sol)

#         ang_between_projections = angle_between_points((x1, y1), solution_tangency_point, (proj_x, proj_y))

#         if ang_between_projections <= 180:
#             valid_solutions.append(sol)


#     return valid_solutions


# def tangency_correct_side_dl(solutions: list[Disk], d: Disk, e: Edge, dl: bool = True, mid_disk: Disk = None):
#     (x1, y1), r1 = (d.center, d.radius)
#     (x2, y2), (x3, y3) = (e.start, e.end)
#     valid_solutions = []
#     angle = 0.0
#     for sol in solutions:
#         (x, y), r = sol.center, sol.radius

#         tangent_point_proj_ed = intersection_point_proj_dl(d, e)
#         if mid_disk == None:
#             tangent_point = tangent_point_ld(e, d) # Should be equal to tangent_point_ed if no mid_disk
#             extra_angle = 0
#         else:
#             tangent_point = tangent_point_dd(mid_disk, d)

#             if dl: dir = 1
#             else: dir = -1
#             if cross_product((x1, y1), tangent_point, tangent_point_proj_ed)*dir > 0:
#                 extra_angle = angle_between_points((x1, y1), tangent_point, tangent_point_proj_ed, not dl) #HERE
#             else:
#                 extra_angle = -angle_between_points((x1, y1), tangent_point, tangent_point_proj_ed, dl) # HERE
#         sol_tangent_point = tangent_point_dd(d, sol)

#         angle = angle_between_points((x1, y1), tangent_point, sol_tangent_point, not dl) #HERE

#         feasible_angle = 180 + extra_angle

#         if angle < feasible_angle:
#             valid_solutions.append((sol, angle))
#     if len(valid_solutions) > 0:
#         valid_solution = min(valid_solutions, key=lambda x: x[1])[0]
#     else:
#         return []
#     return [valid_solution]

# def tangency_correct_side_ld(solutions: list[Disk], e: Edge, d: Disk, mid_disk: Disk = None):
#     return tangency_correct_side_dl(solutions, d, e, False, mid_disk)


# def tangency_correct_side_dd(solutions: list[Disk], d1: Disk, d2: Disk, mid_elm = None):
#     (x1, y1), r1 = (d1.center, d1.radius)
#     (x2, y2) = d2.center
#     valid_solutions = []
#     angle = 0.0
#     for sol in solutions:
#         (x, y), r = sol.center, sol.radius

#         tangent_point_proj_d_d = intersection_point_proj_dd(d1, d2)
#         if mid_elm == None:
#             tangent_point = tangent_point_dd(d1, d2) # Should be equal to tangent_point_d_d if no mid_elm
#             extra_angle = 0
#         else:
#             if isinstance(mid_elm, Edge):
#                 tangent_point = tangent_point_ld(mid_elm, d1)
#             else:
#                 tangent_point = tangent_point_dd(mid_elm, d1)

#             if cross_product((x1, y1), tangent_point, tangent_point_proj_d_d) > 0:
#                 extra_angle = angle_between_points((x1, y1), tangent_point, tangent_point_proj_d_d, False) #HERE
#             else:
#                 extra_angle = -angle_between_points((x1, y1), tangent_point, tangent_point_proj_d_d, True) #HERE

#         # First circle - should be sufficient
#         new_tangent_point = tangent_point_dd(d1, sol)

#         angle = angle_between_points((x1, y1), tangent_point, new_tangent_point, False) #HERE
#         feasible_angle = 180 + extra_angle

#         if angle < feasible_angle:
#             valid_solutions.append((sol, angle))

#     if len(valid_solutions) > 0:
#         valid_solution = min(valid_solutions, key=lambda x: x[1])[0]
#     else:
#         return []
#     return [valid_solution]

# def tangency_correct_side_ll(solutions: list[Disk], e1: Edge, e2: Edge, mid_elm: Disk):
#     (x1, y1), (x2, y2) = (e1.start, e1.end)
#     (x3, y3), (x4, y4) = (e2.start, e2.end)
#     (x5, y5) = mid_elm.center

#     valid_solutions = []
#     for sol in solutions:
#         (x, y) = sol.center
#         proj_center_onto_line1 = projection_onto_line((x5, y5), e1)
#         proj_sol_onto_e1 = projection_onto_line((x, y), e1)

#         if proj_center_onto_line1 > proj_sol_onto_e1:
#             valid_solutions.append(sol)

#     return valid_solutions

# def tangency_correct_side_ldl(solutions: list[Disk], e1: Edge, d1: Disk, e2: Edge):
#     (x1, y1) = d1.center
#     (x2, y2), (x3, y3) = (e1.start,e1.end)
#     (x4, y4), (x5, y5) = (e2.start,e2.end)


#     valid_solutions = []
#     for sol in solutions:
#         (x, y) = sol.center

#         proj_onto_line2 = projection_onto_line((x1, y1), e2)
#         proj_x2 = (1 - proj_onto_line2) * x4 + proj_onto_line2 * x5
#         proj_y2 = (1 - proj_onto_line2) * y4 + proj_onto_line2 * y5

#         solution_tangency_point = tangent_point_dd(d1, sol)

#         ang_between_projections = angle_between_points((x1, y1), (proj_x2, proj_y2), solution_tangency_point)

#         if ang_between_projections <= 180:
#             valid_solutions.append(sol)

#     return valid_solutions
