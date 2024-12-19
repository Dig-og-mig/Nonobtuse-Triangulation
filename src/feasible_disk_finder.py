from disk import Disk
from edge import Edge
from circle_tools import angle_between_points, tangent_point_between_two_elements, dist_point_to_point
import mpmath as mp

tolerance = 1e-9


def feasible_disk_finder(solutions, solution_tup, e1, e2, mid_elm=None):

    # Unpack
    solution_disk, solution_radius, solution_dist = solution_tup

    if solution_disk == None:
        solution_radius = mp.inf
        solution_dist = mp.inf

    found_new_solution = False

    for sol in solutions:

        dist = solution_max_distance(e1, e2, mid_elm, sol)
        # if math.isclose(dist, 0, abs_tol=tolerance):
        #     continue

        if not mid_elm == None and dist <= solution_dist:
            found_new_solution = True
            # Solution tuple
            solution_disk = sol
            solution_radius = sol.radius
            solution_dist = dist

        elif mid_elm == None and sol.radius <= solution_radius:
            found_new_solution = True
            # Solution tuple
            solution_disk = sol
            solution_radius = sol.radius
            solution_dist = dist

    return found_new_solution, (solution_disk, solution_radius, solution_dist)


def solution_max_distance(e1, e2, mid_elm, sol):
    if mid_elm == None:
        return mp.inf

    tangent_point_e1_sol = tangent_point_between_two_elements(e1, sol)
    tangent_point_e2_sol = tangent_point_between_two_elements(e2, sol)

    tangent_point_e1_mid = tangent_point_between_two_elements(e1, mid_elm)
    tangent_point_e2_mid = tangent_point_between_two_elements(e2, mid_elm)

    if isinstance(e1, Disk) and isinstance(e2, Disk):
        angle_e1 = angle_between_points(
            e1.center, tangent_point_e1_mid, tangent_point_e1_sol, False)  # HERE
        angle_e2 = angle_between_points(
            e2.center, tangent_point_e2_mid, tangent_point_e2_sol, True)  # HERE

        dist_e1 = angle_e1/360 * 2 * mp.pi * e1.radius
        dist_e2 = angle_e2/360 * 2 * mp.pi * e2.radius

        dist = max(dist_e1, dist_e2)
        return dist
    elif isinstance(e1, Disk) and isinstance(e2, Edge):
        angle_e1 = angle_between_points(
            e1.center, tangent_point_e1_mid, tangent_point_e1_sol, False)  # HERE
        dist_e1 = angle_e1/360 * 2 * mp.pi * e1.radius
        dist_e2 = dist_point_to_point(
            tangent_point_e2_sol, tangent_point_e2_mid)  # e2 == edge

        dist = max(dist_e1, dist_e2)
        return dist
    elif isinstance(e1, Edge) and isinstance(e2, Disk):
        angle_e2 = angle_between_points(
            e2.center, tangent_point_e2_mid, tangent_point_e2_sol, True)  # HERE
        dist_e1 = dist_point_to_point(
            tangent_point_e1_sol, tangent_point_e1_mid)  # e1 == edge
        dist_e2 = angle_e2/360 * 2 * mp.pi * e2.radius

        dist = max(dist_e1, dist_e2)
        return dist
    elif isinstance(e1, Edge) and isinstance(e2, Edge):
        dist_e1 = dist_point_to_point(
            tangent_point_e1_sol, tangent_point_e1_mid)  # e1 == edge
        dist_e2 = dist_point_to_point(
            tangent_point_e2_sol, tangent_point_e2_mid)  # e2 == edge

        dist = max(dist_e1, dist_e2)
        return dist
