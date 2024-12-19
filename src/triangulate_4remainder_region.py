from disk import Disk
from edge import Edge
from triangulation_tools import find_inscribed_circle_to_Rp, find_perpendicular_line_l1_p1, find_tangent_points_of_Rp
from special_cases_checks import center_is_inside_convex_hull, arcs_less_than_180
from circle_tools import dist_point_to_point, dist_point_to_line, projection_point_on_line, tangent_point_between_two_elements, cross_product
import math
import plotting
import numpy as np

tolerance = 1e-9


def triangulate_4remainder_region(boundary: list, c_star=None):
    count_edges = 0
    for elm in boundary:
        if isinstance(elm, Edge):
            count_edges += 1
    if count_edges > 2:
        raise ValueError("More than 2 edges in 4remainder region")

    if not c_star == None:
        # tangent_0_1 = tangent point between c* and c1
        # tangent_1_2 = tangent point between c1 and c4
        # tangent_2_3 = tangent point between c4 and c3
        # tangent_3_0 = tangent point between c3 and c*
        c_star, c1, c4, c3 = boundary
        tangent_0_1, tangent_1_2, tangent_2_3, tangent_3_0 = find_tangent_points_of_Rp(
            boundary)

        # Line for mutual chord between c* and c4
        if isinstance(c4, Edge):
            if math.isclose(c4.end[0], c4.start[0], abs_tol=tolerance):
                line_common_chord = (True, c4.start[0])
            else:
                l_m_chord = (c4.end[1] - c4.start[1]) / \
                    (c4.end[0] - c4.start[0])
                l_b_chord = c4.start[1] - l_m_chord * c4.start[0]
                line_common_chord = (False, (l_m_chord, l_b_chord))
        else:  # c4 is a Disk
            A = 2 * (c_star.center[0] - c4.center[0])
            B = 2 * (c_star.center[1] - c4.center[1])
            C = (c4.center[0]**2 - c_star.center[0]**2) + (c4.center[1] **
                                                           2 - c_star.center[1]**2) + (c_star.radius**2 - c4.radius**2)
            if math.isclose(B, 0, abs_tol=tolerance):
                line_common_chord = (True, -C/A)
            else:
                l_m_chord = -A / B
                l_b_chord = -C / B
                line_common_chord = (False, (l_m_chord, l_b_chord))

        # tangent line for c*, c1:
        if isinstance(c1, Edge):
            if math.isclose(c1.end[0], c1.start[0], abs_tol=tolerance):
                line_c_1 = (True, c1.start[0])
            else:
                l_m_c_1 = (c1.end[1] - c1.start[1]) / \
                    (c1.end[0] - c1.start[0])
                l_b_c_1 = c1.start[1] - l_m_c_1 * c1.start[0]
                line_c_1 = (False, (l_m_c_1, l_b_c_1))
        else:
            if math.isclose(c_star.center[1], c1.center[1], abs_tol=tolerance):
                line_c_1 = (True, tangent_0_1)
            else:
                l_star_m, l_star_b = find_perpendicular_line_l1_p1(
                    (c_star.center, c1.center), tangent_0_1)
                line_c_1 = (False, (l_star_m, l_star_b))

        # tangent line for c*, c3:
        if isinstance(c3, Edge):
            if math.isclose(c3.end[0], c3.start[0], abs_tol=tolerance):
                line_c_2 = (True, c3.start[0])
            else:
                l_m_c_2 = (c3.end[1] - c3.start[1]) / \
                    (c3.end[0] - c3.start[0])
                l_b_c_2 = c3.start[1] - l_m_c_2 * c3.start[0]
                line_c_2 = (False, (l_m_c_2, l_b_c_2))
        else:
            if math.isclose(c_star.center[1], c3.center[1], abs_tol=tolerance):
                line_c_2 = (True, tangent_3_0)
            else:
                l_star_m, l_star_b = find_perpendicular_line_l1_p1(
                    (c_star.center, c3.center), tangent_3_0)
                line_c_2 = (False, (l_star_m, l_star_b))

        # Find c_1
        if not line_common_chord[0] and not line_c_1[0]:  # Normal case
            c_1_x = (line_common_chord[1][1] - line_c_1[1][1]) / \
                (line_c_1[1][0] - line_common_chord[1][0])
            c_1_y = line_c_1[1][0] * c_1_x + line_c_1[1][1]
        # line_common_chord is vertical
        elif line_common_chord[0] and not line_c_1[0]:
            c_1_x = line_common_chord[1]
            c_1_y = line_c_1[1][0] * c_1_x + line_c_1[1][1]
        elif not line_common_chord[0] and line_c_1[0]:  # line_c_1 is vertical
            c_1_x = line_c_1[1]
            c_1_y = line_common_chord[1][0] * \
                c_1_x + line_common_chord[1][1]
        else:  # Both lines are vertical
            raise ValueError(
                "Both common chord and the tangent line between c* and c1 are vertical")
        c_1 = (c_1_x, c_1_y)

        # Find c_2
        if not line_common_chord[0] and not line_c_2[0]:  # Normal case
            c_2_x = (line_common_chord[1][1] - line_c_2[1][1]) / \
                (line_c_2[1][0] - line_common_chord[1][0])
            c_2_y = line_c_2[1][0] * c_2_x + line_c_2[1][1]
        # line_common_chord is vertical
        elif line_common_chord[0] and not line_c_2[0]:
            c_2_x = line_common_chord[1]
            c_2_y = line_c_2[1][0] * c_2_x + line_c_2[1][1]
        elif not line_common_chord[0] and line_c_2[0]:  # line_c_2 is vertical
            c_2_x = line_c_2[1]
            c_2_y = line_common_chord[1][0] * \
                c_2_x + line_common_chord[1][1]
        else:  # Both lines are vertical
            raise ValueError(
                "Both common chord and the tangent line between c* and c3 are vertical")
        c_2 = (c_2_x, c_2_y)

        # Find c_mid (perpendicular line from c* to the common chord)
        if not line_common_chord[0]:  # normal case
            c_mid = projection_point_on_line(
                c_star.center, Edge(start=c_1, end=c_2))
        else:  # common chord is vertical
            c_mid = (line_common_chord[1], c_star.center[1])

        ### MAKE EDGES ###
        edges_to_plot = []
        # Edges convex hull
        e1 = Edge(start=c_star.center, end=tangent_0_1)
        edges_to_plot.append(e1)
        if isinstance(c1, Disk):
            e2 = Edge(start=tangent_0_1, end=c1.center)
            e3 = Edge(start=c1.center, end=tangent_1_2)
            edges_to_plot += [e2, e3]
        if isinstance(c4, Disk):
            e4 = Edge(start=tangent_1_2, end=c4.center)
            e5 = Edge(start=c4.center, end=tangent_2_3)
            edges_to_plot += [e4, e5]
        if isinstance(c3, Disk):
            e6 = Edge(start=tangent_2_3, end=c3.center)
            e7 = Edge(start=c3.center, end=tangent_3_0)
            edges_to_plot += [e6, e7]
        e8 = Edge(start=tangent_3_0, end=c_star.center)
        edges_to_plot.append(e8)

        # Edges from c*
        e9 = Edge(start=c_star.center, end=c_1)
        e10 = Edge(start=c_star.center, end=c_mid)
        e11 = Edge(start=c_star.center, end=c_2)
        edges_to_plot += [e9, e10, e11]

        # Edges from c1
        if isinstance(c1, Disk):
            e12 = Edge(start=c1.center, end=c_1)
            edges_to_plot.append(e12)

        # Edges from c4
        if isinstance(c4, Disk):
            e13 = Edge(start=c4.center, end=c_1)
            e14 = Edge(start=c4.center, end=c_mid)
            e15 = Edge(start=c4.center, end=c_2)
            edges_to_plot += [e13, e14, e15]

        # Edges from c3
        if isinstance(c3, Disk):
            e16 = Edge(start=c3.center, end=c_2)
            edges_to_plot.append(e16)

        # Edges from tangent c*,c1
        if isinstance(c1, Disk):
            e17 = Edge(start=tangent_0_1, end=c_1)
            edges_to_plot.append(e17)

        # Edges from tangent c1,c4
        if isinstance(c1, Disk):
            e18 = Edge(start=tangent_1_2, end=c_1)
            edges_to_plot.append(e18)

        # Edges from tangent c4,c3
        if isinstance(c3, Disk):
            e19 = Edge(start=tangent_2_3, end=c_2)
            edges_to_plot.append(e19)

        # Edges from tangent c3,c*
        if isinstance(c3, Disk):
            e20 = Edge(start=tangent_3_0, end=c_2)
            edges_to_plot.append(e20)

        # Remaining edges from c_1
        e21 = Edge(start=c_1, end=c_mid)
        e22 = Edge(start=c_mid, end=c_2)
        edges_to_plot += [e21, e22]

        plotting.plot_edges(edges_to_plot, color='black')

        return edges_to_plot, 12, 1

    d = find_inscribed_circle_to_Rp(boundary)
    c = d.center
    center_inside, elm_center_outside = center_is_inside_convex_hull(
        c, boundary)
    small_arcs, big_arc = arcs_less_than_180(boundary)
    if (center_inside and small_arcs):
        if isinstance(boundary[0], Disk):
            p0 = boundary[0].center
        else:
            p0 = projection_point_on_line(c, boundary[0])
        if isinstance(boundary[1], Disk):
            p1 = boundary[1].center
        else:
            p1 = projection_point_on_line(c, boundary[1])
        if isinstance(boundary[2], Disk):
            p2 = boundary[2].center
        else:
            p2 = projection_point_on_line(c, boundary[2])
        if isinstance(boundary[3], Disk):
            p3 = boundary[3].center
        else:
            p3 = projection_point_on_line(c, boundary[3])

        e0_helper = Edge(start=p0, end=c)
        e1_helper = Edge(start=p1, end=c)
        e2_helper = Edge(start=p2, end=c)
        e3_helper = Edge(start=p3, end=c)

        tangent_0_1, tangent_1_2, tangent_2_3, tangent_3_0 = find_tangent_points_of_Rp(
            boundary)
        q1 = projection_point_on_line(tangent_0_1, e1_helper)
        e1 = Edge(start=tangent_0_1, end=q1)
        e2 = Edge(start=q1, end=tangent_1_2)
        e3 = Edge(start=p1, end=q1)
        e4 = Edge(start=q1, end=c)

        q2 = projection_point_on_line(tangent_1_2, e2_helper)
        e5 = Edge(start=tangent_1_2, end=q2)
        e6 = Edge(start=q2, end=tangent_2_3)
        e7 = Edge(start=p2, end=q2)
        e8 = Edge(start=q2, end=c)

        q3 = projection_point_on_line(tangent_2_3, e3_helper)
        e9 = Edge(start=tangent_2_3, end=q3)
        e10 = Edge(start=q3, end=tangent_3_0)
        e11 = Edge(start=p3, end=q3)
        e12 = Edge(start=q3, end=c)

        q4 = projection_point_on_line(tangent_3_0, e0_helper)
        e13 = Edge(start=tangent_3_0, end=q4)
        e14 = Edge(start=q4, end=tangent_0_1)
        e15 = Edge(start=p0, end=q4)
        e16 = Edge(start=q4, end=c)

        e17 = Edge(start=tangent_0_1, end=c)
        e18 = Edge(start=tangent_1_2, end=c)
        e19 = Edge(start=tangent_2_3, end=c)
        e20 = Edge(start=tangent_3_0, end=c)

        e21 = Edge(start=p0, end=tangent_0_1)
        e22 = Edge(start=tangent_0_1, end=p1)
        e23 = Edge(start=p1, end=tangent_1_2)
        e24 = Edge(start=tangent_1_2, end=p2)
        e25 = Edge(start=p2, end=tangent_2_3)
        e26 = Edge(start=tangent_2_3, end=p3)
        e27 = Edge(start=p3, end=tangent_3_0)
        e28 = Edge(start=tangent_3_0, end=p0)

        plotting.plot_edges([e1, e2, e3, e4, e5, e6, e7, e8, e9, e10, e11, e12, e13,
                            e14, e15, e16, e17, e18, e19, e20, e21, e22, e23, e24, e25, e26, e27, e28])
        return [e1, e2, e3, e4, e5, e6, e7, e8, e9, e10, e11, e12, e13, e14, e15, e16, e17, e18, e19, e20, e21, e22, e23, e24, e25, e26, e27, e28], 16, 1

    elif (not small_arcs):
        c3 = big_arc  # has to be a disk
        index_c3 = boundary.index(big_arc)
        c2 = boundary[(index_c3+1) % 4]
        c1 = boundary[(index_c3+2) % 4]
        c4 = boundary[(index_c3+3) % 4]

        if isinstance(c1, Disk):
            line_c3_c1 = Edge(start=c3.center, end=c1.center)
            dist_c3_c1 = dist_point_to_point(c3.center, c1.center)
            c_star_radius = (dist_c3_c1 - c3.radius - c1.radius) / 2
            ratio = (c3.radius + c_star_radius) / dist_c3_c1
        else:
            line_c3_c1 = Edge(
                start=c3.center, end=projection_point_on_line(c3.center, c1))
            dist_c3_c1 = dist_point_to_point(c3.center, line_c3_c1.end)
            c_star_radius = (dist_c3_c1 - c3.radius) / 2
            ratio = (c3.radius + c_star_radius) / dist_c3_c1

        xp = (1 - ratio) * line_c3_c1.start[0] + ratio * line_c3_c1.end[0]
        yp = (1 - ratio) * line_c3_c1.start[1] + ratio * line_c3_c1.end[1]
        c_star = Disk(center=(xp, yp), radius=c_star_radius)

        plotting.plot_disk(c_star, color='cyan')

        # c3, c_star, c1, c4
        overlapping_star = None
        if isinstance(c4, Disk):
            if (c4.radius + c_star.radius) > dist_point_to_point(c4.center, c_star.center):
                overlapping_star = c_star
        else:  # c4 is an edge
            if (c_star.radius) > dist_point_to_line(c_star.center, c4):
                overlapping_star = c_star
        edges1, t1, _ = triangulate_4remainder_region(
            [c_star, c1, c4, c3], c_star=overlapping_star)

        # c3, c2, c1, c_star
        overlapping_star = None
        if isinstance(c2, Disk):
            if (c2.radius + c_star.radius) > dist_point_to_point(c2.center, c_star.center):
                overlapping_star = c_star
        else:  # c2 is an edge
            if (c_star.radius) > dist_point_to_line(c_star.center, c2):
                overlapping_star = c_star
        edges2, t2, _ = triangulate_4remainder_region(
            [c_star, c3, c2, c1], c_star=overlapping_star)
        return edges1+edges2, t1+t2, 2
    elif (not center_inside):
        # Find tangency points:
        # boundary.index(elm_center_outside)
        new_boundary_order = [
            boundary[(i + boundary.index(elm_center_outside)) % 4] for i in range(4)]
        t41, t34, t23, t12 = find_tangent_points_of_Rp(new_boundary_order)

        c1 = new_boundary_order[0]
        c4 = new_boundary_order[1]
        c3 = new_boundary_order[2]
        c2 = new_boundary_order[3]

        # Find cr and cl

        # FIRST CASE - c_star is places 'vertically' between c1 and c3
        # - note that c3 will always be a disk
        if isinstance(c1, Disk):
            line_c3_c1 = Edge(start=c3.center, end=c1.center)
            dist_c3_c1 = dist_point_to_point(c3.center, c1.center)
            c_star_radius = (dist_c3_c1 - c3.radius - c1.radius) / 2
            ratio = (c3.radius + c_star_radius) / dist_c3_c1
        else:
            line_c3_c1 = Edge(
                start=c3.center, end=projection_point_on_line(c3.center, c1))
            dist_c3_c1 = dist_point_to_point(c3.center, line_c3_c1.end)
            c_star_radius = (dist_c3_c1 - c3.radius) / 2
            ratio = (c3.radius + c_star_radius) / dist_c3_c1

        xp = (1 - ratio) * line_c3_c1.start[0] + ratio * line_c3_c1.end[0]
        yp = (1 - ratio) * line_c3_c1.start[1] + ratio * line_c3_c1.end[1]
        c_star = Disk(center=(xp, yp), radius=c_star_radius)
        ts1 = tangent_point_between_two_elements(c_star, c1)
        ts3 = tangent_point_between_two_elements(c_star, c3)

        # cl:
        (x1, y1) = ts3
        (x2, y2) = t23
        (x3, y3) = t12
        M11 = 2*(x3-x1)
        M12 = 2*(y3-y1)
        M21 = 2*(x3-x2)
        M22 = 2*(y3-y2)

        b1 = x3**2 - x1**2 + y3**2 - y1**2
        b2 = x3**2 - x2**2 + y3**2 - y2**2

        if (math.isclose(M11, 0, abs_tol=tolerance)):
            clx = (M12*b2 - M22*b1) / (M12*M21)
            cly = b1/M12
        elif (math.isclose(M21, 0, abs_tol=tolerance)):
            clx = -(M12*b2 - M22*b1) / (M22*M11)
            cly = b2/M22
        else:
            clx = - (M12*b2 - M22*b1) / (M22*M11-M21*M12)
            cly = (b2*M11 - M21*b1) / (M22*M11-M21*M12)
        clr = dist_point_to_point((clx, cly), ts3)

        # cr:
        (x1, y1) = ts3
        (x2, y2) = t34
        (x3, y3) = t41
        M11 = 2*(x3-x1)
        M12 = 2*(y3-y1)
        M21 = 2*(x3-x2)
        M22 = 2*(y3-y2)

        b1 = x3**2 - x1**2 + y3**2 - y1**2
        b2 = x3**2 - x2**2 + y3**2 - y2**2

        if (math.isclose(M11, 0, abs_tol=tolerance)):
            crx = (M12*b2 - M22*b1) / (M12*M21)
            cry = b1/M12
        elif (math.isclose(M21, 0, abs_tol=tolerance)):
            crx = -(M12*b2 - M22*b1) / (M22*M11)
            cry = b2/M22
        else:
            crx = - (M12*b2 - M22*b1) / (M22*M11-M21*M12)
            cry = (b2*M11 - M21*b1) / (M22*M11-M21*M12)
        crr = dist_point_to_point((crx, cry), ts3)
        plotting.plot_disk(c_star, color='red')

        # CHECK IF ONE OF THEM IS OUTSIDE, AND THEN SWEEP
        cr_outside = cross_product(t41, t34, (crx, cry)) < 0
        cl_outside = cross_product(t23, t12, (crx, cry)) < 0

        if cr_outside:
            (x1, y1) = t41
            (x2, y2) = t34

            # cr:
            M11 = y1-y2
            M12 = x2-x1
            M21 = 2*(x2-x1)
            M22 = 2*(y2-y1)

            b1 = -(y2-y1)*x1 + (x2-x1)*y1
            b2 = x2**2 - x1**2 + y2**2 - y1**2

            if math.isclose(M22*M11-M21*M12, 0, abs_tol=tolerance):
                raise ValueError("This should not happen, M22: {}, M11: {}, M21: {}, M12: {}".format(
                    M22, M11, M21, M12))
            crx = - (M12*b2 - M22*b1) / (M22*M11-M21*M12)
            cry = (b2*M11 - M21*b1) / (M22*M11-M21*M12)
            crr = dist_point_to_point((crx, cry), t41)

            intersect_tup = circle_intersect(
                crx, cry, crr, new_boundary_order[2].center[0], new_boundary_order[2].center[1], new_boundary_order[2].radius)

            if math.isclose(intersect_tup[0][0], t34[0], abs_tol=tolerance) and math.isclose(intersect_tup[0][1], t34[1], abs_tol=tolerance):
                ts3 = (intersect_tup[1][0], intersect_tup[1][1])
            else:
                ts3 = (intersect_tup[0][0], intersect_tup[0][1])

            # cl:
            (x1, y1) = ts3
            (x2, y2) = t23
            (x3, y3) = t12
            M11 = 2*(x3-x1)
            M12 = 2*(y3-y1)
            M21 = 2*(x3-x2)
            M22 = 2*(y3-y2)

            b1 = x3**2 - x1**2 + y3**2 - y1**2
            b2 = x3**2 - x2**2 + y3**2 - y2**2

            if (math.isclose(M11, 0, abs_tol=tolerance)):
                clx = (M12*b2 - M22*b1) / (M12*M21)
                cly = b1/M12
            elif (math.isclose(M21, 0, abs_tol=tolerance)):
                clx = -(M12*b2 - M22*b1) / (M22*M11)
                cly = b2/M22
            else:
                clx = - (M12*b2 - M22*b1) / (M22*M11-M21*M12)
                cly = (b2*M11 - M21*b1) / (M22*M11-M21*M12)
            clr = dist_point_to_point((clx, cly), ts3)

        elif cl_outside:
            (x1, y1) = t23
            (x2, y2) = t12

            # cr:
            M11 = y1-y2
            M12 = x2-x1
            M21 = 2*(x2-x1)
            M22 = 2*(y2-y1)

            b1 = -(y2-y1)*x1 + (x2-x1)*y1
            b2 = x2**2 - x1**2 + y2**2 - y1**2

            if math.isclose(M22*M11-M21*M12, 0, abs_tol=tolerance):
                raise ValueError("This should not happen, M22: {}, M11: {}, M21: {}, M12: {}".format(
                    M22, M11, M21, M12))
            clx = - (M12*b2 - M22*b1) / (M22*M11-M21*M12)
            cly = (b2*M11 - M21*b1) / (M22*M11-M21*M12)
            clr = dist_point_to_point((crx, cry), t41)

            intersect_tup = circle_intersect(
                clx, cly, clr, new_boundary_order[2].center[0], new_boundary_order[2].center[1], new_boundary_order[2].radius)

            if math.isclose(intersect_tup[0][0], t34[0], abs_tol=tolerance) and math.isclose(intersect_tup[0][1], t34[1], abs_tol=tolerance):
                ts3 = (intersect_tup[1][0], intersect_tup[1][1])
            else:
                ts3 = (intersect_tup[0][0], intersect_tup[0][1])

            # cl:
            (x1, y1) = ts3
            (x2, y2) = t41
            (x3, y3) = t34
            M11 = 2*(x3-x1)
            M12 = 2*(y3-y1)
            M21 = 2*(x3-x2)
            M22 = 2*(y3-y2)

            b1 = x3**2 - x1**2 + y3**2 - y1**2
            b2 = x3**2 - x2**2 + y3**2 - y2**2

            if (math.isclose(M11, 0, abs_tol=tolerance)):
                crx = (M12*b2 - M22*b1) / (M12*M21)
                cry = b1/M12
            elif (math.isclose(M21, 0, abs_tol=tolerance)):
                crx = -(M12*b2 - M22*b1) / (M22*M11)
                cry = b2/M22
            else:
                crx = - (M12*b2 - M22*b1) / (M22*M11-M21*M12)
                cry = (b2*M11 - M21*b1) / (M22*M11-M21*M12)
            crr = dist_point_to_point((crx, cry), ts3)

        cl = Disk(center=(clx, cly), radius=clr)
        cr = Disk(center=(crx, cry), radius=crr)

        plotting.plot_elements([cl, cr], color='darkorange', linestyle='--')

        # triangulation:

        # Find ts1
        ts1_tup = circle_intersect(clx, cly, clr, crx, cry, crr)
        if math.isclose(ts1_tup[0][0], ts3[0], abs_tol=tolerance) and math.isclose(ts1_tup[0][1], ts3[1], abs_tol=tolerance):
            ts1 = (ts1_tup[1][0], ts1_tup[1][1])
        else:
            ts1 = (ts1_tup[0][0], ts1_tup[0][1])

        edges = []

        # C3 - always a disk
        c3_center = new_boundary_order[2].center
        e1 = Edge(start=c3_center, end=t23)
        e2 = Edge(start=c3_center, end=ts3)
        e3 = Edge(start=c3_center, end=t34)
        e4 = Edge(start=c3_center, end=projection_point_on_line(
            c3_center, Edge(start=ts3, end=t34)))
        e5 = Edge(start=c3_center, end=projection_point_on_line(
            c3_center, Edge(start=t23, end=ts3)))
        e6 = Edge(start=t23, end=projection_point_on_line(
            t23, Edge(start=c3_center, end=(clx, cly))))
        e7 = Edge(start=t34, end=projection_point_on_line(
            t34, Edge(start=c3_center, end=(crx, cry))))
        edges += [e1, e2, e3, e4, e5, e6, e7]

        # C1
        if isinstance(new_boundary_order[0], Disk):
            c1_center = new_boundary_order[0].center
            e6 = Edge(start=c1_center, end=t41)
            e7 = Edge(start=c1_center, end=ts1)
            e8 = Edge(start=c1_center, end=t12)
            e9 = Edge(start=c1_center, end=projection_point_on_line(
                c1_center, Edge(start=t41, end=ts1)))
            e10 = Edge(start=c1_center, end=projection_point_on_line(
                c1_center, Edge(start=ts1, end=t12)))
            e11 = Edge(start=ts1, end=projection_point_on_line(
                ts1, Edge(start=c1_center, end=(clx, cly))))
            e12 = Edge(start=ts1, end=projection_point_on_line(
                ts1, Edge(start=c1_center, end=(crx, cry))))
            e13 = Edge(start=t12, end=projection_point_on_line(
                cl.center, Edge(start=t12, end=ts1)))
            e14 = Edge(start=projection_point_on_line(
                cl.center, Edge(start=t12, end=ts1)), end=ts1)
            e15 = Edge(start=ts1, end=projection_point_on_line(
                cr.center, Edge(start=ts1, end=t41)))
            e16 = Edge(start=projection_point_on_line(
                cr.center, Edge(start=ts1, end=t41)), end=t41)
            edges += [e6, e7, e8, e9, e10, e11, e12, e13, e14, e15, e16]
        else:
            e6 = Edge(start=t12, end=ts1)
            e7 = Edge(start=ts1, end=t41)
            e8 = Edge(start=ts1, end=projection_point_on_line(
                (clx, cly), Edge(start=t12, end=t41)))
            e9 = Edge(start=projection_point_on_line(
                (clx, cly), Edge(start=t12, end=t41)), end=t12)
            e10 = Edge(start=ts1, end=projection_point_on_line(
                (crx, cry), Edge(start=t12, end=t41)))
            e11 = Edge(start=projection_point_on_line(
                (crx, cry), Edge(start=t12, end=t41)), end=t41)
            edges += [e6, e7, e8, e9, e10, e11]

        # C4
        if isinstance(new_boundary_order[1], Disk):
            c4_center = new_boundary_order[1].center
            e11 = Edge(start=c4_center, end=t34)
            e12 = Edge(start=c4_center, end=t41)
            e13 = Edge(start=c4_center, end=projection_point_on_line(
                c4_center, Edge(start=t34, end=t41)))
            e14 = Edge(start=projection_point_on_line(
                (crx, cry), Edge(start=t41, end=t34)), end=t41)
            e15 = Edge(start=projection_point_on_line(
                (crx, cry), Edge(start=t41, end=t34)), end=t34)
            edges += [e11, e12, e13, e14, e15]
        else:
            e11 = Edge(start=(crx, cry), end=t34)
            e12 = Edge(start=(crx, cry), end=t41)
            e13 = Edge(start=projection_point_on_line(
                (crx, cry), Edge(start=t41, end=t34)), end=t41)
            e14 = Edge(start=projection_point_on_line(
                (crx, cry), Edge(start=t41, end=t34)), end=t34)
            edges += [e11, e12, e13, e14]

        # C2
        if isinstance(new_boundary_order[3], Disk):
            c2_center = new_boundary_order[3].center
            e14 = Edge(start=c2_center, end=t23)
            e15 = Edge(start=c2_center, end=t12)
            e16 = Edge(start=c2_center, end=projection_point_on_line(
                c2_center, Edge(start=t23, end=t12)))
            e17 = Edge(start=projection_point_on_line(
                (clx, cly), Edge(start=t23, end=t12)), end=t23)
            e18 = Edge(start=projection_point_on_line(
                (clx, cly), Edge(start=t23, end=t12)), end=t12)
            edges += [e14, e15, e16, e17, e18]
        else:
            e14 = Edge(start=projection_point_on_line(
                (clx, cly), Edge(start=t23, end=t12)), end=t23)
            e15 = Edge(start=projection_point_on_line(
                (clx, cly), Edge(start=t23, end=t12)), end=t12)
            edges += [e14, e15]

        # cl
        e17 = Edge(start=cl.center, end=t23)
        e18 = Edge(start=cl.center, end=t12)
        e19 = Edge(start=cl.center, end=ts1)
        e20 = Edge(start=cl.center, end=ts3)
        e21 = Edge(start=cl.center, end=projection_point_on_line(
            cl.center, Edge(start=t23, end=t12)))
        e22 = Edge(start=cl.center, end=projection_point_on_line(
            cl.center, Edge(start=t12, end=ts1)))
        e23 = Edge(start=cl.center, end=projection_point_on_line(
            cl.center, Edge(start=ts1, end=ts3)))
        e24 = Edge(start=cl.center, end=projection_point_on_line(
            cl.center, Edge(start=ts3, end=t23)))
        edges += [e17, e18, e19, e20, e21, e22, e23, e24]

        # cr
        e25 = Edge(start=cr.center, end=t34)
        e26 = Edge(start=cr.center, end=t41)
        e27 = Edge(start=cr.center, end=ts1)
        e28 = Edge(start=cr.center, end=ts3)
        e29 = Edge(start=cr.center, end=projection_point_on_line(
            cr.center, Edge(start=t34, end=t41)))
        e30 = Edge(start=cr.center, end=projection_point_on_line(
            cr.center, Edge(start=t41, end=ts1)))
        e31 = Edge(start=cr.center, end=projection_point_on_line(
            cr.center, Edge(start=ts1, end=ts3)))
        e32 = Edge(start=cr.center, end=projection_point_on_line(
            cr.center, Edge(start=ts3, end=t34)))
        edges += [e25, e26, e27, e28, e29, e30, e31, e32]

        # Remaining edges
        e33 = Edge(start=ts3, end=projection_point_on_line(
            ts3, Edge(start=cl.center, end=cr.center)))
        e34 = Edge(start=ts1, end=projection_point_on_line(
            ts1, Edge(start=cl.center, end=cr.center)))
        e35 = Edge(start=ts3, end=projection_point_on_line(
            ts3, Edge(start=cl.center, end=c3.center)))
        e36 = Edge(start=ts3, end=projection_point_on_line(
            ts3, Edge(start=cr.center, end=c3.center)))
        edges += [e33, e34, e35, e36]

        plotting.plot_edges(edges, color='black')
        return edges, 28, 1
    else:
        raise ValueError("Case not implemented ??")


def circle_intersect(x0, y0, r0, x1, y1, r1):
    c0 = np.array([x0, y0])
    c1 = np.array([x1, y1])
    v = c1 - c0
    d = np.linalg.norm(v)

    if d > r0 + r1 or d == 0:
        return None

    u = v/np.linalg.norm(v)
    xvec = c0 + (d**2 - r1**2 + r0**2)*u/(2*d)

    uperp = np.array([u[1], -u[0]])
    a = ((-d+r1-r0)*(-d-r1+r0)*(-d+r1+r0)*(d+r1+r0))**0.5/d
    res = (xvec + a*uperp/2, xvec - a*uperp/2)
    return res
