from models.edge import Edge
from disk_packing.circle_tools import tangent_point_ld, tangent_point_dd
from triangulation.triangulation_tools import find_perpendicular_point_l1_p1_to_l2, find_inscribed_circle_to_Rp
import misc.plotting


def triangulate_3remainder_region(boundary: list):
    edges = list(filter(lambda x: isinstance(x, Edge), boundary))
    if len(edges) > 1:
        raise ValueError("More than 1 edge in 3remainder region")
    elif len(edges) == 1:
        if isinstance(boundary[0], Edge):
            e = boundary[0]
            d1 = boundary[1]
            d2 = boundary[2]

        elif isinstance(boundary[1], Edge):
            e = boundary[1]
            d1 = boundary[0]
            d2 = boundary[2]

        else:
            e = boundary[2]
            d1 = boundary[0]
            d2 = boundary[1]

        p = tangent_point_dd(d1, d2)
        ep1 = tangent_point_ld(e, d1)
        ep2 = tangent_point_ld(e, d2)
        s = find_perpendicular_point_l1_p1_to_l2(
            (d1.center, d2.center), p, (ep1, ep2))

        # Create new edges
        e1 = Edge(start=d1.center, end=s)
        e2 = Edge(start=d2.center, end=s)
        e3 = Edge(start=p, end=s)
        e4 = Edge(start=d1.center, end=p)
        e5 = Edge(start=p, end=d2.center)
        e6 = Edge(start=d1.center, end=ep1)
        e7 = Edge(start=d2.center, end=ep2)
        misc.plotting.plot_edges([e1, e2, e3, e4, e5, e6, e7])
        return [e1, e2, e3, e4, e5, e6, e7], 4

    else:  # len(edges) == 0
        d = find_inscribed_circle_to_Rp(boundary)
        c = d.center
        e1 = Edge(start=boundary[0].center, end=c)
        e2 = Edge(start=boundary[1].center, end=c)
        e3 = Edge(start=boundary[2].center, end=c)
        tangent_c0_c1 = tangent_point_dd(boundary[0], boundary[1])
        tangent_c1_c2 = tangent_point_dd(boundary[1], boundary[2])
        tangent_c2_c0 = tangent_point_dd(boundary[2], boundary[0])
        e4 = Edge(start=tangent_c0_c1, end=c)
        e5 = Edge(start=tangent_c1_c2, end=c)
        e6 = Edge(start=tangent_c2_c0, end=c)

        e7 = Edge(start=boundary[0].center, end=tangent_c0_c1)
        e8 = Edge(start=tangent_c0_c1, end=boundary[1].center)
        e9 = Edge(start=boundary[1].center, end=tangent_c1_c2)
        e10 = Edge(start=tangent_c1_c2, end=boundary[2].center)
        e11 = Edge(start=boundary[2].center, end=tangent_c2_c0)
        e12 = Edge(start=tangent_c2_c0, end=boundary[0].center)
        misc.plotting.plot_edges(
            [e1, e2, e3, e4, e5, e6, e7, e8, e9, e10, e11, e12])
        return [e1, e2, e3, e4, e5, e6, e7, e8, e9, e10, e11, e12], 6
