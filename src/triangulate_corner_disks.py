from disk import Disk
from edge import Edge
from circle_tools import tangent_point_ld, tangent_point_dd


def triangulate_corner_disk_convex(d: Disk, e1: Edge, e2: Edge):
    # Unpack edges
    t0 = Edge(start=e1.end, end=d.center)
    t1 = Edge(start=tangent_point_ld(e1, d), end=d.center)
    t2 = Edge(start=tangent_point_ld(e2, d), end=d.center)

    return [t0, t1, t2]


def triangulate_corner_disk_concave(d1: Disk, d2: Disk, e1: Edge, e2: Edge):
    # Unpack edges
    t0 = Edge(start=e1.end, end=d1.center)
    t1 = Edge(start=e1.end, end=d2.center)

    t2 = Edge(start=d1.center, end=d2.center)

    t3 = Edge(start=e1.end, end=tangent_point_dd(d1, d2))

    t4 = Edge(start=tangent_point_ld(e1, d1), end=d1.center)
    t5 = Edge(start=tangent_point_ld(e2, d2), end=d2.center)

    return [t0, t1, t2, t3, t4, t5]
