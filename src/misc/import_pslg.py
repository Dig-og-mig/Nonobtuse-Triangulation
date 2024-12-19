import json
from misc.instance import Cgshop2025Instance
from models.pslg import PSLG
from models.edge import Edge
from models.point import Point
import mpmath as mp


def import_pslg(path):
    f = open(path)
    js = json.load(f)
    f.close()
    instance = Cgshop2025Instance(**js)
    edges = []
    for i in range(len(instance.region_boundary)):
        edge = instance.region_boundary[i]
        x_i = instance.points_x[edge]
        y_i = instance.points_y[edge]
        if i == len(instance.region_boundary) - 1:
            x_j = instance.points_x[instance.region_boundary[0]]
            y_j = instance.points_y[instance.region_boundary[0]]
        else:
            x_j = instance.points_x[instance.region_boundary[i+1]]
            y_j = instance.points_y[instance.region_boundary[i+1]]
        edges.append(Edge(start=[x_i, y_i], end=[x_j, y_j]))
    temp_holes = []
    for i in range(len(instance.additional_constraints)):
        e = instance.additional_constraints[i]
        x1 = instance.points_x[e[0]]
        y1 = instance.points_y[e[0]]
        x2 = instance.points_x[e[1]]
        y2 = instance.points_y[e[1]]
        temp_holes.append(Edge(start=[x1, y1], end=[x2, y2]))
    holes = []

    def helper(e):
        for i in range(len(holes)):
            for e2 in holes[i]:
                if e.start == e2.start or e.end == e2.start or e.start == e2.end or e.end == e2.end:
                    return holes[i].append(e)

        return holes.append([e])

    for e in temp_holes:
        if holes == []:
            holes.append([e])
            continue
        helper(e)

    def helper2(i):
        x = mp.mpf(instance.points_x[i])
        y = mp.mpf(instance.points_y[i])
        for e in edges:
            if e.start == (x, y) or e.end == (x, y):
                return
        for h in holes:
            if isinstance(h[0], Point):
                continue
            for e in h:
                if e.start == (x, y) or e.end == (x, y):
                    return
        holes.append([Point(point=[x, y])])

    for i in range(instance.num_points):
        helper2(i)

    # Merge holes with boundary

    def helper3(i):
        for e2 in holes[i]:
            if isinstance(e2, Point):
                continue
            for e in edges:
                if isinstance(e, Point):
                    continue
                if e.start == e2.start or e.end == e2.start or e.start == e2.end or e.end == e2.end:
                    edges.extend(holes[i])
                    pops.append(i)
                    return
    pops = []
    for i in range(len(holes)):
        helper3(i)

    for i in sorted(pops, reverse=True):
        del holes[i]

    return PSLG(boundary=edges, holes=holes)
