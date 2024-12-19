from edge import Edge
from disk import Disk
from point import Point
import math
from pydantic import BaseModel, Field
from typing import Union
import plotting
from polygon import Polygon
import copy

from circle_tools import angle_between_edges, dist_point_to_line, intersection_point_proj_dd, dist_point_to_point

tolerance = 1e-9


class PSLG(BaseModel):
    boundary: list[Edge] = Field(..., description="A list of edges that form the boundary of the polygon.",
                                 min_items=3, default_factory=list)
    holes: list[list[Union[Edge, Point]]] = Field(
        ..., description="A list of lists of edges that form the holes of the polygon.", default_factory=list)
    edges: list[Edge] = Field(
        ..., description="A list of edges that form the triangulation of the polygon.", default_factory=list)
    disks: list[Disk] = Field(
        default_factory=list, description="A list of disks used for the triangulation")

    n: int = 0

    n_corner_edge: int = 0
    n_triangles: int = 0
    n_disks: int = 0
    n_remainder: int = 0
    n_regions: int = 0

    def plot_boundary(self):
        plotting.plot_pslg(self.boundary)

    def plot_pslg(self):
        plotting.plot_pslg(self.boundary)
        for hole in self.holes:
            plotting.plot_pslg(hole)

    ### Appending edges and disks to the polygon ###
    # Edges
    def append_edge_to_edges(self, e: Edge):
        """
        Function to append a single edge to the polygon

        Args:
            e: Edge to be appended
        """
        self.edges.append(e)

    def append_edges_to_edges(self, edges: list[Edge]):
        """
        Function to append a list of edges to the polygon

        Args:
            edges: List of edges to be appended
        """
        for edge in edges:
            self.append_edge_to_edges(edge)

    # Disks
    def append_disk_to_disks(self, d: Disk):
        """
        Function to append a single disk to the polygon

        Args:
            d: Disk to be appended
        """
        self.disks.append(d)
        self.n_disks += 1

    def append_disks_to_disks(self, disks: list[Disk]):
        """
        Function to append a list of disks to the polygon

        Args:
            disks: List of disks to be appended
        """
        for disk in disks:
            self.append_disk_to_disks(disk)

    def compute_area_of_boundary(self, boundary: list[Edge]):
        """
        Function to compute the area of the boundary of the polygon

        Args:
            boundary: List of edges that form the boundary of the polygon
        """
        area = 0
        for edge in boundary:
            area += edge.start[0] * edge.end[1] - edge.end[0] * edge.start[1]
        return area / 2

    def triangulate(self):
        holes = len(self.holes)
        polygons_to_triangulate = self.make_polygons()

        n_triangles_temp = self.n_triangles
        n_disks_temp = self.n_disks
        n_remainder_temp = self.n_remainder
        n_regions_temp = self.n_regions
        N = 0

        self.n_corner_edge += self.n

        # polygons_to_triangulate[1].disk_packing()
        for polygon in polygons_to_triangulate:
            N += len(polygon.boundary) + \
                sum([len(hole) for hole in polygon.holes])
            polygon.disk_packing()
            n_triangles_temp += polygon.n_triangles
            n_disks_temp += polygon.n_disks
            n_remainder_temp += polygon.n_remainder
            n_regions_temp += polygon.n_regions
            print("Polygon:")
            print("   -Number of disks: ", polygon.n_disks)
            print("   -Number of triangles: ", polygon.n_triangles)
            print("   -Number of remainder regions: ", polygon.n_remainder)
        print("Total:")
        print("")
        print("Number of elements in the PSLG: ", self.n)
        print("Number of corner and edge disks: ", self.n_corner_edge)
        print("Number of arcs: ", N)
        print("Number of disks: ", n_disks_temp)
        print("Number of triangles: ", n_triangles_temp)
        print("Number of remainder regions: ", n_remainder_temp)
        print("Number of regions: ", n_regions_temp)
        print("Number of clusters: ", holes)

    def make_polygons(self):
        self.n = len(self.boundary) + sum([len(hole) for hole in self.holes])
        self.add_disks_to_edges()
        list_of_boundaries = self.make_sub_problems()
        positive_boundaries = []
        negative_boundaries = []
        for boundary in list_of_boundaries:
            if len(boundary) == 1 and isinstance(boundary[0], Point):
                A = -1
                negative_boundaries.append((boundary, A))
            else:
                A = self.compute_area_of_boundary(boundary)
                if A > 0:
                    positive_boundaries.append((boundary, A))
                else:
                    negative_boundaries.append((boundary, A))

        outer_boundary = min(negative_boundaries, key=lambda x: x[1])
        negative_boundaries.remove(outer_boundary)

        holes = []
        for hole in negative_boundaries:
            if len(hole[0]) == 1 and isinstance(hole[0][0], Point):
                d = hole[0][0].disk.model_copy()
                d.hole = True
                disk_hole = [d]
            else:
                disk_hole = []
                for edge in hole[0]:
                    d = edge.disk_start.model_copy()
                    d.hole = True
                    disk_hole.append(d)
                    for elm in edge.disks_start_end:
                        d2 = elm.model_copy()
                        d2.hole = True
                        disk_hole.append(d2)

                d1 = disk_hole[0]
                d2 = disk_hole[-1]
                if not math.isclose(dist_point_to_point(d1.center, d2.center), 0, abs_tol=1e-9):
                    ValueError("The hole is not closed")
            holes.append((disk_hole, hole[1]))

        polygons = []
        for boundary, A in positive_boundaries:
            disk_boundary = []
            for edge in boundary:
                disk_boundary.append(edge.disk_start.model_copy())
                disk_boundary += edge.disks_start_end.copy()

            d1 = disk_boundary[0]
            d2 = disk_boundary[-1]
            if not math.isclose(dist_point_to_point(d1.center, d2.center), 0, abs_tol=1e-9):
                ValueError("The boundary is not closed")

            holes_to_add = []
            for hole in holes:
                point = hole[0][0].center
                inside = self.point_in_polygon(point, disk_boundary)
                if not math.isclose(-hole[1], A) and inside:
                    holes_to_add.append(hole[0])
            polygon = Polygon(boundary=disk_boundary, holes=holes_to_add)
            polygons.append(polygon)

        return polygons

    def point_in_polygon(self, point, boundary):
        num_vertices = len(boundary)
        x, y = point
        inside = False

        v1 = boundary[0]
        x1, y1 = v1.center[0], v1.center[1]

        for i in range(1, num_vertices + 1):
            v2 = boundary[i % num_vertices]
            x2, y2 = v2.center[0], v2.center[1]

            # Check if the edge is not horizontal to avoid division by zero
            # This is safe to ignore as we do not have to consider horizontal intersection
            # with a horizontal line.
            if y1 != y2:
                # Determine if point is within the y range of the edge
                if min(y1, y2) < y <= max(y1, y2):
                    # Calculate x-intersection of the ray to the edge
                    x_intersection = (y - y1) * (x2 - x1) / (y2 - y1) + x1
                    # Check if point is exactly on the edge (for exclusion)
                    if x == x_intersection:
                        return False  # The point lies on the boundary

                    # Check if point is to the left of the x-intersection
                    if x < x_intersection:
                        # Flip the inside flag
                        inside = not inside

            v1, x1, y1 = v2, x2, y2

        return inside


# positive orientation remember this word

    def add_disks_to_edges(self):
        all_edges = self.boundary.copy()
        for hole in self.holes:
            all_edges += hole.copy()
        for edge in self.boundary:
            self.add_disks_to_edge(edge, all_edges)
            if isinstance(edge, Point):
                self.n_corner_edge += 1
            else:
                self.n_corner_edge += len(edge.disks_start_end)
        for hole in self.holes:
            for edge in hole:
                self.add_disks_to_edge(edge, all_edges)
                if isinstance(edge, Point):
                    self.n_corner_edge += 1
                else:
                    self.n_corner_edge += len(edge.disks_start_end)

    def add_disks_to_edge(self, edge: Edge, all_edges: list[Edge]):
        if isinstance(edge, Point):
            smallestDistance = math.inf
            corner_point = edge.point
            for elm in all_edges:
                if isinstance(elm, Point):
                    if elm.point == corner_point:
                        continue
                    dist = dist_point_to_point(corner_point, elm.point)
                else:
                    if elm.start == corner_point:
                        dist = dist_point_to_point(corner_point, elm.end)
                    elif elm.end == corner_point:
                        dist = dist_point_to_point(corner_point, elm.start)
                    else:
                        dist = dist_point_to_line(corner_point, elm)

                # If the distance is smaller - update the smallest distance
                if dist < smallestDistance:
                    smallestDistance = dist

            radius = smallestDistance * (9/20)
            center = corner_point
            d_start = Disk(center=center, radius=radius)
            plotting.plot_disk(d_start, color="green")
            edge.disk = d_start
        if isinstance(edge, Edge):
            if edge.disk_start == None:
                smallestDistance = math.inf
                corner_point = edge.start
                for elm in all_edges:

                    if isinstance(elm, Point):
                        dist = dist_point_to_point(corner_point, elm.point)
                    else:
                        if elm.start == corner_point:
                            dist = dist_point_to_point(corner_point, elm.end)
                        elif elm.end == corner_point:
                            dist = dist_point_to_point(corner_point, elm.start)
                        else:
                            dist = dist_point_to_line(corner_point, elm)

                    # If the distance is smaller - update the smallest distance
                    if dist < smallestDistance:
                        smallestDistance = dist

                radius = smallestDistance * (9/20)
                center = corner_point
                d_start = Disk(center=center, radius=radius)
                plotting.plot_disk(d_start, color="green")
                edge.disk_start = d_start
            if edge.disk_end == None:
                smallestDistance = math.inf
                corner_point = edge.end
                for elm in all_edges:
                    if isinstance(elm, Point):
                        dist = dist_point_to_point(corner_point, elm.point)

                    else:
                        if elm.start == corner_point:
                            dist = dist_point_to_point(corner_point, elm.end)
                        elif elm.end == corner_point:
                            dist = dist_point_to_point(corner_point, elm.start)
                        else:
                            dist = dist_point_to_line(corner_point, elm)

                    # If the distance is smaller - update the smallest distance
                    if dist < smallestDistance:
                        smallestDistance = dist

                radius = smallestDistance * (9/20)
                center = corner_point
                d_end = Disk(center=center, radius=radius)
                plotting.plot_disk(d_end, color="green")
                edge.disk_end = d_end

            # This is the dynamic start (starting from the intersection of the corner disk)
            start = intersection_point_proj_dd(edge.disk_start, edge.disk_end)
            # And this is the end (ending at the intersection of the corner disk)
            end = intersection_point_proj_dd(edge.disk_end, edge.disk_start)

            # Find smallest distance to edge
            smallestDistance = math.inf
            for elm in all_edges:
                if elm == edge:
                    continue
                if isinstance(elm, Point):
                    dist1 = dist_point_to_point(start, elm.point)
                    dist2 = dist_point_to_point(end, elm.point)
                    dist3 = dist_point_to_line(
                        elm.point, Edge(start=start, end=end))
                    min_dist = min(dist1, dist2, dist3)
                    if min_dist < smallestDistance:
                        smallestDistance = min_dist
                    continue

                else:
                    dist1 = dist_point_to_line(start, elm)
                    dist2 = dist_point_to_line(end, elm)
                    dist3 = dist_point_to_line(
                        elm.start, Edge(start=start, end=end))
                    dist4 = dist_point_to_line(
                        elm.end, Edge(start=start, end=end))

                    min_dist = min(dist1, dist2, dist3, dist4)
                    if min_dist < smallestDistance:
                        smallestDistance = min_dist

            radius = smallestDistance * (9/20)
            if radius > dist_point_to_point(start, end) / 2:
                radius = dist_point_to_point(start, end) / 2

            disk_list_start_end = []

            count = 0
            while dist_point_to_point(end, start) >= 2 * radius:
                count += 1
                helper_disk = Disk(center=start, radius=radius)
                center = intersection_point_proj_dd(helper_disk, edge.disk_end)
                d = Disk(center=center, radius=radius)
                plotting.plot_disk(d, color="blue")
                disk_list_start_end.append(d)

                start = intersection_point_proj_dd(d, edge.disk_end)

            if not math.isclose(dist_point_to_point(end, start), 0, abs_tol=1e-9):
                helper_disk = Disk(
                    center=start, radius=dist_point_to_point(end, start)/2)
                center = intersection_point_proj_dd(helper_disk, edge.disk_end)
                d = Disk(center=center, radius=dist_point_to_point(end, start)/2)
                plotting.plot_disk(d, color="blue")
                disk_list_start_end.append(d)

            edge.disks_start_end = disk_list_start_end
            edge.disks_end_start = disk_list_start_end.copy()[::-1]

    def make_sub_problems(self):
        sub_problems = []
        boundary = self.boundary.copy()
        holes = self.holes.copy()
        sub_problems += self.make_sub_problems_from_list(boundary)
        for hole in holes:
            sub_problems += self.make_sub_problems_from_list(hole)
        return sub_problems

    def make_sub_problems_from_list(self, boundary: list[Edge]):

        zero_times = boundary.copy()
        one_times = []
        two_times = []

        sub_problems = []

        count = 0

        while (zero_times or one_times):  # and count < 100:
            count += 1
            if zero_times:
                e = zero_times.pop(0)
                if isinstance(e, Point):
                    e.visited = 2
                    two_times.append(e)
                    sub_problems.append([e])
                    continue
                one_times.append(Edge(start=e.end, end=e.start, visited=1, disk_start=e.disk_end.model_copy(
                ), disk_end=e.disk_start.model_copy(), disks_start_end=e.disks_end_start.copy(), disks_end_start=e.disks_start_end.copy()))
                e.visited = 1

                if Edge(start=e.start, end=e.end, visited=0, disk_start=e.disk_start, disk_end=e.disk_end, disks_start_end=e.disks_start_end, disks_end_start=e.disks_end_start) in boundary:
                    index = boundary.index(Edge(start=e.start, end=e.end, visited=0, disk_start=e.disk_start,
                                           disk_end=e.disk_end, disks_start_end=e.disks_start_end, disks_end_start=e.disks_end_start))
                    boundary[index].visited = 1
                if Edge(start=e.end, end=e.start, visited=0, disk_start=e.disk_end, disk_end=e.disk_start, disks_start_end=e.disks_end_start, disks_end_start=e.disks_start_end) in boundary:
                    index = boundary.index(Edge(start=e.end, end=e.start, visited=0, disk_start=e.disk_end,
                                           disk_end=e.disk_start, disks_start_end=e.disks_end_start, disks_end_start=e.disks_start_end))
                    boundary[index].visited = 1

            else:
                e = one_times.pop(0)
                two_times.append(e)
                e.visited = 2

                if Edge(start=e.start, end=e.end, visited=1, disk_start=e.disk_start, disk_end=e.disk_end, disks_start_end=e.disks_start_end, disks_end_start=e.disks_end_start) in boundary:
                    index = boundary.index(Edge(start=e.start, end=e.end, visited=1, disk_start=e.disk_start,
                                           disk_end=e.disk_end, disks_start_end=e.disks_start_end, disks_end_start=e.disks_end_start))
                    boundary[index].visited = 2
                if Edge(start=e.end, end=e.start, visited=1, disk_start=e.disk_end, disk_end=e.disk_start, disks_start_end=e.disks_end_start, disks_end_start=e.disks_start_end) in boundary:
                    index = boundary.index(Edge(start=e.end, end=e.start, visited=1, disk_start=e.disk_end,
                                           disk_end=e.disk_start, disks_start_end=e.disks_end_start, disks_end_start=e.disks_start_end))
                    boundary[index].visited = 2

            edge_iterator = e.model_copy()
            first = True

            sub_problem = []
            while (edge_iterator != e or first):  # and count < 100:
                count += 1
                sub_problem.append(edge_iterator)

                first = False
                e_end = edge_iterator.end
                nexts = []

                for elm in boundary:
                    if elm.start == edge_iterator.start and elm.end == edge_iterator.end:
                        continue
                    if elm.start == edge_iterator.end and elm.end == edge_iterator.start:
                        continue

                    if e_end == elm.start:
                        nexts.append(elm)
                    if e_end == elm.end:
                        nexts.append(Edge(start=elm.end, end=elm.start, visited=elm.visited, disk_start=elm.disk_end.model_copy(
                        ), disk_end=elm.disk_start.model_copy(), disks_start_end=elm.disks_end_start.copy(), disks_end_start=elm.disks_start_end.copy()))

                if len(nexts) == 0:
                    edge_iterator = Edge(start=edge_iterator.end, end=edge_iterator.start, visited=1, disk_start=edge_iterator.disk_end.model_copy(
                    ), disk_end=edge_iterator.disk_start.model_copy(), disks_start_end=edge_iterator.disks_end_start.copy(), disks_end_start=edge_iterator.disks_start_end.copy())
                elif len(nexts) == 1:
                    edge_iterator = nexts[0].model_copy()
                else:  # len(nexts) > 1:
                    elm_with_smallest_angle = min(
                        nexts, key=lambda x: angle_between_edges(edge_iterator, x))
                    edge_iterator = elm_with_smallest_angle.model_copy()

                if edge_iterator.start == e.start and edge_iterator.end == e.end:
                    continue

                if edge_iterator.visited == 0:

                    if Edge(start=edge_iterator.start, end=edge_iterator.end, visited=0, disk_start=edge_iterator.disk_start, disk_end=edge_iterator.disk_end, disks_start_end=edge_iterator.disks_start_end, disks_end_start=edge_iterator.disks_end_start) in zero_times:
                        zero_times.remove(Edge(start=edge_iterator.start, end=edge_iterator.end, visited=0, disk_start=edge_iterator.disk_start.model_copy(
                        ), disk_end=edge_iterator.disk_end.model_copy(), disks_start_end=edge_iterator.disks_start_end.copy(), disks_end_start=edge_iterator.disks_end_start.copy()))
                    if Edge(start=edge_iterator.end, end=edge_iterator.start, visited=0, disk_start=edge_iterator.disk_start, disk_end=edge_iterator.disk_end, disks_start_end=edge_iterator.disks_start_end, disks_end_start=edge_iterator.disks_end_start) in zero_times:
                        zero_times.remove(Edge(start=edge_iterator.end, end=edge_iterator.start, visited=0, disk_start=edge_iterator.disk_start.model_copy(
                        ), disk_end=edge_iterator.disk_end.model_copy(), disks_start_end=edge_iterator.disks_start_end.copy(), disks_end_start=edge_iterator.disks_end_start.copy()))

                    if Edge(start=edge_iterator.start, end=edge_iterator.end, visited=0, disk_start=edge_iterator.disk_start, disk_end=edge_iterator.disk_end, disks_start_end=edge_iterator.disks_start_end, disks_end_start=edge_iterator.disks_end_start) in boundary:
                        index = boundary.index(Edge(start=edge_iterator.start, end=edge_iterator.end, visited=0, disk_start=edge_iterator.disk_start.model_copy(
                        ), disk_end=edge_iterator.disk_end.model_copy(), disks_start_end=edge_iterator.disks_start_end.copy(), disks_end_start=edge_iterator.disks_end_start.copy()))
                        boundary[index].visited = 1
                    if Edge(start=edge_iterator.end, end=edge_iterator.start, visited=0, disk_start=edge_iterator.disk_end, disk_end=edge_iterator.disk_start, disks_start_end=edge_iterator.disks_end_start, disks_end_start=edge_iterator.disks_start_end) in boundary:
                        index = boundary.index(Edge(start=edge_iterator.end, end=edge_iterator.start, visited=0, disk_start=edge_iterator.disk_end.model_copy(
                        ), disk_end=edge_iterator.disk_start.model_copy(), disks_start_end=edge_iterator.disks_end_start.copy(), disks_end_start=edge_iterator.disks_start_end.copy()))
                        boundary[index].visited = 1

                    edge_iterator.visited = 1
                    one_times.append(Edge(start=edge_iterator.end, end=edge_iterator.start, visited=1, disk_start=edge_iterator.disk_end.model_copy(
                    ), disk_end=edge_iterator.disk_start.model_copy(), disks_start_end=edge_iterator.disks_end_start.copy(), disks_end_start=edge_iterator.disks_start_end.copy()))
                elif edge_iterator.visited == 1:

                    if Edge(start=edge_iterator.start, end=edge_iterator.end, visited=1, disk_start=edge_iterator.disk_start, disk_end=edge_iterator.disk_end, disks_start_end=edge_iterator.disks_start_end, disks_end_start=edge_iterator.disks_end_start) in one_times:
                        one_times.remove(Edge(start=edge_iterator.start, end=edge_iterator.end, visited=1, disk_start=edge_iterator.disk_start.model_copy(
                        ), disk_end=edge_iterator.disk_end.model_copy(), disks_start_end=edge_iterator.disks_start_end.copy(), disks_end_start=edge_iterator.disks_end_start.copy()))
                    if Edge(start=edge_iterator.end, end=edge_iterator.start, visited=1, disk_start=edge_iterator.disk_end, disk_end=edge_iterator.disk_start, disks_start_end=edge_iterator.disks_end_start, disks_end_start=edge_iterator.disks_start_end) in one_times:
                        one_times.remove(Edge(start=edge_iterator.end, end=edge_iterator.start, visited=1, disk_start=edge_iterator.disk_end.model_copy(
                        ), disk_end=edge_iterator.disk_start.model_copy(), disks_start_end=edge_iterator.disks_end_start.copy(), disks_end_start=edge_iterator.disks_start_end.copy()))

                    if Edge(start=edge_iterator.start, end=edge_iterator.end, visited=1, disk_start=edge_iterator.disk_start, disk_end=edge_iterator.disk_end, disks_start_end=edge_iterator.disks_start_end, disks_end_start=edge_iterator.disks_end_start) in boundary:
                        index = boundary.index(Edge(start=edge_iterator.start, end=edge_iterator.end, visited=1, disk_start=edge_iterator.disk_start.model_copy(
                        ), disk_end=edge_iterator.disk_end.model_copy(), disks_start_end=edge_iterator.disks_start_end.copy(), disks_end_start=edge_iterator.disks_end_start.copy()))
                        boundary[index].visited = 2
                    if Edge(start=edge_iterator.end, end=edge_iterator.start, visited=1, disk_start=edge_iterator.disk_end, disk_end=edge_iterator.disk_start, disks_start_end=edge_iterator.disks_end_start, disks_end_start=edge_iterator.disks_start_end) in boundary:
                        index = boundary.index(Edge(start=edge_iterator.end, end=edge_iterator.start, visited=1, disk_start=edge_iterator.disk_end.model_copy(
                        ), disk_end=edge_iterator.disk_start.model_copy(), disks_start_end=edge_iterator.disks_end_start.copy(), disks_end_start=edge_iterator.disks_start_end.copy()))
                        boundary[index].visited = 2

                    edge_iterator.visited = 2
                    two_times.append(edge_iterator)

            sub_problems.append(copy.deepcopy(sub_problem))

        return sub_problems
