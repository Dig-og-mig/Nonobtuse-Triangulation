from .edge import Edge
from .disk import Disk
from pydantic import BaseModel, Field
from typing import Union
import math
import misc.plotting
import copy

from triangulation.triangulate_corner_disks import triangulate_corner_disk_convex, triangulate_corner_disk_concave
from disk_packing.circle_tools import tangent_point_dd, tangent_point_ld, angle_between_points, cross_product, dist_point_to_point, dist_point_to_line
from disk_packing.corner_disks import disk_at_convex_corner, disk_at_concave_corner
from disk_packing.find_tangent_solution import find_tangent_solution
from disk_packing.feasible_disk_finder import feasible_disk_finder
from triangulation.triangulate_3remainder_region import triangulate_3remainder_region
from triangulation.triangulate_4remainder_region import triangulate_4remainder_region
from triangulation.triangulation_tools import find_tangent_points_of_Rp

tolerance = 1e-9


class Polygon(BaseModel):
    boundary: list[Union[Edge, Disk]] = Field(..., description="A list of edges that form the boundary of the polygon.",
                                              min_items=3, default_factory=list)
    holes: list[list[Union[Edge, Disk]]] = Field(
        ..., description="A list of lists of edges that form the holes of the polygon.", default_factory=list)
    edges: list[Edge] = Field(
        ..., description="A list of edges that form the triangulation of the polygon.", default_factory=list)
    disks: list[Disk] = Field(
        default_factory=list, description="A list of disks used for the triangulation")

    n: int = 0

    n_corner: int = 0
    n_triangles: int = 0
    n_disks: int = 0
    n_remainder: int = 0
    n_regions: int = 0

    def plot_polygon(self):
        misc.plotting.plot_polygon(self.boundary)
        for hole in self.holes:
            misc.plotting.plot_polygon(hole)

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

    ### Triangulate corner disks ###

    def triangulate_corner_disk_convex(self, d: Disk, e1: Edge, e2: Edge):
        """
        Function to triangulate a convex corner disk

        Args:
            d: Disk at the corner
            e1: First edge at the corner
            e2: Second edge at the corner
        """
        edges = triangulate_corner_disk_convex(d, e1, e2)
        for edge in edges:
            self.append_edge_to_edges(edge)
            misc.plotting.plot_edge(edge)
        self.n_triangles += 2

    def triangulate_corner_disk_concave(self, d1: Disk, d2: Disk, e1: Edge, e2: Edge):
        """
        Function to triangulate a concave corner disk

        Args:
            d1: First disk at the corner
            d2: Second disk at the corner
            e1: First edge at the corner
            e2: Second edge at the corner
        """
        edges = triangulate_corner_disk_concave(d1, d2, e1, e2)
        for edge in edges:
            self.append_edge_to_edges(edge)
            misc.plotting.plot_edge(edge)
        self.n_triangles += 4

    def triangulate_3remainder_region(self, boundary: list):
        """
        Function to triangulate a 3 remainder region

        Args:
            boundary: List of elements that form the boundary of the region
        """
        edges, t = triangulate_3remainder_region(boundary)
        self.n_triangles += t
        self.n_remainder += 1
        self.n_regions += 1
        self.append_edges_to_edges(edges)

    def triangulate_4remainder_region(self, boundary: list, c_star=None):
        """
        Function to triangulate a 4 remainder region

        Args:
            boundary: List of elements that form the boundary of the region
            c_star: Center of the added circle in case of a disk with more than 180degrees arch
        """
        # Check if it is a 3-remainder
        a1 = boundary[0]
        a2 = boundary[1]
        a3 = boundary[2]
        a4 = boundary[3]

        if isinstance(a1, Disk) and isinstance(a3, Disk):
            if math.isclose(a1.radius + a3.radius, dist_point_to_point(a1.center, a3.center), abs_tol=tolerance):
                self.triangulate_3remainder_region([a1, a2, a3])
                self.triangulate_3remainder_region([a1, a3, a4])
                return
        if isinstance(a2, Disk) and isinstance(a4, Disk):
            if math.isclose(a2.radius + a4.radius, dist_point_to_point(a2.center, a4.center), abs_tol=tolerance):
                self.triangulate_3remainder_region([a1, a2, a4])
                self.triangulate_3remainder_region([a2, a3, a4])
                return
        if isinstance(a1, Disk) and isinstance(a3, Edge):
            if math.isclose(a1.radius, dist_point_to_line(a1.center, a3), abs_tol=tolerance):
                self.triangulate_3remainder_region([a1, a2, a3])
                self.triangulate_3remainder_region([a1, a3, a4])
                return
        if isinstance(a2, Disk) and isinstance(a4, Edge):
            if math.isclose(a2.radius, dist_point_to_line(a2.center, a4), abs_tol=tolerance):
                self.triangulate_3remainder_region([a1, a2, a4])
                self.triangulate_3remainder_region([a2, a3, a4])
                return
        if isinstance(a1, Edge) and isinstance(a3, Disk):
            if math.isclose(a3.radius, dist_point_to_line(a3.center, a1), abs_tol=tolerance):
                self.triangulate_3remainder_region([a1, a2, a3])
                self.triangulate_3remainder_region([a1, a3, a4])
                return
        if isinstance(a2, Edge) and isinstance(a4, Disk):
            if math.isclose(a4.radius, dist_point_to_line(a4.center, a2), abs_tol=tolerance):
                self.triangulate_3remainder_region([a1, a2, a4])
                self.triangulate_3remainder_region([a2, a3, a4])
                return

        edges, t, r = triangulate_4remainder_region(boundary, c_star)
        self.n_triangles += t
        self.n_remainder += 1
        self.n_regions += r
        self.append_edges_to_edges(edges)

    ###### Disk packing methods ######

    def disk_packing(self):
        """
        Function to pack disks in the polygon

        :return: Void
        """
        self.pack_disks_in_polygon(self.boundary)

    def pack_disks_in_polygon(self, boundary: list):
        """
        Function to pack disks in a subregion of the polygon

        :param boundary: List of elements that form the boundary of the subregion
        :return: Void
        """
        # Base case: 3 remainder region
        if len(boundary) == 3:
            hole_in_boundary = False
            for hole in self.holes:
                point = None
                if isinstance(hole[0], Disk):
                    point = hole[0].center
                else:
                    point = hole[0].start
                if self.point_in_polygon(point, boundary):
                    hole_in_boundary = True
                    break
            if not hole_in_boundary:
                return self.triangulate_3remainder_region(boundary)

        # Base case: 4 remainder region
        elif len(boundary) == 4:
            hole_in_boundary = False
            for hole in self.holes:
                point = None
                if isinstance(hole[0], Disk):
                    point = hole[0].center
                else:
                    point = hole[0].start
                if self.point_in_polygon(point, boundary):
                    hole_in_boundary = True
                    break
            if not hole_in_boundary:
                return self.triangulate_4remainder_region(boundary)

        # Recursive case
        else:
            # Find 3 non-consecutive elements of the boundary that gives a feasible tangent disk
            neighbor1 = boundary[-1]
            e1 = boundary[0]
            e2 = boundary[1]
            neighbor2 = boundary[2]

            e1_tup = (neighbor1, e1, e2)
            e2_tup = (e1, e2, neighbor2)

            (solution_disk, min_tangent_elements) = self.find_min_tangent_disk(
                e1_tup, e2_tup, boundary)

            # Safe case: At least 3 elements are non-consecutive
            if len(min_tangent_elements) > 1:
                self.append_disk_to_disks(solution_disk)
                misc.plotting.plot_disk(solution_disk, color='magenta')

                # Split the region
                regions = self.split_region(
                    solution_disk, boundary, e1, e2, min_tangent_elements)

                # Recursively pack the regions
                for region in regions:
                    self.pack_disks_in_polygon(region)

            # Non-possible case: No tangent elements
            elif len(min_tangent_elements) == 0:
                raise ValueError("No tangent elements found - general case")

            # Unsafe case: The single tangent element might be a neighbor
            else:
                # Case: Choose two non-neighbors; neighbor1 and e2
                if min_tangent_elements[0] == neighbor1:
                    temp_boundary = boundary.copy()
                    temp_boundary.remove(e1)
                    neighbor1_tup = (boundary[-2], neighbor1, e1)
                    (solution_disk, min_tangent_elements) = self.find_min_tangent_disk(
                        neighbor1_tup, e2_tup, temp_boundary, solution_disk)  # e1)

                    if (len(min_tangent_elements) == 0):
                        raise ValueError(
                            "No tangent elements found - Case: neighbor1 and e2")

                    # Split the region
                    regions = self.split_region(
                        solution_disk, temp_boundary, neighbor1, e2, min_tangent_elements, True, e1)

                # Case: Choose two non-neighbors; e1 and neighbor2
                elif min_tangent_elements[0] == neighbor2:
                    temp_boundary = boundary.copy()
                    temp_boundary.remove(e2)
                    neighbor2_tup = (e2, neighbor2, boundary[3])
                    (solution_disk, min_tangent_elements) = self.find_min_tangent_disk(
                        e1_tup, neighbor2_tup, temp_boundary, solution_disk)  # e2)

                    if (len(min_tangent_elements) == 0):
                        raise ValueError(
                            "No tangent elements found - Case: e1 and neighbor2")

                    # Split the region
                    regions = self.split_region(
                        solution_disk, temp_boundary, e1, neighbor2, min_tangent_elements, True, e2)

                # Case: The single tangent element is safe
                else:
                    regions = self.split_region(
                        solution_disk, boundary, e1, e2, min_tangent_elements, False, None)

                self.append_disk_to_disks(solution_disk)
                misc.plotting.plot_disk(solution_disk, color='magenta')

                # Recursively pack the regions
                for region in regions:
                    self.pack_disks_in_polygon(region)

    def find_min_tangent_disk(self, e1_tup, e2_tup, boundary: list, mid_elm=None):
        """
        Function to find a non-intersecting disk that is tangent to e1 and e2 and the boundary
         - If mid_elm is given, the elements e1 and e2 are non-consecutive

        Args:
            e1: First element
            e2: Second element
            boundary: List of elements that form the boundary of the region
            mid_elm: Middle element

        Returns:
            Tuple: Solution disk and the tangent elements
        """
        solution_disk = None
        solution_radius = math.inf
        tangent_elements = []
        solution_dist = math.inf

        e1 = e1_tup[1]
        e2 = e2_tup[1]

        for boundary_list in self.holes + [boundary]:
            for index in range(len(boundary_list)):
                elm = boundary_list[index]
                n1 = boundary_list[(index - 1) % len(boundary_list)]
                n2 = boundary_list[(index + 1) % len(boundary_list)]
                elm_tup = (n1, elm, n2)

                if elm == e1 or elm == e2:
                    continue
                else:
                    solutions = find_tangent_solution(
                        e1_tup, e2_tup, elm_tup)
                    if len(solutions) == 0:
                        continue

                    found_new_solution, solution_tuple = feasible_disk_finder(
                        solutions, (solution_disk, solution_radius, solution_dist), e1, e2, mid_elm)
                    if found_new_solution:
                        d = solution_tuple[0]
                        if solution_disk != None and math.isclose(d.radius, solution_radius, abs_tol=tolerance) and math.isclose(d.center[0], solution_disk.center[0], abs_tol=tolerance) and math.isclose(d.center[1], solution_disk.center[1], abs_tol=tolerance):
                            if not elm in tangent_elements:
                                tangent_elements.append(elm)
                                continue
                            else:
                                raise ValueError(
                                    "Element is represented twice in the boundary list - error.")
                        else:
                            solution_disk, solution_radius, solution_dist = solution_tuple
                            tangent_elements = [elm]

        return solution_disk, tangent_elements

    def split_region(self, disk: Disk, boundary: list, e1, e2, tangent_elements: list, lost_element=False, element=None):
        if len(tangent_elements) > 0:
            tangent_elements_ordered = tangent_elements
            tangent_elements_ordered += [e1, e2]
            # tangent_elements_ordered = tangent_elements + [e1, e2]

            def angle_for_element(disk, e1, elm):
                if isinstance(e1, Disk):
                    tangent_e1 = tangent_point_dd(e1, disk)
                else:
                    tangent_e1 = tangent_point_ld(e1, disk)

                if isinstance(elm, Disk):
                    tangent_elm = tangent_point_dd(elm, disk)
                else:
                    tangent_elm = tangent_point_ld(elm, disk)

                return angle_between_points(disk.center, tangent_e1, tangent_elm, False)

            tangent_elements_ordered = sorted(
                tangent_elements_ordered, key=lambda x: angle_for_element(disk, e1, x))

            # Split regions
            regions = []
            i = 0
            while i < len(tangent_elements_ordered):
                region = []

                if e1 == tangent_elements_ordered[i] and e2 == tangent_elements_ordered[(i+1) % len(tangent_elements_ordered)]:
                    region.append(e1.model_copy())
                    if lost_element:  # then append the lost element
                        region.append(element.model_copy())
                    region.append(e2.model_copy())
                    i = i + 1
                else:
                    def find_instance_index(lst, target_instance):
                        for i, item in enumerate(lst):
                            if item is target_instance:
                                return i
                        return -1  # If the instance isn't found

                    _from = find_instance_index(
                        boundary, tangent_elements_ordered[i])

                    # Collect elements in holes till next non-hole element
                    j = i + 1
                    holes = []
                    while j < len(tangent_elements_ordered):
                        # if the tangent_element is part of a hole
                        if (tangent_elements_ordered[j].hole):
                            # Find hole:
                            found_hole = None
                            for hole in self.holes:
                                if tangent_elements_ordered[j] in hole:
                                    found_hole = hole
                                    break
                            # if we already visisted the hole we will assume that this tangency point is non-existing.
                            if not found_hole == None:
                                self.holes.remove(found_hole)
                                # Index for element in hole:
                                tangent_element_index = find_instance_index(
                                    found_hole, tangent_elements_ordered[j])
                                # found_hole.index(tangent_elements_ordered[j])
                                hole_to_be_appended = []
                                if len(found_hole) == 1:
                                    elm_to_be_appended = found_hole[0].model_copy(
                                    )
                                    hole_to_be_appended.append(
                                        elm_to_be_appended)
                                else:  # +1 to include the tangent element twice
                                    for k in range(len(found_hole) + 1):
                                        found_hole[(tangent_element_index + k) %
                                                   len(found_hole)].hole = False
                                        elm_to_be_appended = found_hole[(
                                            tangent_element_index + k) % len(found_hole)].model_copy()
                                        hole_to_be_appended.append(
                                            elm_to_be_appended)

                                # holes = [disk.model_copy()] + hole_to_be_appended + holes
                                holes.append(hole_to_be_appended)
                            j = j + 1
                        else:
                            break
                    hole_elements = []
                    for i in range(len(holes)-1, -1, -1):
                        hole_elements.append(disk.model_copy())
                        for elm in holes[i]:
                            hole_elements.append(elm)
                    for elm in hole_elements:
                        elm.hole = False

                    # Next non-hole element
                    # _to = boundary.index(tangent_elements_ordered[(j) % len(tangent_elements_ordered)])
                    _to = find_instance_index(
                        boundary, tangent_elements_ordered[(j) % len(tangent_elements_ordered)])
                    i = j

                    if _from < _to:
                        for k in range(_from, _to + 1):
                            region.append(copy.deepcopy(boundary[k]))
                    else:
                        for k in range(_from, len(boundary)):
                            region.append(copy.deepcopy(boundary[k]))
                        for k in range(0, _to + 1):
                            region.append(copy.deepcopy(boundary[k]))

                    region += hole_elements

                region.append(disk.model_copy())
                regions.append(region)

            return regions

        else:
            raise ValueError(
                "No tangent elements found - something went wrong")

    def point_in_polygon(self, point, boundary):
        t_points = find_tangent_points_of_Rp(boundary)
        num_vertices = len(t_points)
        x, y = point
        inside = False

        x1, y1 = t_points[0]

        for i in range(1, num_vertices + 1):
            # v2 = boundary[i % num_vertices]
            # x2, y2 = v2.center[0], v2.center[1]

            x2, y2 = t_points[i % num_vertices]

            if y1 != y2:
                if min(y1, y2) < y <= max(y1, y2):
                    x_intersection = (y - y1) * (x2 - x1) / (y2 - y1) + x1
                    if x == x_intersection:
                        return False

                    if x < x_intersection:
                        inside = not inside

            x1, y1 = x2, y2

        return inside


###### CORNER DISK METHODS ######


    def add_corner_disk(self, edge: Edge, next_edge: Edge):
        # Unpack edges
        (x1, y1), (x2, y2) = (edge.start, edge.end)
        (x4, y4) = next_edge.end

        # list of disks to be added
        disks = []

        # Cross product
        crossproduct = cross_product((x2, y2), (x1, y1), (x4, y4))
        # v1x * v2y - v1y * v2x

        # all elements in polygon
        all_elements_in_polygon = self.boundary + \
            [item for sublist in self.holes for item in sublist]
        if crossproduct < 0 and not math.isclose(crossproduct, 0, abs_tol=tolerance):
            d = disk_at_convex_corner(
                edge, next_edge, all_elements_in_polygon)  # self.boundary)
            self.append_disk_to_disks(d)
            self.n_disks -= 1
            self.n_corner += 1

            disks.append(d)

            misc.plotting.plot_disk(d, color="green")
            self.triangulate_corner_disk_convex(
                d, edge, next_edge)
        else:
            d1, d2 = disk_at_concave_corner(
                edge, next_edge, all_elements_in_polygon)  # self.boundary)
            self.append_disk_to_disks(d1)
            self.append_disk_to_disks(d2)
            self.n_corner += 2
            self.n_disks -= 2
            disks.append(d1)
            disks.append(d2)

            misc.plotting.plot_disk(d1, color="green")
            misc.plotting.plot_disk(d2, color="green")
            self.triangulate_corner_disk_concave(
                d1, d2, edge, next_edge)
        return disks

    def add_corner_disks(self):

        new_boundary = []

        # On boundary
        iterator = iter(self.boundary)
        length = len(self.boundary)
        edge = next(iterator)
        for i in range(length):
            new_boundary.append(edge)
            if (i == length-1):
                next_edge = self.boundary[0]
            else:
                next_edge = next(iterator)

            disks = self.add_corner_disk(edge, next_edge)
            new_boundary += disks

            edge = next_edge

        # On holes
        new_holes = []
        for hole in self.holes:
            new_hole = []
            iterator = iter(hole)
            length = len(hole)
            edge = next(iterator)
            for i in range(length):
                new_hole.append(edge)
                if (i == length-1):
                    next_edge = hole[0]
                else:
                    next_edge = next(iterator)

                disks = self.add_corner_disk(edge, next_edge)
                new_hole += disks

                edge = next_edge
            # self.holes[self.holes.index(hole)] = new_hole
            new_holes.append(new_hole)

        self.boundary = new_boundary
        self.holes = new_holes
