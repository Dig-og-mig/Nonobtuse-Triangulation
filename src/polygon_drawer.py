import matplotlib.pyplot as plt
import numpy as np
from polygon import Polygon
from edge import Edge
import plotting
from enum import Enum
import traceback
import math

from timeit import default_timer as timer
from datetime import timedelta


class State(Enum):
    BOUNDARY = 1
    HOLE = 2
    ADDING_CORNER_DISKS = 3
    ADDED_CORNER_DISKS = 4
    DISK_PACKING = 5
    TRIANGULATED = 6
    FINISHED = 7


class Polygon_Drawer(object):
    p: Polygon = None
    cid_button = None
    cid_key = None
    first_point = None
    i = 0
    state = State.BOUNDARY
    state_text = None
    figure_text = None
    error_text = None
    hole_i = 0
    delay = 0.00000001
    add_corner_disks = None
    disk_pack = None
    file_saved = False
    redos = []
    disks = []
    edges = []
    draw = None

    def __init__(self, add_corner_disks=True, disk_pack=True, polygon=None, draw=True):
        plt.draw()
        self.add_corner_disks = add_corner_disks
        self.disk_pack = disk_pack
        self.draw = draw
        self.init_state()
        self.cid_button = plt.connect('button_press_event',
                                      lambda event: on_click(event=event, self=self))
        self.cid_key = plt.connect('key_press_event',
                                   lambda event: on_key_press(event=event, self=self))
        if polygon:
            self.p = polygon
            self.state = State.HOLE
            self.hole_i = len(self.p.holes)
            self.p.plot_polygon()
            self.p.n = len(self.p.boundary)
            for hole in self.p.holes:
                self.p.n += len(hole)
        else:
            self.p = Polygon()

        while True:
            self.set_fig_text()
            try:
                if self.state == State.FINISHED:
                    if not self.file_saved:
                        js = self.p.model_dump_json(indent=4)
                        f = open("polygons/last_instance.json", "w")
                        f.write(js)
                        f.close()
                        self.file_saved = True
                        print("Polygon saved to polygons/last_instance.json")
                        if not self.draw:
                            plotting.toggle_draw()

                    if self.add_corner_disks:
                        self.state = State.ADDING_CORNER_DISKS
                        self.set_fig_text()
                        start = timer()
                        self.p.add_corner_disks()
                        end = timer()
                        print("Time to add corner disks: ",
                              timedelta(seconds=end-start))
                        self.state = State.ADDED_CORNER_DISKS
                        self.set_fig_text()

                if self.state == State.ADDED_CORNER_DISKS and self.disk_pack:
                    self.state = State.DISK_PACKING
                    self.set_fig_text()
                    start = timer()
                    self.p.disk_packing()
                    end = timer()
                    print("Time to pack disks: ", timedelta(seconds=end-start))
                    self.state = State.TRIANGULATED
                    self.set_fig_text()
                    print("Number of triangles: ", self.p.n_triangles)
                    print("Number of disks: ", self.p.n_disks)
                    print("Number of elements in polygon: ", self.p.n)
                    print("Number of holes: ", self.hole_i)
                    print("Number of corner disks: ", self.p.n_corner)
                    print("Number of elements in boundary (total): ",
                          self.p.n + self.p.n_corner)
                    print("Number of remainder regions: ", self.p.n_remainder)
                    print("Number of regions: ", self.p.n_regions)

            except Exception as e:
                self.error_text.set_text(str(e))
                traceback.print_exc()

            plt.waitforbuttonpress()

    def point_in_polygon(self, point):
        num_vertices = len(self.p.boundary)
        x, y = point
        inside = False

        v1 = self.p.boundary[0]
        x1, y1 = v1.start[0], v1.start[1]

        for i in range(1, num_vertices + 1):
            v2 = self.p.boundary[i % num_vertices]
            x2, y2 = v2.start[0], v2.start[1]

            # Check if the edge is not horizontal to avoid division by zero
            # This is safe to ignore as we do not have to consider horizontal intersection
            # with a horizontal line.
            if y1 != y2:
                # Determine if point is within the y range of the edge
                if min(y1, y2) < y <= max(y1, y2):
                    # Calculate x-intersection of the ray to the edge
                    x_intersection = (y - y1) * (x2 - x1) / (y2 - y1) + x1
                    # Check if point is to the left of the x-intersection
                    if x <= max(x1, x2) and (x1 == x2 or x <= x_intersection):
                        # Flip the inside flag
                        inside = not inside

            v1, x1, y1 = v2, x2, y2

        return inside

    def intersects(self, e1):
        for e2 in self.p.boundary:
            p = self.find_intersection_l1_l2(e1, e2)
            if p is None:
                continue
            proj1 = self.projection_scalar_for_line(p, e1)
            proj2 = self.projection_scalar_for_line(p, e2)
            if 0 < proj1 < 1 and 0 < proj2 < 1 and not math.isclose(proj1, 1, abs_tol=1e-9) and not math.isclose(proj2, 1, abs_tol=1e-9) and not math.isclose(proj1, 0, abs_tol=1e-9) and not math.isclose(proj2, 0, abs_tol=1e-9):
                return True

        for h in self.p.holes:
            for e2 in h:
                p = self.find_intersection_l1_l2(e1, e2)
                if p is None:
                    continue
                proj1 = self.projection_scalar_for_line(p, e1)
                proj2 = self.projection_scalar_for_line(p, e2)
                if 0 < proj1 < 1 and 0 < proj2 < 1 and not math.isclose(proj1, 1, abs_tol=1e-9) and not math.isclose(proj2, 1, abs_tol=1e-9) and not math.isclose(proj1, 0, abs_tol=1e-9) and not math.isclose(proj2, 0, abs_tol=1e-9):
                    return True
        return False

    def find_intersection_l1_l2(self, l1, l2):
        """
        Function to find the intersection point between the perpendicular line to l1 passing through p and l2.

        Args:
            l1 (tuple): The first line
            l2 (tuple): The second line

        Returns:
            tuple: The intersection point
        """
        # Unpack
        (x1, y1), (x2, y2) = l1.start, l1.end
        (x3, y3), (x4, y4) = l2.start, l2.end

        # Calculate slope and intercept for L1
        if (math.isclose(x2, x1, abs_tol=1e-9)):
            a_l1 = math.inf
        else:
            a_l1 = (y2 - y1) / (x2 - x1)
        b1 = y1 - a_l1 * x1
        # Calculate slope and intercept for L2
        if (math.isclose(x4, x3, abs_tol=1e-9)):
            a_l2 = math.inf
        else:
            a_l2 = (y4 - y3) / (x4 - x3)
        b2 = y3 - a_l2 * x3
        if math.isclose(a_l1, a_l2, abs_tol=1e-9):
            return None
        # Find intersection between the perpendicular line and L2
        if a_l2 == math.inf:
            x_intersect = x3
            y_intersect = a_l1 * x_intersect + b1
        elif a_l1 == math.inf:
            x_intersect = x1
            y_intersect = a_l2 * x_intersect + b2
        else:
            x_intersect = (b2 - b1) / (a_l1 - a_l2)
            y_intersect = a_l2 * x_intersect + b2

        return [x_intersect, y_intersect]

    def projection_scalar_for_line(self, p: tuple[float, float], e: Edge):
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

    def area(self, boundary: list[Edge], last_edge: Edge):
        A = 0
        for e in boundary:
            A += e.start[0] * e.end[1] - e.end[0] * e.start[1]
        A += last_edge.start[0] * last_edge.end[1] - \
            last_edge.end[0] * last_edge.start[1]
        return (A / 2)

    def is_ccw(self, boundary: list[Edge], last_edge: Edge):
        return self.area(boundary, last_edge) > 0

    def too_close_to_prev_point(self, point):
        if self.i == 0:
            return False
        if self.i == 1:
            prev_point = self.first_point
        elif self.state == State.HOLE:
            prev_point = self.p.holes[self.hole_i][-1].end
        else:
            prev_point = self.p.boundary[-1].end
        return np.linalg.norm(np.array(prev_point) - np.array(point)) < 1e-3

    def init_state(self):
        self.p = Polygon()
        self.file_saved = False
        self.first_point = None
        self.i = 0
        self.state = State.BOUNDARY
        self.redos = []
        self.undos = []
        self.disks = []
        self.edges = []
        self.hole_i = 0
        self.state_text = plt.figtext(
            0.5, 0.95, "", fontsize=8,
            horizontalalignment='center',
            verticalalignment='center')
        self.figure_text = plt.figtext(
            0.50, 0.065, "", fontsize=8,
            horizontalalignment='center',
            verticalalignment='center')
        self.error_text = plt.figtext(
            0.5, 0.90, "", fontsize=8,
            horizontalalignment='center',
            verticalalignment='center',
            color='red')
        self.set_fig_text()

    def undo(self):
        if self.state == State.HOLE:
            if self.i > 1:
                self.redos.append(self.p.holes[self.hole_i].pop())
                self.i -= 1
                plt.gca().lines[-1].remove()
                self.p.n -= 1
            elif self.i == 1:
                self.i -= 1
                self.redos.append(self.first_point)
                self.first_point = None
                plt.gca().collections[-1].remove()
                self.p.holes.pop()
            elif self.i == 0:
                if self.hole_i > 0:
                    self.hole_i -= 1
                    self.i = len(self.p.holes[self.hole_i])
                    self.redos.append(self.p.holes[self.hole_i].pop())
                    self.first_point = self.p.holes[self.hole_i][0].start
                else:
                    self.state = State.BOUNDARY
                    self.i = len(self.p.boundary)
                    edge = self.p.boundary.pop()
                    self.redos.append(edge)
                    self.first_point = self.p.boundary[0].start

                plt.gca().lines[-1].remove()
                self.p.n -= 1

        elif self.state == State.BOUNDARY:
            if self.i > 1:
                plt.gca().lines[-1].remove()
                self.p.n -= 1
                self.redos.append(self.p.boundary.pop())
                self.i -= 1
            elif self.i == 1:
                self.i -= 1
                self.redos.append(self.first_point)
                self.first_point = None
                plt.gca().collections[-1].remove()

    def redo(self):
        if not self.redos:
            return

        if self.i == 0:
            self.first_point = self.redos.pop()
            self.i += 1
            plotting.plot_point(self.first_point)
            if self.state == State.HOLE:
                self.p.holes.append([])
            return
        if self.state == State.HOLE:
            edge = self.redos.pop()
            self.p.holes[self.hole_i].append(edge)
            self.i += 1
            self.p.n += 1
            plotting.plot_edge(edge)
            if edge.end == self.first_point:
                self.i = 0
                self.hole_i += 1
                self.first_point = None

        elif self.state == State.BOUNDARY:
            edge = self.redos.pop()
            self.p.boundary.append(edge)
            self.i += 1
            self.p.n += 1
            plotting.plot_edge(edge)
            if edge.end == self.first_point:
                self.i = 0
                self.state = State.HOLE
                self.first_point = None

    def toggle_disks(self):
        if self.state != State.TRIANGULATED:
            return
        print("Toggling disks")
        if self.disks:
            for disk in self.disks:
                plt.gca().add_patch(disk)
            self.disks = []
        else:
            for disk in plt.gca().patches:
                self.disks.append(disk)
                disk.remove()

    def toggle_lines(self):
        if self.state != State.TRIANGULATED and self.state != State.ADDED_CORNER_DISKS:
            return
        print("Toggling lines")
        if self.edges:
            for edge in self.edges:
                plt.gca().add_line(edge)
            self.edges = []
        else:
            for edge in plt.gca().lines[self.p.n:]:
                self.edges.append(edge)
                edge.remove()

    def set_fig_text(self):
        s = self.state.name.lower()
        if self.state == State.BOUNDARY or self.state == State.HOLE:
            self.state_text.set_text(
                "You are drawing a " + s + " \n"
                + "Add corner disks: " + str(self.add_corner_disks) + "\n"
                + "Disk pack: " + str(self.disk_pack) + "\n"
                + "n = " + str(self.p.n) + "\n"
            )
        elif self.state == State.ADDED_CORNER_DISKS or self.state == State.TRIANGULATED:
            self.state_text.set_text(
                "You have " + s + " \n"
                + "Add corner disks: " + str(self.add_corner_disks) + "\n"
                + "Disk pack: " + str(self.disk_pack) + "\n"
                + "n = " + str(self.p.n) + "\n"
            )
        else:
            self.state_text.set_text(
                "You are " + s + " \n"
                + "Add corner disks: " + str(self.add_corner_disks) + "\n"
                + "Disk pack: " + str(self.disk_pack) + "\n"
                + "n = " + str(self.p.n) + "\n"
            )
        if self.i == 0 and self.state == State.HOLE:
            self.figure_text.set_text("Left click to place point in a new hole \n"
                                      + "Right click to finish polygon \n"
                                      + "Use left arrow key to undo last point \n"
                                      + "Press 'c' to toggle corner disks \n"
                                      + "Press 'd' to toggle disk packing \n"
                                      + "Press 'r' to reset")
        elif self.state == State.BOUNDARY or self.state == State.HOLE:
            self.figure_text.set_text("Left Click to place point in the " + s + "\n"
                                      + "Right click to finish the " + s + "\n"
                                      + "Use left arrow key to undo last point \n"
                                      + "Press 'c' to toggle corner disks \n"
                                      + "Press 'd' to toggle disk packing \n"
                                      + "Press 'r' to reset")
        elif self.state == State.TRIANGULATED:
            self.figure_text.set_text("Press 'p' to toggle disks \n"
                                      + "Press 'r' to reset")
        else:
            self.figure_text.set_text("Press 'c' to toggle corner disks \n"
                                      + "Press 'd' to toggle disk packing \n"
                                      + "Press 'r' to reset")


def on_key_press(event, self):
    self.error_text.set_text("")
    if event.key == 'x':
        plotting.toggle_delay()
    if self.state == State.ADDING_CORNER_DISKS or self.state == State.DISK_PACKING:
        return
    if event.key == 'p':
        self.toggle_disks()
        return
    elif event.key == 'l':
        self.toggle_lines()
        return

    elif event.key == 'r':
        print("Resetting")
        plotting.clear_plot()
        self.init_state()
        return

    elif event.key == 'd':
        self.disk_pack = not self.disk_pack
    elif event.key == 'c':
        self.add_corner_disks = not self.add_corner_disks

    if self.state == State.TRIANGULATED or self.state == State.FINISHED:
        return

    if event.key == 'left':
        print("Undoing")
        self.undo()
    elif event.key == 'right':
        # Redo
        print("Redoing")
        if len(self.redos) > 0:
            self.redo()

    return


def on_click(event, self):
    # Right click
    self.error_text.set_text("")

    if event.xdata is None or event.ydata is None:
        return
    if self.state == State.ADDING_CORNER_DISKS or self.state == State.TRIANGULATED or self.state == State.FINISHED:
        return
    if event.button == 3:
        if self.state == State.HOLE:

            if self.i > 0:
                if self.i <= 2:
                    self.error_text.set_text(
                        "Hole must have at least 3 edges, try again")
                    return
                edge = Edge(
                    start=self.p.holes[self.hole_i][self.i-2].end, end=self.first_point, hole=True)
                if not self.is_ccw(self.p.holes[self.hole_i], edge):
                    self.p.holes[self.hole_i].append(edge)
                    self.i = 0
                    self.hole_i += 1
                    self.first_point = None
                    self.redos = []
                    self.p.n += 1
                    plotting.plot_edge(edge)
                else:
                    self.error_text.set_text(
                        "Hole is not clockwise, try again")
                    return
            else:
                self.state = State.FINISHED
                self.p.n = len(self.p.boundary)
                for hole in self.p.holes:
                    self.p.n += len(hole)

        else:
            if self.i > 0:
                if self.i <= 2:
                    self.error_text.set_text(
                        "Boundary must have at least 3 edges, try again")
                    return
            edge = Edge(
                start=self.p.boundary[self.i-2].end, end=self.first_point)
            if self.is_ccw(self.p.boundary, edge):
                self.p.boundary.append(edge)
                self.state = State.HOLE
                self.first_point = None
                self.i = 0
                self.redos = []
                self.p.n += 1
                plotting.plot_edge(edge)
            else:
                self.error_text.set_text(
                    "Boundary is not counter clockwise, try again")
                return
        return
    # Left click
    if event.button == 1:
        x, y = event.xdata, event.ydata
        if self.too_close_to_prev_point([x, y]):
            self.error_text.set_text(
                "Point is too close to previous point, try again")
            return
        if self.state == State.HOLE:
            if not self.point_in_polygon([x, y]):
                self.error_text.set_text(
                    "Point is not inside the polygon, try again")
                return
            if self.first_point is None:
                self.p.holes.append([])
                self.first_point = [x, y]
                self.redos = []
                plotting.plot_point([x, y])
                plt.draw()
                self.i += 1
            else:
                if self.i == 1:
                    edge = Edge(start=self.first_point, end=[x, y], hole=True)
                else:
                    edge = Edge(
                        start=self.p.holes[self.hole_i][self.i-2].end, end=[x, y], hole=True)
                if self.intersects(edge):
                    self.error_text.set_text(
                        "Edge intersects with another edge, try again")
                    return
                self.p.holes[self.hole_i].append(edge)
                self.redos = []
                self.p.n += 1
                plotting.plot_edge(edge)
                plt.draw()
                self.i += 1

        elif self.state == State.BOUNDARY:
            if self.first_point is None:
                self.first_point = [x, y]
                self.redos = []
                plotting.plot_point([x, y])
                plt.draw()
                self.i += 1
            else:
                edge = None
                if self.i == 1:
                    edge = Edge(start=self.first_point, end=[x, y])
                else:
                    edge = Edge(
                        start=self.p.boundary[self.i-2].end, end=[x, y])
                if self.intersects(edge):
                    self.error_text.set_text(
                        "Edge intersects with another edge, try again")
                    return
                self.p.boundary.append(edge)
                self.redos = []
                self.p.n += 1
                plotting.plot_edge(edge)
                plt.draw()
                self.i += 1
