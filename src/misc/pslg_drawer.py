from models.edge import Edge
from models.pslg import PSLG
from models.point import Point
from enum import Enum
from copy import deepcopy
from timeit import default_timer as timer
from datetime import timedelta
import misc.plotting
import matplotlib.pyplot as plt
import numpy as np
import traceback
import mpmath as mp
import math


class State(Enum):
    BOUNDARY = 1
    CLUSTER = 2
    TRIANGULATING = 3
    TRIANGULATED = 4
    FINISHED = 5


class PSLG_Drawer(object):
    n_prime = 0
    p: PSLG = None
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
    disks = []
    edges = []
    draw = None

    def __init__(self, triangulate=True, p=None, draw=True):
        plt.draw()
        self.add_corner_disks = triangulate
        self.draw = draw
        self.init_state()
        self.cid_button = plt.connect('button_press_event',
                                      lambda event: on_click(event=event, self=self))
        self.cid_key = plt.connect('key_press_event',
                                   lambda event: on_key_press(event=event, self=self))
        if p:
            self.p = p
            self.state = State.CLUSTER
            self.hole_i = len(self.p.holes)
            self.p.plot_pslg()
            self.p.n = len(self.p.boundary)
            for hole in self.p.holes:
                self.p.n += len(hole)
        else:
            self.p = PSLG()

        while True:
            self.set_fig_text()
            try:
                if self.state == State.FINISHED:
                    if not self.file_saved:
                        js = self.p.model_dump_json(indent=4)
                        f = open("../pslgs/last_instance.json", "w")
                        f.write(js)
                        f.close()
                        self.file_saved = True
                        print("Pslg saved to pslgs/last_instance.json")
                        if not self.draw:
                            misc.plotting.toggle_draw()

                    if self.add_corner_disks:
                        self.state = State.TRIANGULATING
                        self.set_fig_text()
                        start = timer()
                        self.p.triangulate()
                        end = timer()
                        print("Triangulation took: ",
                              timedelta(seconds=end-start))
                        self.state = State.TRIANGULATED
                        self.set_fig_text()

            except Exception as e:
                self.error_text.set_text(str(e))
                traceback.print_exc()

            plt.waitforbuttonpress()

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
        elif a_l1 == mp.inf:
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

    def intersects(self, e1):
        for e2 in self.p.boundary:
            p = self.find_intersection_l1_l2(e1, e2)
            if p is None:
                continue
            proj1 = self.projection_scalar_for_line(p, e1)
            proj2 = self.projection_scalar_for_line(p, e2)
            # if 0 < proj1 < 1 and 0 < proj2 < 1:
            #     return True
            if 0 < proj1 < 1 and 0 < proj2 < 1 and not math.isclose(proj1, 1, abs_tol=1e-9) and not math.isclose(proj2, 1, abs_tol=1e-9) and not math.isclose(proj1, 0, abs_tol=1e-9) and not math.isclose(proj2, 0, abs_tol=1e-9):
                return True

        for h in self.p.holes:
            for e2 in h:
                p = self.find_intersection_l1_l2(e1, e2)
                if p is None:
                    continue
                proj1 = self.projection_scalar_for_line(p, e1)
                proj2 = self.projection_scalar_for_line(p, e2)
                # if 0 < proj1 < 1 and 0 < proj2 < 1:
                #     return True
                if 0 < proj1 < 1 and 0 < proj2 < 1 and not math.isclose(proj1, 1, abs_tol=1e-9) and not math.isclose(proj2, 1, abs_tol=1e-9) and not math.isclose(proj1, 0, abs_tol=1e-9) and not math.isclose(proj2, 0, abs_tol=1e-9):
                    return True
        return False

    def intersection_split(self, edge: Edge, merge=False):
        intersection = False
        i = 0
        holes2 = deepcopy(self.p.holes)
        for h in self.p.holes:
            h2 = deepcopy(h)
            for e in h:
                if isinstance(e, Point):
                    continue
                if e == edge:
                    break
                p = self.find_intersection_l1_l2(edge, e)
                if p is None:
                    continue
                t = self.projection_scalar_for_line(p, e)
                t_1 = self.projection_scalar_for_line(p, edge)

                if (math.isclose(t, 0, abs_tol=1e-9) or math.isclose(t, 1, abs_tol=1e-9)) and (math.isclose(t_1, 0, abs_tol=1e-9) or math.isclose(t_1, 1, abs_tol=1e-9)):
                    continue
                elif (math.isclose(t, 0, abs_tol=1e-9) or math.isclose(t, 1, abs_tol=1e-9)) and 0 < t_1 < 1:
                    continue
                elif (math.isclose(t_1, 0, abs_tol=1e-9) or math.isclose(t_1, 1, abs_tol=1e-9)) and 0 < t < 1:
                    continue
                elif 0 < t < 1 and 0 < t_1 < 1:
                    new_es = self.split_edge(e, p)
                    h2.remove(e)
                    h2.extend(new_es)
                    self.p.n += 1
                    new_edges = self.split_edge(edge, p)
                    self.n_prime += 1
                    misc.plotting.plot_point(p)
                    intersection = True
                    holes2[i] = deepcopy(h2)
                    self.p.holes = deepcopy(holes2)
                    if h2 != self.p.holes[self.hole_i]:
                        for e in self.p.holes[i]:
                            if e not in self.p.holes[self.hole_i]:
                                self.p.holes[self.hole_i] = [
                                    e] + self.p.holes[self.hole_i]
                        self.p.holes.pop(i)
                        self.hole_i -= 1
                        self.i = len(self.p.holes[self.hole_i])+1
                    else:
                        self.i += 1
                    for e in new_edges:
                        self.intersection_split(e)
                    return True
            i += 1

        boundary2 = deepcopy(self.p.boundary)
        for e in self.p.boundary:
            if e == edge:
                break
            p = self.find_intersection_l1_l2(edge, e)
            if p is None:
                continue
            t = self.projection_scalar_for_line(p, e)
            t_1 = self.projection_scalar_for_line(p, edge)

            if (math.isclose(t, 0, abs_tol=1e-9) or math.isclose(t, 1, abs_tol=1e-9)) and (math.isclose(t_1, 0, abs_tol=1e-9) or math.isclose(t_1, 1, abs_tol=1e-9)):
                continue
            elif (math.isclose(t, 0, abs_tol=1e-9) or math.isclose(t, 1, abs_tol=1e-9)) and 0 < t_1 < 1:
                continue
            elif (math.isclose(t_1, 0, abs_tol=1e-9) or math.isclose(t_1, 1, abs_tol=1e-9)) and 0 < t < 1:
                continue
            elif 0 < t < 1 and 0 < t_1 < 1:
                new_boundary = self.split_edge(e, p)
                misc.plotting.plot_point(p)
                self.p.n += 1
                self.n_prime += 1
                boundary2.remove(e)
                boundary2.extend(new_boundary)
                new_edges = self.split_edge(edge, p)
                self.p.boundary = deepcopy(boundary2)
                for e2 in new_edges:
                    self.intersection_split(e2, True)
                return True

        if not intersection:
            self.p.holes[self.hole_i].append(edge)
            self.p.n += 1
            self.i += 1
            if merge:
                self.p.boundary.append(edge)

            misc.plotting.plot_edge(edge)
            plt.draw()
        return intersection

    def split_edge(self, edge: Edge, point):
        e1 = Edge(start=edge.start, end=point)
        e2 = Edge(start=point, end=edge.end)
        return [e1, e2]

    def area(self, boundary: list[Edge], last_edge: Edge):
        A = 0
        for e in boundary:
            A += e.start[0] * e.end[1] - e.end[0] * e.start[1]
        A += last_edge.start[0] * last_edge.end[1] - \
            last_edge.end[0] * last_edge.start[1]
        return (A / 2)

    def too_close_to_prev_point(self, point):
        if self.i == 0:
            return False
        if self.i == 1:
            prev_point = self.first_point
        elif self.state == State.CLUSTER:
            prev_point = Point(point=self.p.holes[self.hole_i][-1].end)
        else:
            prev_point = Point(point=(self.p.boundary[-1].end))
        return np.linalg.norm(np.array(prev_point.point) - np.array(point.point)) < 1e-3

    def init_state(self):
        self.p = PSLG()
        self.n_prime = 0
        self.file_saved = False
        self.first_point = None
        self.i = 0
        self.state = State.BOUNDARY
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
        if self.state != State.TRIANGULATED and self.state != State.TRIANGULATING:
            return
        print("Toggling lines")
        if self.edges:
            for edge in self.edges:
                plt.gca().add_line(edge)
            self.edges = []
        else:
            for edge in plt.gca().lines[self.p.n - self.n_prime:]:
                self.edges.append(edge)
                edge.remove()

    def set_fig_text(self):
        s = self.state.name.lower()
        if self.state == State.BOUNDARY or self.state == State.CLUSTER:
            self.state_text.set_text(
                "You are drawing a " + s + " \n"
                + "Drawing delay: " + str(misc.plotting.is_delay()) + "\n"
                + "n = " + str(self.p.n) + "\n"
            )
        elif self.state == State.TRIANGULATED:
            self.state_text.set_text(
                "You have " + s + " \n"
                + "Drawing delay: " + str(misc.plotting.is_delay()) + "\n"
                + "n = " + str(self.p.n) + "\n"
            )
        else:
            self.state_text.set_text(
                "You are " + s + " \n"
                + "Drawing delay: " + str(misc.plotting.is_delay()) + "\n"
                + "n = " + str(self.p.n) + "\n"
            )
        if self.i == 0 and self.state == State.CLUSTER:
            self.figure_text.set_text("Left click to place point in a new " + s + "\n"
                                      + "Right click to close and finish pslg \n"
                                      + "Press 'x' to toggle drawing delay \n"
                                      + "Press 'r' to reset \n")
        elif self.state == State.CLUSTER:
            self.figure_text.set_text("Left click to place point in the " + s + "\n"
                                      + "Right click to close and finish the " + s + "\n"
                                      + "Press space bar to finish the " + s + "\n"
                                      + "Press 'x' to toggle drawing delay \n"
                                      + "Press 'r' to reset")
        elif self.state == State.BOUNDARY:
            self.figure_text.set_text("Left click to place point in the " + s + "\n"
                                      + "Right click to close and finish the " + s + "\n"
                                      + "Press 'x' to toggle drawing delay \n"
                                      + "Press 'r' to reset")
        elif self.state == State.TRIANGULATED:
            self.figure_text.set_text("Press 'x' to toggle drawing delay \n"
                                      + "Press 'r' to reset")


def on_key_press(event, self):
    self.error_text.set_text("")
    if event.key == 'x':
        misc.plotting.toggle_delay()
    if self.state == State.TRIANGULATING:
        return
    if event.key == ' ':
        if self.state == State.CLUSTER:

            if self.i > 0:
                self.i = 0
                self.hole_i += 1
                self.first_point = None
                hole = self.p.holes[self.hole_i-1]
                for e in hole:
                    if e in self.p.boundary:
                        for e in hole:
                            if e not in self.p.boundary:
                                self.p.boundary.append(e)
                        self.p.holes.remove(hole)
                        self.hole_i -= 1
                        break
            else:
                self.state = State.FINISHED
                self.p.n = len(self.p.boundary)
                for hole in self.p.holes:
                    self.p.n += len(hole)

    if event.key == 'p':
        self.toggle_disks()
        return
    elif event.key == 'l':
        self.toggle_lines()
        return

    elif event.key == 'r':
        print("Resetting")
        misc.plotting.clear_plot()
        self.init_state()
        return

    elif event.key == 't':
        self.triangulate = not self.triangulate

    if self.state == State.TRIANGULATED or self.state == State.FINISHED:
        return
    return


def on_click(event, self):
    self.error_text.set_text("")
    if event.xdata is None or event.ydata is None:
        return
    if self.state == State.TRIANGULATING or self.state == State.TRIANGULATED or self.state == State.FINISHED:
        return
    # Right click
    if event.button == 3:
        if self.state == State.CLUSTER:
            if self.i > 0:
                if self.i <= 2:
                    self.error_text.set_text(
                        "Cycle must have at least 3 edges, try again")
                    return
                edge = Edge(
                    start=self.p.holes[self.hole_i][self.i-2].end, end=self.first_point.point, hole=True)
                self.intersection_split(edge)
                self.i = 0
                self.hole_i += 1
                self.first_point = None
                self.p.n += 1
                hole = self.p.holes[self.hole_i-1]
                for e in hole:
                    if e in self.p.boundary:
                        for e in hole:
                            if e not in self.p.boundary:
                                self.p.boundary.append(e)
                        self.p.holes.remove(hole)
                        self.hole_i -= 1
                        break
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
                start=self.p.boundary[self.i-2].end, end=self.first_point.point)
            self.p.boundary.append(edge)
            self.state = State.CLUSTER
            self.first_point = None
            self.i = 0
            self.p.n += 1
            misc.plotting.plot_edge(edge)

        return
    # Left click
    if event.button == 1:
        x, y = event.xdata, event.ydata
        if self.too_close_to_prev_point(Point(point=[x, y])):
            self.error_text.set_text(
                "Point is too close to previous point, try again")
            return
        if self.state == State.CLUSTER:
            if self.first_point is None:
                self.first_point = Point(point=[x, y])
                self.p.holes.append([self.first_point])
                misc.plotting.plot_point(self.first_point)
                plt.draw()
                self.i += 1
                self.p.n += 1
            else:
                if self.i == 1:
                    self.p.holes[self.hole_i].remove(self.first_point)
                    edge = Edge(start=self.first_point.point,
                                end=[x, y], hole=True)
                    self.p.n -= 1
                else:
                    edge = Edge(
                        start=self.p.holes[self.hole_i][self.i-2].end, end=[x, y], hole=True)

                _ = self.intersection_split(
                    edge)

        elif self.state == State.BOUNDARY:
            if self.first_point is None:
                self.first_point = Point(point=[x, y])
                misc.plotting.plot_point(self.first_point)
                plt.draw()
                self.i += 1
            else:
                edge = None
                if self.i == 1:
                    edge = Edge(start=self.first_point.point, end=[x, y])
                else:
                    edge = Edge(
                        start=self.p.boundary[self.i-2].end, end=[x, y])

                if self.intersects(edge):
                    self.error_text.set_text(
                        "Edge intersects with another edge, try again")
                    return

                self.p.boundary.append(edge)
                self.p.n += 1
                misc.plotting.plot_edge(edge)
                plt.draw()
                self.i += 1
