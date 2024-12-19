import matplotlib.pyplot as plt
from models.edge import Edge
from models.disk import Disk
from models.point import Point


def handle_close(evt):
    raise SystemExit('Closed figure, exit program.')


plotting_delay = None
delay_amount = None
should_draw = None


def toggle_draw():
    global should_draw
    should_draw = not should_draw


def plotting_setup(xlim=(-20, 20), ylim=(-20, 20), delay=True, delay_amnt=0.000000000001, draw=True):
    global plotting_delay
    global delay_amount
    global should_draw
    plotting_delay = delay
    delay_amount = delay_amnt
    should_draw = draw
    if not plt.fignum_exists(1):
        plt.figure(figsize=(6, 6))
    plt.subplots_adjust(bottom=0.17)
    plt.xlim(xlim)
    plt.ylim(ylim)
    plt.gca().set_aspect('equal')
    if plt.rcParams['keymap.pan']:
        plt.rcParams['keymap.pan'].remove('p')
        plt.rcParams['keymap.back'].remove('left')
        plt.rcParams['keymap.back'].remove('c')
        plt.rcParams['keymap.forward'].remove('right')
        plt.rcParams['keymap.yscale'].remove('l')
    plt.connect('close_event', handle_close)


def toggle_delay():
    global plotting_delay
    plotting_delay = not plotting_delay
    print('Toggling delay: ', plotting_delay)


def is_delay():
    return plotting_delay


def clear_plot():
    x_lim = plt.gca().get_xlim()
    y_lim = plt.gca().get_ylim()
    plt.clf()
    plotting_setup(x_lim, y_lim, delay=plotting_delay)

# def plot_polygon(p : polygon.Polygon, delay=0.1):
#     for e in p.boundary:
#         plot_edge(e, delay = delay)


def plot_polygon(edges: list[Edge]):
    if edges:
        plot_point(edges[0].start)
    for e in edges:
        plot_edge(e)


def plot_edge(e: Edge, color='black'):
    if not should_draw:
        return
    plt.plot([e.start[0], e.end[0]], [e.start[1], e.end[1]],
             color=color, marker='o', markersize=0.1, linewidth=0.4)
    # plt.show(block=False)
    plt.draw()
    if plotting_delay:
        plt.pause(delay_amount)


def plot_edges(edges: list[Edge], color='black'):
    if not should_draw:
        return
    for edge in edges:
        plot_edge(edge, color)


def plot_line(p1, p2, color='black'):
    if not should_draw:
        return
    plt.plot([p1[0], p2[0]], [p1[1], p2[1]], color=color, markersize=0.5)
    # plt.show(block=False)
    plt.draw()
    if plotting_delay:
        plt.pause(delay_amount)


def plot_disk(c: Disk, color='black', linestyle='-'):
    if not should_draw:
        return
    circle = plt.Circle(c.center, c.radius, color=color,
                        fill=False, linewidth=0.4, linestyle=linestyle)
    plt.gca().add_patch(circle)
    plt.draw()
    if plotting_delay:
        plt.pause(delay_amount)


def plot_elements(list_of_elements, color='black', linestyle='-'):
    if not should_draw:
        return
    for element in list_of_elements:
        if isinstance(element, Disk):
            plot_disk(element, color, linestyle)
        elif isinstance(element, Edge):
            plot_edge(element, color)
        else:
            plot_point(element, color)


def plot_disks(disks: list[Disk], color='black', linestyle='-'):
    if not should_draw:
        return
    for disk in disks:
        plot_disk(disk, color, linestyle)


def plot_point(p, color='black'):
    if not should_draw:
        return
    if isinstance(p, Point):
        plt.scatter(p.point[0], p.point[1], color=color, s=0.1)
    else:
        plt.scatter(p[0], p[1], color=color, s=0.1)
    # plt.show(block=False)
    plt.draw()
    if plotting_delay:
        plt.pause(delay_amount)


def plot_pslg(edges: list[Edge]):
    if edges:
        if isinstance(edges[0], Edge):
            plot_point(edges[0].start)
        else:
            plot_point(edges[0].point)
    for e in edges:
        if isinstance(e, Edge):
            plot_edge(e)
        else:
            plot_point(e.point)


# plotting_setup()
# plot_circle((0, 0), 10)
# plot_circle((15, 0), 5)
# plot_circle((0, 15), 5)
# plot_point((5, 5))
# plot_point((5, 10))
# plot_line((5, 5), (5, 10))
# plot_line((5, 5), (10, 5))
# plot_line((5, 5), (0, 5))
