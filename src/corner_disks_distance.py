from circle_tools import dist_point_to_point, projection_scalar_for_line, angle_between_points
import math
from edge import Edge

def dist_point_to_line_inside_corner(p: tuple[float, float], e: Edge, e1 : Edge, e2 : Edge):
    """
    Function to find the distance between a point p and an edge e *inside* a corner.
    
    Args:
        p (tuple[float, float]): The point
        e (Edge): The edge
        e1 (Edge): The first edge
        e2 (Edge): The second edge
    
    Returns:
        float: The distance between the point and the edge
    """
    # Unpack
    x, y = p
    (x1, y1), (x2, y2) = (e.start, e.end)

    # If point is an end point of the edge return the length of the edge
    if (x, y) == (x1, y1):
        return dist_point_to_point(p, (x2, y2))
    if (x, y) == (x2, y2):
        return dist_point_to_point(p, (x1, y1))
    
    t = projection_scalar_for_line((x, y), e)
    if t < 0.0:
        xt = x1
        yt = y1
    elif t > 1.0:
        xt = x2
        yt = y2
    else:
        xt = (1 - t) * x1 + t * x2
        yt = (1 - t) * y1 + t * y2

    # Is xt,yt inside the corner?
    angle_e2_e1 = angle_between_points(e1.end, e2.end, e1.start)
    angle_e2_p = angle_between_points(e1.end, e2.end, (xt, yt))
    if angle_e2_p <= angle_e2_e1:
        return dist_point_to_point(p, (xt, yt))
    else:
        return math.inf