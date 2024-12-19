import math
import numpy as np
from models.disk import Disk
from models.edge import Edge
import misc.plotting
from disk_packing.circle_tools import projection_scalar_for_line
from mpmath import mp

tolerance = 1e-9


def circleTangentToCCC(c1: Disk, c2: Disk, c3: Disk, s1: int, s2: int, s3: int):
    # Unpack
    (x1, y1), r1 = (c1.center, c1.radius)
    (x2, y2), r2 = (c2.center, c2.radius)
    (x3, y3), r3 = (c3.center, c3.radius)

    x1 = mp.mpf(x1)
    y1 = mp.mpf(y1)
    r1 = mp.mpf(r1)
    x2 = mp.mpf(x2)
    y2 = mp.mpf(y2)
    r2 = mp.mpf(r2)
    x3 = mp.mpf(x3)
    y3 = mp.mpf(y3)
    r3 = mp.mpf(r3)

    # Check if the circles overlap - should never happen though
    dist1 = mp.sqrt((x1-x2)**2 + (y1-y2)**2)
    dist2 = mp.sqrt((x1-x3)**2 + (y1-y3)**2)
    dist3 = mp.sqrt((x2-x3)**2 + (y2-y3)**2)
    if (dist1 < r1 + r2 and not math.isclose(dist1, r1 + r2, abs_tol=tolerance)):
        misc.plotting.plot_disk(c1, color='red')
        misc.plotting.plot_disk(c2, color='red')
        print("c1: ", c1)
        print("c2: ", c2)
        raise ValueError(
            "circleTangentToCCC: Circles overlap1 with: ", dist1 - (r1 + r2))
    if (dist2 < r3 + r1 and not math.isclose(dist2, r3 + r1, abs_tol=tolerance)):
        misc.plotting.plot_disk(c1, color='red')
        misc.plotting.plot_disk(c3, color='red')
        raise ValueError(
            "circleTangentToCCC: Circles overlap2 with: ", dist2 - (r3 + r1))
    if (dist3 < r2 + r3 and not math.isclose(dist3, r2 + r3, abs_tol=tolerance)):
        misc.plotting.plot_disk(c2, color='red')
        misc.plotting.plot_disk(c3, color='red')
        raise ValueError(
            "circleTangentToCCC: Circles overlap3 with: ", dist3 - (r2 + r3))

    if (math.isclose(x1, x2, abs_tol=tolerance) and math.isclose(y1, y2, abs_tol=tolerance)):
        return []
    if (math.isclose(x1, x3, abs_tol=tolerance) and math.isclose(y1, y3, abs_tol=tolerance)):
        return []
    if (math.isclose(x2, x3, abs_tol=tolerance) and math.isclose(y2, y3, abs_tol=tolerance)):
        return []

    M11 = 2*x2-2*x1
    M12 = 2*y2-2*y1
    M13 = 2*s1*r1-2*s2*r2
    b1 = x2**2 + y2**2 - (x1**2 + y1**2) + r1**2 - r2**2

    M21 = 2*x3-2*x1
    M22 = 2*y3-2*y1
    M23 = 2*s1*r1-2*s3*r3
    b2 = x3**2 + y3**2 - (x1**2 + y1**2) + r1**2 - r3**2

    if M11 == 0:  # collinear circles - c1 c2 - vertically
        if M21 == 0:  # collinear circles - c1 c3 - vertically
            return circleTangentToCCCvertical(c1, c2, c3, s1, s2, s3)
    elif M22 == 0:  # collinear circles - c1 c2 - horizontally
        if M12 == 0:  # collinear circles - c1 c3 - horizontally
            return circleTangentToCCChorizontal(c1, c2, c3, s1, s2, s3)

    # Calculate the area of the triangle formed by the centers
    area = x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2)

    # If the area is 0, the points are collinear
    if math.isclose(area, 0, abs_tol=tolerance):
        return circleTangentToCCCcolinear(c1, c2, c3, s1, s2, s3)

    M13p = (M12*M23-M13*M22)/(M22*M11-M21*M12)
    Mb1p = - (M12*b2-M22*b1)/(M22*M11-M21*M12)
    M23p = - (M23*M11-M21*M13)/(M22*M11-M21*M12)
    Mb2p = (b2*M11-M21*b1)/(M22*M11-M21*M12)

    # Second degree equation
    a = (M13p**2 + M23p**2 - 1)
    b = 2*Mb1p*M13p + 2*Mb2p*M23p - 2*M13p*x1 - 2*M23p*y1 + 2*s1*r1
    c = Mb1p**2 + Mb2p**2 + x1**2 + y1**2 - 2*Mb1p*x1 - 2*Mb2p*y1 - s1**2*r1**2

    D = b**2 - 4*a*c
    if math.isclose(D, 0, abs_tol=tolerance):
        D = 0

    if (D < 0):  # and (not math.isclose(D, 0, abs_tol=1e-3))):
        return []
    elif (a == 0):
        rs = np.array([-c / b])
        length = 1
    else:
        rs = np.array([(-b + mp.sqrt(D))/(2*a), (-b - mp.sqrt(D))/(2*a)])
        length = rs.size

    xs = Mb1p + M13p*rs  # should be array
    ys = Mb2p + M23p*rs

    if (len(xs) > 0):
        # return [Disk(center=(xs, ys), radius = abs(rs))] #
        return list(filter(lambda x: x.radius > 0, [Disk(center=(xs[i], ys[i]), radius=rs[i]) for i in range(length)]))
    else:
        raise ValueError("circleTangentToCCC: No solutions found")


def circleTangentToCCCcolinear(c1: Disk, c2: Disk, c3: Disk, s1: int, s2: int, s3: int):
    # Unpack
    (x1, y1), r1 = (c1.center, c1.radius)
    (x2, y2), r2 = (c2.center, c2.radius)
    (x3, y3), r3 = (c3.center, c3.radius)

    (x1p, y1p) = (x1-x1, y1-y1)
    (x2p, y2p) = (x2-x1, y2-y1)
    (x3p, y3p) = (x3-x1, y3-y1)

    # rotate to the x-axis
    angle = mp.atan2(y2p, x2p)  # should be same for c3
    x1pp = x1p  # should be 0
    y1pp = y1p  # should be 0

    x2pp = x2p*mp.cos(-angle) - y2p*mp.sin(-angle)
    y2pp = x2p*mp.sin(-angle) + y2p*mp.cos(-angle)

    x3pp = x3p*mp.cos(-angle) - y3p*mp.sin(-angle)
    y3pp = x3p*mp.sin(-angle) + y3p*mp.cos(-angle)

    im_solutions = circleTangentToCCChorizontal(Disk(center=(x1pp, y1pp), radius=r1), Disk(
        center=(x2pp, y2pp), radius=r2), Disk(center=(x3pp, y3pp), radius=r3), s1, s2, s3)

    solutions = []
    for sol in im_solutions:
        (x, y), r = (sol.center, sol.radius)
        xs = x*mp.cos(angle) - y*mp.sin(angle)
        ys = x*mp.sin(angle) + y*mp.cos(angle)
        xs = xs + x1
        ys = ys + y1
        solutions.append(Disk(center=(xs, ys), radius=r))

    return solutions


def circleTangentToCCCvertical(c1: Disk, c2: Disk, c3: Disk, s1: int, s2: int, s3: int):
    # Unpack
    (x1, y1), r1 = (c1.center, c1.radius)
    (x2, y2), r2 = (c2.center, c2.radius)
    (x3, y3), r3 = (c3.center, c3.radius)

    if r1 == r2 and r2 == r3:
        return []

    # M11 = 2*x2-2*x1
    M12 = 2*y2-2*y1
    M13 = 2*s1*r1-2*s2*r2
    b1 = x2**2 + y2**2 - (x1**2 + y1**2) + r1**2 - r2**2

    # M21 = 2*x3-2*x1
    M22 = 2*y3-2*y1
    M23 = 2*s1*r1-2*s3*r3
    b2 = x3**2 + y3**2 - (x1**2 + y1**2) + r1**2 - r3**2

    Mb1p = - (M13*b2-M23*b1)/(M23*M12-M22*M13)
    Mb2p = (b2*M12-M22*b1)/(M23*M12-M22*M13)

    ys = Mb1p
    rs = Mb2p

    if (math.isclose((rs-s1*r1)**2, (ys-y1)**2, abs_tol=tolerance)):
        xs1 = x1
        xs2 = x1
    else:
        xs1 = mp.sqrt((rs - s1*r1)**2 - (ys - y1)**2) + x1
        xs2 = -mp.sqrt((rs - s1*r1)**2 - (ys - y1)**2) + x1

    xs = np.array([xs1, xs2])
    ys = np.array([ys, ys])
    rs = np.array([rs, rs])

    # [Disk(center=(xs[i], ys[i]), radius=rs[i]) for i in range(2)]
    return list(filter(lambda x: x.radius > 0, [Disk(center=(xs[i], ys[i]), radius=rs[i]) for i in range(2)]))


def circleTangentToCCChorizontal(c1: Disk, c2: Disk, c3: Disk, s1: int, s2: int, s3: int):
    # Unpack
    (x1, y1), r1 = (c1.center, c1.radius)
    (x2, y2), r2 = (c2.center, c2.radius)
    (x3, y3), r3 = (c3.center, c3.radius)

    x1 = mp.mpf(x1)
    y1 = mp.mpf(y1)
    r1 = mp.mpf(r1)
    x2 = mp.mpf(x2)
    y2 = mp.mpf(y2)
    r2 = mp.mpf(r2)
    x3 = mp.mpf(x3)
    y3 = mp.mpf(y3)
    r3 = mp.mpf(r3)

    if r1 == r2 and r2 == r3:
        return []

    M11 = 2*x2-2*x1
    # M12 = 2*y2-2*y1
    M13 = 2*s1*r1-2*s2*r2
    b1 = x2**2 + y2**2 - (x1**2 + y1**2) + r1**2 - r2**2

    M21 = 2*x3-2*x1
    # M22 = 2*y3-2*y1
    M23 = 2*s1*r1-2*s3*r3
    b2 = x3**2 + y3**2 - (x1**2 + y1**2) + r1**2 - r3**2

    if (M23*M11-M21*M13 == 0.0):
        div = mp.mpf(10) ** -mp.dps
    else:
        div = (M23*M11-M21*M13)

    Mb1p = - (M13*b2-M23*b1)/div
    Mb2p = (b2*M11-M21*b1)/div

    xs = Mb1p
    rs = Mb2p

    if (math.isclose((rs-s1*r1)**2, (xs-x1)**2, abs_tol=tolerance)):
        ys1 = y1
        ys2 = y1
    else:
        ys1 = mp.sqrt((rs - s1*r1)**2 - (xs - x1)**2) + y1
        ys2 = -mp.sqrt((rs - s1*r1)**2 - (xs - x1)**2) + y1

    xs = np.array([xs, xs])
    ys = np.array([ys1, ys2])
    rs = np.array([rs, rs])

    return list(filter(lambda x: x.radius > 0, [Disk(center=(xs[i], ys[i]), radius=rs[i]) for i in range(2)]))


def circleTangentToCCL(c1: Disk, c2: Disk, L1: Edge, s1: int, s2: int):
    # Unpack
    (x1, y1), r1 = (c1.center, c1.radius)
    (x2, y2), r2 = (c2.center, c2.radius)
    (x3, y3), (x4, y4) = (L1.start, L1.end)

    if (math.isclose(x1, x2) and math.isclose(y1, y2)):
        return []
    if (math.isclose(x1, x3) and math.isclose(y1, y3)):
        return []
    if (math.isclose(x2, x3) and math.isclose(y2, y3)):
        return []

    M11 = 2*x2-2*x1
    M12 = 2*y2-2*y1
    M13 = 2*s1*r1-2*s2*r2
    b1 = x2**2 + y2**2 - (x1**2 + y1**2) + r1**2 - r2**2

    # Equation for a line
    A = y3 - y4
    B = x4 - x3
    C = (y4 - y3)*x3 - (x4 - x3)*y3
    # Denominator
    K = mp.sqrt((y4-y3)**2 + (x4-x3)**2)

    if A == 0:  # horizontal line
        if M11 == 0:  # collinear circles vertically
            return circleTangentToCCLvertical(c1, c2, L1, s1, s2)
    elif B == 0:  # vertical line
        if M12 == 0:  # collinear circles horizontally
            return circleTangentToCCLhorizontal(c1, c2, L1, s1, s2)

    projection_onto_line_c1 = (
        (x1 - x3) * (x4 - x3) + (y1 - y3) * (y4 - y3)) / ((x4 - x3)**2 + (y4 - y3)**2)
    projection_onto_line_c2 = (
        (x2 - x3) * (x4 - x3) + (y2 - y3) * (y4 - y3)) / ((x4 - x3)**2 + (y4 - y3)**2)

    if math.isclose(projection_onto_line_c1, projection_onto_line_c2, rel_tol=1e-9):
        return circleTangentToCCLcolProj(c1, c2, L1, s1, s2)

    # L1 positive
    M21 = A
    M22 = B
    M23 = -K
    b2 = -C

    res1 = circleTangentToCCL_helper(
        c1, s1, M11, M12, M13, b1, M21, M22, M23, b2)
    res = list(filter(lambda tup: (
        A*tup.center[0] + B*tup.center[1] + C)*tup.radius > 0, res1))

    # L1 negative
    M21 = -A
    M22 = -B
    M23 = -K
    b2 = C

    res2 = circleTangentToCCL_helper(
        c1, s1, M11, M12, M13, b1, M21, M22, M23, b2)
    res += list(filter(lambda tup: (-A *
                tup.center[0] - B*tup.center[1] - C)*tup.radius > 0, res2))

    return res


def circleTangentToCCL_helper(c1, s1, M11, M12, M13, b1, M21, M22, M23, b2):
    # Unpack
    (x1, y1), r1 = (c1.center, c1.radius)

    M13p = (M12*M23-M13*M22)/(M22*M11-M21*M12)
    Mb1p = - (M12*b2-M22*b1)/(M22*M11-M21*M12)
    M23p = - (M23*M11-M21*M13)/(M22*M11-M21*M12)
    Mb2p = (b2*M11-M21*b1)/(M22*M11-M21*M12)

    # Second degree equation
    a = (M13p**2 + M23p**2 - 1)
    b = 2*Mb1p*M13p + 2*Mb2p*M23p - 2*M13p*x1 - 2*M23p*y1 + 2*s1*r1
    c = Mb1p**2 + Mb2p**2 + x1**2 + y1**2 - 2*Mb1p*x1 - 2*Mb2p*y1 - s1**2*r1**2

    D = b**2 - 4*a*c
    if math.isclose(D, 0, abs_tol=tolerance):
        D = 0
    if (D < 0):
        return []
    elif (math.isclose(a, 0, abs_tol=tolerance)):
        if (math.isclose(b, 0, abs_tol=tolerance)):
            return []
        rs = np.array([-c / b])
        length = 1
    else:
        rs = np.array([(-b + mp.sqrt(D))/(2*a), (-b - mp.sqrt(D))/(2*a)])
        length = rs.size

    xs = Mb1p + M13p*rs  # should be array
    ys = Mb2p + M23p*rs

    return list(filter(lambda x: x.radius > 0, [Disk(center=(xs[i], ys[i]), radius=rs[i]) for i in range(length)]))


# CHANGE THIS WITH ROTATING TO THE X-AXIS AND BACK AND USE THE HORIZONTAL FUNCTION circleTangentToCCLhorizontal
def circleTangentToCCLcolProj(c1, c2, L1, s1, s2):
    # Unpack
    (x1, y1), r1 = (c1.center, c1.radius)
    (x2, y2), r2 = (c2.center, c2.radius)
    (x3, y3), (x4, y4) = (L1.start, L1.end)

    (x1p, y1p) = (x1-x1, y1-y1)
    (x2p, y2p) = (x2-x1, y2-y1)
    (x3p, y3p) = (x3-x1, y3-y1)
    (x4p, y4p) = (x4-x1, y4-y1)

    # rotate to the x-axis
    angle = mp.atan2(y2p, x2p)  # should be same for L1

    x1pp = x1p  # should be 0
    y1pp = y1p  # should be 0

    x2pp = x2p*mp.cos(-angle) - y2p*mp.sin(-angle)
    y2pp = x2p*mp.sin(-angle) + y2p*mp.cos(-angle)

    x3pp = x3p*mp.cos(-angle) - y3p*mp.sin(-angle)
    y3pp = x3p*mp.sin(-angle) + y3p*mp.cos(-angle)
    x4pp = x4p*mp.cos(-angle) - y4p*mp.sin(-angle)
    y4pp = x4p*mp.sin(-angle) + y4p*mp.cos(-angle)

    im_solutions = circleTangentToCCLhorizontal(
        Disk(center=(x1pp, y1pp), radius=r1),
        Disk(center=(x2pp, y2pp), radius=r2),
        Edge(start=(x3pp, y3pp), end=(x4pp, y4pp)), s1, s2)

    solutions = []
    for sol in im_solutions:
        (x, y), r = (sol.center, sol.radius)
        xs = x*mp.cos(angle) - y*mp.sin(angle)
        ys = x*mp.sin(angle) + y*mp.cos(angle)
        xs = xs + x1
        ys = ys + y1
        solutions.append(Disk(center=(xs, ys), radius=r))

    return solutions


def circleTangentToCCLvertical(c1, c2, L1, s1, s2):
    # Unpack
    (x1, y1), r1 = (c1.center, c1.radius)
    (x2, y2), r2 = (c2.center, c2.radius)
    (_, y3), (_, _) = (L1.start, L1.end)

    if (y1 > y3 > y2) or (y1 < y3 < y2):
        return []  # Then the line is between the two circles

    # M11 = 2*x2-2*x1
    M12 = 2*y2-2*y1
    M13 = 2*s1*r1-2*s2*r2
    b1 = x2**2 + y2**2 - (x1**2 + y1**2) + r1**2 - r2**2

    if math.isclose(M12, -M13):
        rs = (y3 - b1/M12)/2
    else:
        rs = (b1 - M12*y3)/(M12 + M13)
    rs = np.array([rs, rs])

    ys = np.array([y3 + rs[0], y3 - rs[0]])

    # if math.isclose((rs - s1*r1)**2 - (ys - y1)**2, 0, abs_tol=tolerance):
    #     xs1 = x1
    #     xs2 = x1
    # else:
    #     xs1 = np.sqrt((rs - s1*r1)**2 - (ys - y1)**2) + x1
    #     xs2 = -np.sqrt((rs - s1*r1)**2 - (ys - y1)**2) + x1

    solutions = []

    for i in range(ys.shape[0]):
        if math.isclose((rs[i] - s1*r1)**2 - (ys[i] - y1)**2, 0, abs_tol=tolerance):
            solutions.append(Disk(center=(x1, ys[i]), radius=rs))
        elif (rs[i] - s1*r1)**2 - (ys[i] - y1)**2 < 0:
            continue
        else:
            xs1 = mp.sqrt((rs[i] - s1*r1)**2 - (ys[i] - y1)**2) + x1
            xs2 = -mp.sqrt((rs[i] - s1*r1)**2 - (ys[i] - y1)**2) + x1
            solutions.append(Disk(center=(xs1, ys[i]), radius=rs[i]))
            solutions.append(Disk(center=(xs2, ys[i]), radius=rs[i]))

    return list(filter(lambda x: x.radius > 0, solutions))


def circleTangentToCCLhorizontal(c1, c2, L1, s1, s2):
    # Unpack
    (x1, y1), r1 = (c1.center, c1.radius)
    (x2, y2), r2 = (c2.center, c2.radius)
    (x3, _), (_, _) = (L1.start, L1.end)

    if (x1 > x3 > x2) or (x1 < x3 < x2):
        return []  # Then the line is between the two circles

    M11 = 2*x2-2*x1
    M13 = 2*s1*r1-2*s2*r2
    b1 = x2**2 + y2**2 - (x1**2 + y1**2) + r1**2 - r2**2

    if math.isclose(M11, -M13, abs_tol=tolerance):
        rs = (x3 - b1/M11)/2
    else:
        rs = (b1 - M11*x3)/(M11 + M13)
    rs = np.array([rs, rs])

    xs = np.array([x3 + rs[0], x3 - rs[0]])

    solutions = []

    for i in range(xs.shape[0]):
        if math.isclose((rs[i] - s1*r1)**2 - (xs[i] - x1)**2, 0, abs_tol=tolerance):
            solutions.append(Disk(center=(y1, xs[i]), radius=rs))
        elif (rs[i] - s1*r1)**2 - (xs[i] - x1)**2 < 0:
            continue
        else:
            ys1 = mp.sqrt((rs[i] - s1*r1)**2 - (xs[i] - x1)**2) + y1
            ys2 = -mp.sqrt((rs[i] - s1*r1)**2 - (xs[i] - x1)**2) + y1
            solutions.append(Disk(center=(xs[i], ys1), radius=rs[i]))
            solutions.append(Disk(center=(xs[i], ys2), radius=rs[i]))

    return list(filter(lambda x: x.radius > 0, solutions))

    # ys1 = np.array([y1, y1])
    # ys2 = np.array([y1, y1])
    # temp_array = (rs - s1*r1)**2 - (xs - x1)**2

    # if not math.isclose(temp_array[0], 0, abs_tol=tolerance):
    #     ys1[0] = np.sqrt((rs[0] - s1*r1)**2 - (xs[0] - x1)**2) + y1
    #     ys2[0] = -np.sqrt((rs[0] - s1*r1)**2 - (xs[0] - x1)**2) + y1

    # if not math.isclose(temp_array[1], 0, abs_tol=tolerance):
    #     ys1[1] = np.sqrt((rs[1] - s1*r1)**2 - (xs[1] - x1)**2) + y1
    #     ys2[1] = -np.sqrt((rs[1] - s1*r1)**2 - (xs[1] - x1)**2) + y1

    # rs = np.array([rs[0], rs[0], rs[0], rs[0]])
    # xs = np.array([xs[0], xs[1], xs[0], xs[1]])
    # ys = np.array([ys1[0], ys1[1], ys2[0], ys2[1]])

    # # res = [Disk(center=(xs[i], ys[i]), radius=abs(rs[i])) for i in range(4)]

    # return list(filter(lambda x: x.radius > 0, [Disk(center=(xs[i], ys[i]), radius=rs[i]) for i in range(4)]))


def circleTangentToCLL(c1: Disk, L1: Edge, L2: Edge, s1: int):
    # Unpack
    (x3, y3), (x4, y4) = (L1.start, L1.end)
    (x5, y5), (x6, y6) = (L2.start, L2.end)

    def is_parallel(L1, L2):  # Check if slopes are equal
        (x3, y3), (x4, y4) = (L1.start, L1.end)
        (x5, y5), (x6, y6) = (L2.start, L2.end)
        return math.isclose((x4-x3)*(y6-y5), (x6-x5)*(y4-y3), abs_tol=tolerance)

    # Line 1
    if is_parallel(L1, L2):
        A1 = y3 - y4
        B1 = x4 - x3
        C1 = (y4 - y3)*x3 - (x4 - x3)*y3
        C2 = (y6 - y5)*x5 - (x6 - x5)*y5
        A2 = y5 - y6
        B2 = x6 - x5

        if (x3 == x4):
            rs = abs(x3 - x5)/2
        else:
            # find projection, and angle to the projection, and rotate the other coordinates accordingly
            proj_onto_line = projection_scalar_for_line((0, 0), L1)
            xp = (1 - proj_onto_line) * x3 + proj_onto_line * x4
            yp = (1 - proj_onto_line) * y4 + proj_onto_line * y4
            if (xp == 0):
                if (yp > 0):
                    angle = mp.pi/2
                else:
                    angle = mp.pi*3/2
            else:
                angle = mp.atan(yp/xp)
            x3p = x3*mp.cos(angle) - y3*mp.sin(angle)
            x5p = x5*mp.cos(angle) - y5*mp.sin(angle)

            rs = abs(x3p - x5p)/2

        return circleTangentToCcenterL(c1, L1, L2, s1, rs)
    else:
        A1 = y3 - y4
        B1 = x4 - x3
        C1 = (y4 - y3)*x3 - (x4 - x3)*y3
        K1 = mp.sqrt((y4-y3)**2 + (x4-x3)**2)

        A2 = y5 - y6
        B2 = x6 - x5
        C2 = (y6 - y5)*x5 - (x6 - x5)*y5
        K2 = mp.sqrt((y6-y5)**2 + (x6-x5)**2)

    # L1 positive, L2 positive
    M11 = A1
    M12 = B1
    M13 = -K1
    b1 = -C1

    M21 = A2
    M22 = B2
    M23 = -K2
    b2 = -C2

    res1 = circleTangentToCLL_helper(
        c1, s1, M11, M12, M13, b1, M21, M22, M23, b2)
    res = list(filter(lambda tup: (A1*tup.center[0] + B1*tup.center[1] + C1)*tup.radius > 0 and (
        A2*tup.center[0] + B2*tup.center[1] + C2)*tup.radius > 0, res1))

    # L1 positive, L2 negative
    M11 = A1
    M12 = B1
    M13 = -K1
    b1 = -C1

    M21 = -A2
    M22 = -B2
    M23 = -K2
    b2 = C2

    res2 = circleTangentToCLL_helper(
        c1, s1, M11, M12, M13, b1, M21, M22, M23, b2)
    res += list(filter(lambda tup: (A1*tup.center[0] + B1*tup.center[1] + C1)*tup.radius >
                0 and (-A2*tup.center[0] - B2*tup.center[1] - C2)*tup.radius > 0, res2))

    # L1 negative, L2 positive
    M11 = -A1
    M12 = -B1
    M13 = -K1
    b1 = C1

    M21 = A2
    M22 = B2
    M23 = -K2
    b2 = -C2

    res3 = circleTangentToCLL_helper(
        c1, s1, M11, M12, M13, b1, M21, M22, M23, b2)
    res += list(filter(lambda tup: (-A1*tup.center[0] - B1*tup.center[1] - C1)*tup.radius > 0 and (
        A2*tup.center[0] + B2*tup.center[1] + C2)*tup.radius > 0, res3))

    # L1 negative, L2 negative
    M11 = -A1
    M12 = -B1
    M13 = -K1
    b1 = C1

    M21 = -A2
    M22 = -B2
    M23 = -K2
    b2 = C2

    res4 = circleTangentToCLL_helper(
        c1, s1, M11, M12, M13, b1, M21, M22, M23, b2)
    res += list(filter(lambda tup: (-A1*tup.center[0] - B1*tup.center[1] - C1)*tup.radius >
                0 and (-A2*tup.center[0] - B2*tup.center[1] - C2)*tup.radius > 0, res4))

    return res


def circleTangentToCLL_helper(c1, s1, M11, M12, M13, b1, M21, M22, M23, b2):
    # Unpack
    (x1, y1), r1 = (c1.center, c1.radius)

    M13p = (M12*M23-M13*M22)/(M22*M11-M21*M12)
    Mb1p = - (M12*b2-M22*b1)/(M22*M11-M21*M12)
    M23p = - (M23*M11-M21*M13)/(M22*M11-M21*M12)
    Mb2p = (b2*M11-M21*b1)/(M22*M11-M21*M12)

    # Second degree equation
    a = (M13p**2 + M23p**2 - 1)
    b = 2*Mb1p*M13p + 2*Mb2p*M23p - 2*M13p*x1 - 2*M23p*y1 + 2*s1*r1
    c = Mb1p**2 + Mb2p**2 + x1**2 + y1**2 - 2*Mb1p*x1 - 2*Mb2p*y1 - s1**2*r1**2

    D = b**2 - 4*a*c
    if math.isclose(D, 0, abs_tol=tolerance):
        D = 0
    if (D < 0):
        return []
    elif (a == 0):
        rs = np.array([-c / b])
        length = 1
    else:
        rs = np.array([(-b + mp.sqrt(D))/(2*a), (-b - mp.sqrt(D))/(2*a)])
        length = rs.size

    xs = Mb1p + M13p*rs  # should be array
    ys = Mb2p + M23p*rs

    # res = [Disk(center=(xs[i], ys[i]), radius=rs[i]) for i in range(length)]

    return list(filter(lambda x: x.radius > 0, [Disk(center=(xs[i], ys[i]), radius=rs[i]) for i in range(length)]))


def circleTangentToCcenterL(c1, L1, L2, s1, rs):
    # Unpack
    (x1, y1), r1 = (c1.center, c1.radius)
    (x3, y3), (x4, y4) = (L1.start, L1.end)
    (x5, y5), (x6, y6) = (L2.start, L2.end)

    new_line = Edge(start=((x3+x5)/2, (y3+y5)/2), end=((x4+x6)/2, (y4+y6)/2))
    (x7, y7), (x8, y8) = (new_line.start, new_line.end)

    A1 = y7 - y8
    B1 = x8 - x7
    C1 = (y8 - y7)*x7 - (x8 - x7)*y7

    if math.isclose(A1, 0, abs_tol=tolerance) and math.isclose(B1, 0, abs_tol=tolerance):
        return []  # This could be two vertical lines, but without being ''close'' to each other

    if math.isclose(A1, 0, abs_tol=tolerance):
        ys = -C1/B1
        a = 1
        b = -2*x1
        c = ys**2 + x1**2 + y1**2 - 2*ys*y1 - (rs**2 + r1**2 - 2*s1*r1*rs)

        D = b**2 - 4*a*c
        if math.isclose(D, 0, abs_tol=tolerance):
            D = 0
        if (D < 0):
            return []
        elif (a == 0):
            xs = np.array([-c / b])
            length = 1
        else:
            xs = np.array([(-b + mp.sqrt(D))/(2*a),
                          (-b - mp.sqrt(D))/(2*a)])
            length = xs.size
        ys = np.array([ys, ys])
    elif math.isclose(B1, 0, abs_tol=tolerance):
        xs = -C1/A1
        a = 1
        b = -2*x1
        c = xs**2 + x1**2 + y1**2 - 2*xs*x1 - (rs**2 + r1**2 - 2*s1*r1*rs)

        D = b**2 - 4*a*c
        if math.isclose(D, 0, abs_tol=tolerance):
            D = 0
        if (D < 0):
            return []
        elif (a == 0):
            ys = np.array([-c / b])
            length = 1
        else:
            ys = np.array([(-b + mp.sqrt(D))/(2*a),
                          (-b - mp.sqrt(D))/(2*a)])
            length = ys.size
        xs = np.array([xs, xs])
    else:
        # second degree equation
        a = (1 + (A1/B1)**2)
        b = 2*((A1/B1)*(C1/B1) - x1 + (A1/B1)*y1)
        c = x1**2 + y1**2 + 2*y1*(C1/B1) + (C1/B1)**2 - \
            (rs**2 + r1**2 - 2*s1*r1*rs)

        D = b**2 - 4*a*c
        if math.isclose(D, 0, abs_tol=tolerance):
            D = 0
        if (D < 0):
            return []
        elif (a == 0):
            xs = np.array([-c / b])
            length = 1
        else:
            xs = np.array([(-b + mp.sqrt(D))/(2*a),
                          (-b - mp.sqrt(D))/(2*a)])
            length = xs.size

        ys = -(A1/B1)*xs - C1/B1

    rs = np.array([rs, rs])
    # res = [Disk(center=(xs[i], ys[i]), radius=rs[i]) for i in range(length)]
    # res = list(filter(lambda tup: (A1*tup.center[0] + B1*tup.center[1] + C1)*tup.radius > 0 and (A2*tup.center[0] + B2*tup.center[1] + C2)*tup.radius > 0, res))

    return list(filter(lambda x: x.radius > 0, [Disk(center=(xs[i], ys[i]), radius=rs[i]) for i in range(length)]))


def circleTangentToLLL(L1: Edge, L2: Edge, L3: Edge):
    # Unpack
    (x1, y1), (x2, y2) = (L1.start, L1.end)
    (x3, y3), (x4, y4) = (L2.start, L2.end)
    (x5, y5), (x6, y6) = (L3.start, L3.end)

    if (math.isclose(x1, x3) and math.isclose(y1, y3) and math.isclose(x2, x4) and math.isclose(y2, y4)):
        return []
    if (math.isclose(x3, x5) and math.isclose(y3, y5) and math.isclose(x4, x6) and math.isclose(y4, y6)):
        return []
    if (math.isclose(x5, x1) and math.isclose(y5, y1) and math.isclose(x6, x2) and math.isclose(y6, y2)):
        return []

    def is_parallel(L1, L2):  # Check if slopes are equal
        return math.isclose((x4-x3)*(y6-y5), (x6-x5)*(y4-y3), abs_tol=tolerance)

    if (is_parallel(L1, L2) and is_parallel(L2, L3)):
        return []

    # Equation for line 1
    A1 = y1 - y2
    B1 = x2 - x1
    C1 = (y2 - y1)*x1 - (x2 - x1)*y1
    # Denominator
    K1 = mp.sqrt((y2-y1)**2 + (x2-x1)**2)
    # Equation for line 2
    A2 = y3 - y4
    B2 = x4 - x3
    C2 = (y4 - y3)*x3 - (x4 - x3)*y3
    # Denominator
    K2 = mp.sqrt((y4-y3)**2 + (x4-x3)**2)
    # Equation for line 13
    A3 = y5 - y6
    B3 = x6 - x5
    C3 = (y6 - y5)*x5 - (x6 - x5)*y5
    # Denominator
    K3 = mp.sqrt((y6-y5)**2 + (x6-x5)**2)

    # L1 pos, pos, pos, pos, neg, neg, neg, neg
    M11 = np.array([A1, A1, A1, A1, -A1, -A1, -A1, -A1])
    M12 = np.array([B1, B1, B1, B1, -B1, -B1, -B1, -B1])
    M13 = np.array([-K1, -K1, -K1, -K1, -K1, -K1, -K1, -K1])
    b1 = np.array([-C1, -C1, -C1, -C1, C1, C1, C1, C1])

    # L2 pos, pos, neg, neg, pos, pos, neg, neg
    M21 = np.array([A2, A2, -A2, -A2, A2, A2, -A2, -A2])
    M22 = np.array([B2, B2, -B2, -B2, B2, B2, -B2, -B2])
    M23 = np.array([-K2, -K2, -K2, -K2, -K2, -K2, -K2, -K2])
    b2 = np.array([-C2, -C2, C2, C2, -C2, -C2, C2, C2])

    # L3 pos, neg, pos, neg, pos, neg, pos, neg
    M31 = np.array([A3, -A3, A3, -A3, A3, -A3, A3, -A3])
    M32 = np.array([B3, -B3, B3, -B3, B3, -B3, B3, -B3])
    M33 = np.array([-K3, -K3, -K3, -K3, -K3, -K3, -K3, -K3])
    b3 = np.array([-C3, C3, -C3, C3, -C3, C3, -C3, C3])

    xs = (M12*M23*b3 - M12*M33*b2 - M13*M22*b3 + M13*M32*b2 + M22*M33*b1 - M23*M32*b1) / \
        (M11*M22*M33 - M11*M23*M32 - M12*M21*M33 +
         M12*M23*M31 + M13*M21*M32 - M13*M22*M31)
    ys = - (M11*M23*b3 - M11*M33*b2 - M13*M21*b3 + M13*M31*b2 + M21*M33*b1 - M23*M31*b1) / \
           (M11*M22*M33 - M11*M23*M32 - M12*M21*M33 +
            M12*M23*M31 + M13*M21*M32 - M13*M22*M31)
    rs = (M11*M22*b3 - M11*M32*b2 - M12*M21*b3 + M12*M31*b2 + M21*M32*b1 - M22*M31*b1) / \
        (M11*M22*M33 - M11*M23*M32 - M12*M21*M33 +
         M12*M23*M31 + M13*M21*M32 - M13*M22*M31)

    # res = [Disk(center=(xs[i], ys[i]), radius=abs(rs[i])) for i in range(8)]

    return list(filter(lambda x: x.radius > 0, [Disk(center=(xs[i], ys[i]), radius=rs[i]) for i in range(8)]))
