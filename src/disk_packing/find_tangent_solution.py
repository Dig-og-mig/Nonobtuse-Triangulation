from .make_disk_ccc import disk_ccc
from .make_disk_ccl import disk_ccl
from .make_disk_clc import disk_clc
from .make_disk_cll import disk_cll
from .make_disk_lll import disk_lll
from .make_disk_lcc import disk_lcc
from .make_disk_lcl import disk_lcl
from .make_disk_llc import disk_llc

from models.disk import Disk
from models.edge import Edge


def find_tangent_solution(elm1_tup, elm2_tup, elm3_tup):
    e1 = elm1_tup[1]
    e2 = elm2_tup[1]
    e3 = elm3_tup[1]

    if isinstance(e1, Disk) and isinstance(e2, Disk) and isinstance(e3, Disk):
        solution = disk_ccc(elm1_tup, elm2_tup, elm3_tup)
    elif isinstance(e1, Disk) and isinstance(e2, Disk) and isinstance(e3, Edge):
        solution = disk_ccl(elm1_tup, elm2_tup, elm3_tup)
    elif isinstance(e1, Disk) and isinstance(e2, Edge) and isinstance(e3, Disk):
        solution = disk_clc(elm1_tup, elm2_tup, elm3_tup)
    elif isinstance(e1, Disk) and isinstance(e2, Edge) and isinstance(e3, Edge):
        solution = disk_cll(elm1_tup, elm2_tup, elm3_tup)
    elif isinstance(e1, Edge) and isinstance(e2, Edge) and isinstance(e3, Edge):
        solution = disk_lll(elm1_tup, elm2_tup, elm3_tup)
    elif isinstance(e1, Edge) and isinstance(e2, Disk) and isinstance(e3, Disk):
        solution = disk_lcc(elm1_tup, elm2_tup, elm3_tup)
    elif isinstance(e1, Edge) and isinstance(e2, Disk) and isinstance(e3, Edge):
        solution = disk_lcl(elm1_tup, elm2_tup, elm3_tup)
    elif isinstance(e1, Edge) and isinstance(e2, Edge) and isinstance(e3, Disk):
        solution = disk_llc(elm1_tup, elm2_tup, elm3_tup)
    else:
        raise ValueError(
            "You are trying to find a tangent disk for an invalid combination of elements.")

    return solution
