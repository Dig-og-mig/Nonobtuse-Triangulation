from pydantic import BaseModel, Field, field_validator
from typing import Optional
from disk import Disk
import mpmath as mp
from mpmath import mpf
from typing import Tuple


class Edge(BaseModel):

    start: Tuple[mpf, mpf] = Field(
        ..., description="The start point of the edge.", min_items=2, max_items=2)
    end: Tuple[mpf, mpf] = Field(
        ..., description="The start point of the edge.", min_items=2, max_items=2)
    next: 'Optional[Edge]' = Field(
        None, description="The next edge in the polygon.")
    # prev : 'Optional[Edge]' = Field(None, description="The previous edge in the polygon.")
    hole: bool = Field(
        False, description="Whether the edge is part of a hole in the polygon.")

    def set_next(self, next_edge):
        self.next = next_edge

    # Used for PSLG only
    visited: int = Field(
        0, description="Whether the edge has been visited in the initialization.")

    disk_start: Optional[Disk] = Field(
        None, description="The disk at the start of the edge.")
    disk_end: 'Optional[Disk]' = Field(
        None, description="The disk at the end of the edge.")

    disks_start_end: list = Field(
        [], description="List of disks going from the start to the end of the edge.")

    disks_end_start: list = Field(
        [], description="List of disks going from the end to the start of the edge.")

    class Config:
        arbitrary_types_allowed = True
        json_encoders = {
            mpf: lambda v: str(v),  # Serialize mpf as a string
        }

    def __eq__(self, other):
        if isinstance(other, Edge):
            return self.start == other.start and self.end == other.end
        return False

    @field_validator("start", mode="before")
    def start_check(cls, v):
        return (mp.mpf(v[0]), mp.mpf(v[1]))

    @field_validator("end", mode="before")
    def end_check(cls, v):
        return (mp.mpf(v[0]), mp.mpf(v[1]))
