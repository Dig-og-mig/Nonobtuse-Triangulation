from pydantic import BaseModel, Field, field_validator
from mpmath import mpf
from typing import Tuple
import mpmath as mp


class Disk(BaseModel):
    center: Tuple[mpf, mpf] = Field(
        ..., description="The center of the circle.", default_factory=tuple)
    # ge = 0, but can't work rn, check solver.py
    radius: mpf = Field(..., description="The radius of the circle.")
    tangents: list = Field(
        default_factory=list, description="A list of points that represent the tangents of the circle.")
    hole: bool = Field(
        False, description="Whether the edge is part of a hole in the polygon.")

    class Config:
        arbitrary_types_allowed = True

    @field_validator("radius", mode="before")
    def radius_check(cls, r):
        return mpf(r)

    @field_validator("center", mode="before")
    def center_check(cls, v):
        return (mp.mpf(v[0]), mp.mpf(v[1]))
