from pydantic import BaseModel, Field, field_validator
from typing import Optional
from .disk import Disk
# from mpmath import mpf
from typing import Tuple
# import mpmath as mp


class Point(BaseModel):

    point: Tuple[float,float] = Field(
        ..., description="The start point of the edge.", min_items=2, max_items=2)
    hole: bool = Field(
        False, description="Whether the pooint is part of a hole in the polygon.")

    # Used for PSLG only
    visited: int = Field(
        0, description="Whether the edge has been visited in the initialization.")

    disk: Optional[Disk] = Field(
        None, description="The disk centered at the point.")

    # class Config:
    #     arbitrary_types_allowed = True
    #     json_encoders = {
    #         mpf: lambda v: str(v)
    #     }

    # @field_validator("point", mode="before")
    # def point(cls, v):
    #     return (mp.mpf(v[0]), mp.mpf(v[1]))
