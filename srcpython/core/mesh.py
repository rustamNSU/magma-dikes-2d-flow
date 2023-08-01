import numpy as np

class Mesh2d:
    def __init__(self, xmin, xmax, width, nx, ny) -> None:
        self.nx = nx
        self.ny = ny
        self.xmin = xmin
        self.xmax = xmax
        self.width = wi