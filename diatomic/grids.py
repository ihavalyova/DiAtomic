import numpy as np
from utils import C_bohr

__all__ = ['Grid']


class Grid:

    def __init__(self, npoints, rgrid, solver='sinc', alpha=0.0, rbar=0.0):

        self.ngrid = npoints
        self.rmin = rgrid[0] / C_bohr
        self.rmax = rgrid[1] / C_bohr
        rbar = rbar / C_bohr
        self.solver = solver.lower()

        self.Gy = np.ones(self.ngrid)
        self.Fy = np.zeros(self.ngrid)
        self.rgrid = self.generate_uniform_grid()

        if alpha > 0.0:
            # mapping is allowed with sinc method only
            self.solver = 'sinc'

            self.rmin = self.get_grid_bounding_values(self.rmin, rbar, alpha)
            self.rmax = self.get_grid_bounding_values(self.rmax, rbar, alpha)

            self.rgrid, ygrid = self.generate_nonuniform_grid(alpha, rbar)

            gy_power1 = np.power(1.0+ygrid, (1.0/alpha)-1.0)
            gy_power2 = np.power(1.0-ygrid, (1.0/alpha)+1.0)
            self.Gy = (2.0*rbar/alpha) * gy_power1 / gy_power2

            fy_power = (np.power((1.0 - np.power(ygrid, 2)), 2))
            self.Fy = (1.0 - (1.0/(alpha**2))) / fy_power

    def get_grid_points(self):

        return self.rgrid * C_bohr

    def get_grid_bounding_values(self, rlimit, rbar, alpha):

        return ((rlimit/rbar)**alpha - 1.0) / ((rlimit/rbar)**alpha + 1.0)

    def generate_uniform_grid(self):

        # FGH Fourier grid
        if self.solver == 'fourier':
            return np.linspace(self.rmin, self.rmax, num=self.ngrid,
                               endpoint=False)

        # FGH Sinc grid and FD5 grid
        return np.linspace(self.rmin, self.rmax, num=self.ngrid, endpoint=True)

    def generate_nonuniform_grid(self, alpha, rbar):

        ystep = (self.rmax - self.rmin) / (self.ngrid - 1)  # / ngrid - 1 ??
        # ygrid = np.ogrid[self.rmin+ystep:self.rmax+ystep:ystep]
        # ygrid = np.ogrid[self.rmin:self.rmax:ystep]
        # ygrid = np.linspace(self.rmin, self.rmax, num=self.ngrid)
        # ygrid = np.arange(self.rmin, self.rmax, step=ystep)
        # ygrid = np.linspace(
        # self.rmin, self.rmax, num=self.ngrid, endpoint=True
        # )

        ygrid = np.empty(self.ngrid)

        for j in range(1, self.ngrid+1):
            ygrid[j-1] = self.rmin + ystep*(j-1.0)

        Ry = rbar * np.power((1.0+ygrid) / (1.0-ygrid), 1.0/alpha)

        print(ygrid)
        print(len(ygrid))

        return Ry, ygrid
