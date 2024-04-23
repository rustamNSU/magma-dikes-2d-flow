import numpy as np
from collections import deque
import scipy.linalg
import scipy.sparse
import scipy.sparse.linalg


class Mesh:
    def __init__(self, xmin, xmax, ymin, ymax, nx, ny) -> None:
        self.nx = nx
        self.ny = ny
        self.xmin = xmin
        self.xmax = xmax
        self.ymin = ymin
        self.ymax = ymax
        self.xb = np.linspace(xmin, xmax, nx + 1)
        self.yb = np.linspace(ymin, ymax, ny + 1)
        self.yc = 0.5 * (self.yb[:-1] + self.yb[1:])
        self.xc = 0.5 * (self.xb[:-1] + self.xb[1:])
        self.xc2d, self.yc2d = np.meshgrid(self.xc, self.yc, indexing='ij')
        

class EasyNS:
    def __init__(self, xmin, xmax, ymin, ymax, nx, ny) -> None:
        self.mesh = Mesh(xmin, xmax, ymin, ymax, nx, ny)
        nx = self.mesh.nx
        ny = self.mesh.ny
        self.p = np.zeros((nx, ny))
        self.u = np.zeros((nx+1, ny))
        self.v = np.zeros((nx, ny+1))
        self.set_ids()
    
    
    def set_ids(self) -> None:
        nx = self.mesh.nx
        ny = self.mesh.ny
        self.N = nx*ny + (nx+1)*ny + nx*(ny+1)
        n = self.N
        start = 0
        stop = nx*ny
        self.pid = np.arange(start=0, stop=stop, step=1).reshape(nx, ny)
        
        start = stop
        stop += (nx+1)*ny
        self.uid = np.arange(start=start, stop=stop, step=1).reshape(nx+1, ny)
        
        start = stop
        stop += nx*(ny+1)
        self.vid = np.arange(start=start, stop=stop, step=1).reshape(nx, ny+1)
    
    
    def assemble_system(self):
        A = np.zeros(self.N, self.N)
        rhs = np.zeros(self.N)
    
    
    def assemble_mass_balance(self, A, rhs):
        nx = self.mesh.nx
        ny = self.mesh.ny
        for i in range(nx):
            for j in range(ny):
                
    # def set_initial_state(self, T0, Cp, rho, mu0, beta0, t0, Tr, alpha, k):
    #     self.w2d  = np.full(self.shape2d, self.yb[1:] - self.yb[:-1])
    #     self.c2d  = np.zeros(self.shape2d)
    #     self.T2d  = np.full(self.shape2d, T0)
    #     self.convectT2d  = np.full(self.shape2d, T0)
    #     self.mu2d = np.full(self.shape2d, mu0)
    #     self.Cp2d = np.full(self.shape2d, Cp)
    #     self.rho2d = np.full(self.shape2d, rho)
    #     self.q2d  = np.zeros(self.shape2d)
    #     self.beta2d = np.full(self.shape2d, beta0)
    #     self.time = t0
    #     self.Tr = Tr
    #     self.alpha = alpha
    #     self.k = k
        
    
    # def define_mobility(self, ix):
    #     mu = self.mu2d[ix, :]
    #     c = np.zeros(self.yc.shape)
    #     yb = self.yb
    #     c[-1] = -yb[-1]**2 / mu[-1] / 2
    #     for iy in range(2, self.ny + 1):
    #         c[-iy] = c[-iy+1] + yb[-iy]**2 / 2 * (1 / mu[-iy+1] - 1 / mu[-iy])
        
    #     qy = (yb[1:]**3 - yb[:-1]**3) / 6 / mu + (yb[1:] - yb[:-1]) * c
    #     self.q2d[ix, :] = qy
    #     self.c2d[ix, :] = c
    #     self.q[ix] = np.sum(qy)
        
        
    # def define_dp(self, ix, Q):
    #     self.dp[ix] = Q / self.q[ix]
    #     self.p[-1] = self.pmax - self.dp[-1] * (self.xmax - self.xc[-1])
    #     for ix in range(2, self.nx + 1):
    #         self.p[-ix] = self.p[-ix + 1] - self.dp[-ix] * (self.xc[-ix+1] - self.xc[-ix])
        
        
    # def solve_convect_heat_equation(self, fc_old):
    #     dt = self.time - fc_old.time
    #     dx = self.xb[1:] - self.xb[:-1]
    #     for ix in range(1, self.nx):
    #         Tl = self.T2d[ix-1, :]
    #         Vl = self.dp[ix-1] * self.q2d[ix-1, :] * dt
    #         ml = Vl * self.rho2d[ix-1, :]
            
    #         Tr = self.T2d[ix, :]
    #         Vr = self.dp[ix] * self.q2d[ix, :] * dt
    #         mr = Vr * self.rho2d[ix, :]
            
    #         Vold = fc_old.w2d[ix, :] * dx[ix]
    #         mold = fc_old.rho2d[ix, :] * Vold
    #         Told = fc_old.T2d[ix, :]
            
    #         Vnew = self.w2d[ix, :] * dx[ix]
    #         Tnew = self.T2d[ix, :]
            
    #         Vall = np.concatenate((Vold, Vl))
    #         Tall = np.concatenate((Told, Tl))
    #         indxs = np.argsort(Tall, kind='stable')    
            
    #         Vstack = deque(Vall[indxs][::-1])
    #         Tstack = deque(Tall[indxs][::-1])
            
    #         i = 1
    #         E = -Vr[-i] * Tr[-i]
    #         V = Vnew[-i] + Vr[-i]
    #         while Vstack:
    #             Vi = Vstack.pop()
    #             Ti = Tstack.pop()
    #             if V > Vi:
    #                 V = V - Vi
    #                 E = E + Ti * Vi
    #             elif i < self.ny:
    #                 E = E + Ti * V
    #                 Tnew[-i] = E / Vnew[-i]
    #                 Vstack.append(Vi - V)
    #                 Tstack.append(Ti)
    #                 i = i + 1
    #                 V = Vnew[-i] + Vr[-i]
    #                 E = -Vr[-i] * Tr[-i]
    #             else:   
    #                 V = V - Vi
    #                 E = E + Ti * Vi
                
    #         Tnew[-i] = E / Vnew[-i]
    #         self.convectT2d[ix, :] = Tnew
    
    
    # def solve_diffusion_heat_equation(self, fc_old):
    #     dt = self.time - fc_old.time
    #     for ix in range(0, self.nx):
    #         w = self.w2d[ix, :]
    #         yc = self.yc2d[ix, :]
    #         rho = self.rho2d[ix, :]
    #         cp = self.Cp2d[ix, :]
    #         k = self.k
    #         Told = self.convectT2d[ix, :]
    #         rhs = Told * rho * cp * w / dt
    #         A = np.diag(rho * cp * w / dt)
            
    #         # Add up stream 
    #         for iy in range(self.ny - 1):
    #             dy = yc[iy+1] - yc[iy]
    #             A[iy, iy+1] += -k / dy
    #             A[iy, iy] += k / dy
            
    #         # Add down stream
    #         for iy in range(1, self.ny):
    #             dy = yc[iy] - yc[iy-1]
    #             A[iy, iy] += k / dy
    #             A[iy, iy-1] += - k / dy
            
    #         # Add boundary condition
    #         A[-1, -1] += self.alpha
    #         rhs[-1] += self.alpha * self.Tr[ix]
            
    #         Tnew = np.linalg.solve(A, rhs)
    #         self.T2d[ix, :] = Tnew
            
    
    # def update_viscosity(self, mu_T):
    #     self.mu2d = mu_T(self.T2d)
        
    
    # def update_beta(self, beta_T):
    #     self.beta2d = beta_T(self.T2d)
    
if __name__ == "__main__":
    model = EasyNS(0, 1, 0, 1, 5, 3)
    a = 1