from collections import deque

class AdaptiveTimestep:
    def __init__(self, dt) -> None:
        self.basic_dt = dt
        self.dt_deque = deque()
        
    def divide_timestep(self, dt_cur):
        self.dt_deque.extend([dt_cur / 2, dt_cur / 2])
    
    def get_timestep(self):
        if self.dt_deque:
            return self.dt_deque.pop()
        else:
            return self.basic_dt