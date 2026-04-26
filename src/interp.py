import numpy as np

class LinearInterpolator: 

    def __init__(self, xs, ys):
        self.xs = np.asarray(xs, dtype=float)
        self.ys = np.asarray(ys, dtype=float)

        if len(self.xs) != len(self.ys):
            raise ValueError("Списки xs и ys должны быть одинаковой длины")
        
        if len(self.xs) < 2:
            raise ValueError("Минимальное количество точек должно быть равно 2 для интерполяции")
        
        if not np.all(np.diff(self.xs) > 0):
            raise ValueError("xs надо отсортировыввать по возрастанию")

    def predict(self, xp: float) -> float:

        """интерп. в точке xp."""
        
        if xp < self.xs[0] or xp > self.xs[-1]:
            raise ValueError(f"Значение xp={xp} выходит за границы интерполяции "
                             f"[{self.xs[0]:.1f}, {self.xs[-1]:.1f}] атм")
        
        return float(np.interp(xp, self.xs, self.ys))