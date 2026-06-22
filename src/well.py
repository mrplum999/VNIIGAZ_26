import math
from src.fluid import Fluid
from src.pipe import Pipe


class Well:
    """
    Модель газовой скважины.
    """

    def __init__(
        self,
        fluid: Fluid,
        k: float,
        h: float,
        re: float,
        rw: float,
        pipe: Pipe = None,
        name: str = "well",
    ):
        self.name = name
        self.fluid = fluid

        self.k = k
        self.h = h
        self.re = re
        self.rw = rw

        self.pipe = pipe

        # Коэффициент Дарси
        self._beta = 0.00852702
        self._ln = math.log(re / rw)

    def q(self, P_res: float, P_bhp: float) -> float:
        """
        Дебит скважины (ст.м³/сут).

        P_res — пластовое давление, атм.
        P_bhp — забойное давление, атм.
        """

        if P_bhp >= P_res or P_res <= 0:
            return 0.0

        mu = self.fluid.mu(P_res)
        Bg = self.fluid.bg(P_res)

        productivity = (
            self._beta * self.k * self.h
        ) / (mu * self._ln)

        q_std = productivity * (P_res - P_bhp) / Bg

        return q_std

    def ipr(self, P_res: float, n_points: int = 50):
        """
        Построение кривой IPR.
        Возвращает список (q, P_bhp).
        """

        points = []

        for i in range(n_points + 1):
            P_bhp = P_res * i / n_points
            q = self.q(P_res, P_bhp)
            points.append((q, P_bhp))

        return points

    def vlp(
        self,
        P_man: float,
        q_max: float = 3000.0,
        n_points: int = 50,
    ):
        """
        Построение кривой VLP.
        Возвращает список (q, P_bhp).
        """

        if self.pipe is None:
            return []

        points = []

        for i in range(n_points + 1):
            q = q_max * i / n_points

            node = self.pipe.dp(P_man, q)

            P_bhp = P_man + node.dP

            points.append((q, P_bhp))

        return points

    def __str__(self):
        return (
            f"Well(name={self.name}, "
            f"k={self.k}, "
            f"h={self.h})"
        )