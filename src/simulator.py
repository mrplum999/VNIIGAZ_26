# simulator.py

from typing import Dict, List, Optional

import pandas as pd
from scipy.optimize import fsolve

from src.compressor import DCS
from src.pipe import Pipe
from src.reservoir import Reservoir
from src.state import NodeState
from src.well import Well


class FieldSimulator:
    """
    Симулятор куста скважин.

    Объединяет:
    - пласт (Reservoir)
    - скважины (Well)
    - шлейф (Pipe)
    - ДКС (DCS)
    """

    def __init__(
        self,
        reservoir: Reservoir,
        wells: List[Well],
        shlyf: Pipe,
        dcs: DCS
    ):
        self.reservoir = reservoir
        self.wells = wells
        self.shlyf = shlyf
        self.dcs = dcs

    def solve(
        self,
        P_res: float,
        x0: Optional[List[float]] = None
    ) -> Dict[str, NodeState]:
        """
        Поиск рабочей точки системы.

        Неизвестные:
        q1, q2, q3, P_man
        """

        def equations(x):
            q1, q2, q3, P_man = x

            q1 = max(q1, 0.0)
            q2 = max(q2, 0.0)
            q3 = max(q3, 0.0)

            P_man = max(P_man, self.dcs.P_in() + 0.3)

            qs = [q1, q2, q3]

            F = []

            for i, well in enumerate(self.wells):
                node = well.pipe.dp(P_man, qs[i])

                P_bhp = P_man + node.dP

                F.append(
                    well.q(P_res, P_bhp) - qs[i]
                )

            q_total = q1 + q2 + q3 + self.dcs.q_ext

            node_shlyf = self.shlyf.dp(
                self.dcs.P_in(),
                q_total
            )

            F.append(
                P_man - (
                    self.dcs.P_in() +
                    node_shlyf.dP
                )
            )

            return F

        if x0 is None:
            x0 = [
                500.0,
                500.0,
                500.0,
                self.dcs.P_in() + 0.5
            ]

        sol = fsolve(
            equations,
            x0,
            xtol=1e-10,
            maxfev=5000
        )

        q1, q2, q3, P_man = sol

        q1 = max(q1, 0.0)
        q2 = max(q2, 0.0)
        q3 = max(q3, 0.0)

        P_man = max(
            P_man,
            self.dcs.P_in()
        )

        result: Dict[str, NodeState] = {}

        for i, well in enumerate(self.wells):

            q = [q1, q2, q3][i]

            node = well.pipe.dp(P_man, q)

            result[f"well_{i+1}"] = NodeState(
                name=f"well_{i+1}",
                P_in=P_man + node.dP,
                P_out=P_man,
                dP=node.dP,
                q_std=q,
                q_res=node.q_res,
                v=node.v,
                rho=node.rho
            )

        q_total = q1 + q2 + q3 + self.dcs.q_ext

        result["shlyf"] = self.shlyf.dp(
            self.dcs.P_in(),
            q_total
        )

        result["dcs"] = NodeState(
            name="dcs",
            P_in=self.dcs.P_in(),
            P_out=self.dcs.P_line,
            dP=self.dcs.P_line - self.dcs.P_in(),
            q_std=q_total,
            q_res=None,
            v=None,
            rho=None
        )

        return result

    def run(
        self,
        N_days: int,
        dt: float = 1.0
    ) -> pd.DataFrame:
        """
        Расчёт разработки по времени.
        """

        data = []

        P_res = self.reservoir.resprops.P

        Gp = 0.0

        last_x = None

        print(
            f"Запуск симуляции на {N_days} суток..."
        )

        for day in range(N_days):

            if day % 30 == 0:
                print(
                    f"День {day:4d} | "
                    f"P_res = {P_res:.2f} атм"
                )

            states = self.solve(
                P_res,
                x0=last_x
            )

            q1 = states["well_1"].q_std
            q2 = states["well_2"].q_std
            q3 = states["well_3"].q_std

            q_total = q1 + q2 + q3

            P_man = states["well_1"].P_out

            last_x = [
                q1,
                q2,
                q3,
                P_man
            ]

            Gp += q_total * dt

            data.append({
                "t": day,
                "P_res": round(P_res, 4),
                "P_man": round(P_man, 4),
                "q1": round(q1, 2),
                "q2": round(q2, 2),
                "q3": round(q3, 2),
                "q_total": round(q_total, 2),
                "Gp": round(Gp / 1000.0, 3)
            })

            P_res = self.reservoir.p2(
                q_total,
                dt
            )

            self.reservoir.resprops.P = P_res

        print(
            f"День {N_days:4d} | "
            f"P_res = {P_res:.2f} атм"
        )

        print("Симуляция завершена.\n")

        return pd.DataFrame(data)