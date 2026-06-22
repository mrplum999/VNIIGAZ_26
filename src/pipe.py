import math
from src.fluid import Fluid
from src.state import NodeState


class Pipe:
    """
    Универсальный класс для расчёта гидравлики труб (НКТ и шлейф).
    """

    def __init__(self, L, D, roughness, fluid: Fluid,
                 vertical_depth: float = 0.0, name: str = "pipe"):
        self.L             = L
        self.D             = D
        self.roughness     = roughness
        self.fluid         = fluid
        self.vertical_depth = vertical_depth
        self.name          = name

    def _friction_factor(self, Re: float, eps_D: float) -> float:
        Re = abs(float(Re))
        if Re < 1e-10:
            return 0.02
        if Re < 2300:
            return 64.0 / Re

        lam = 0.02
        for _ in range(100):
            inner   = (eps_D / 3.7) + (2.51 / (Re * math.sqrt(lam)))
            lam_new = (-2.0 * math.log10(inner)) ** (-2.0)
            if abs(lam_new - lam) < 1e-6:
                return lam_new
            lam = lam_new
        return lam

    def dp(self, P_in: float, q_std: float) -> NodeState:
        if q_std <= 0:
            return NodeState(
                name=self.name, P_in=P_in, P_out=P_in, dP=0.0,
                q_std=0.0, q_res=None, v=None, rho=None
            )

        P_avg = max(P_in, 1.0)

        for _ in range(3):
            P_avg = max(P_avg, 1.0)        
            rho  = self.fluid.ro(P_avg)
            mu   = self.fluid.mu(P_avg)
            Bg   = self.fluid.bg(P_avg)

            q_res = q_std * Bg
            area  = math.pi * (self.D / 2) ** 2
            v     = (q_res / 86400.0) / area

            Re    = abs(float(rho * v * self.D / (mu / 1000.0)))  # ← abs и float
            eps_D = self.roughness / self.D
            lam   = self._friction_factor(Re, eps_D)

            friction    = lam * (self.L / self.D) * (rho * v**2 / 2)
            hydrostatic = rho * 9.81 * self.vertical_depth

            delta_P_atm = (friction + hydrostatic) / 101325.0
            P_avg       = max(P_in - delta_P_atm / 2, 1.0)

        P_out = max(P_in - delta_P_atm, 1.0)

        return NodeState(
            name=self.name,
            P_in=P_in,
            P_out=P_out,
            dP=delta_P_atm,
            q_std=q_std,
            q_res=q_res,
            v=v,
            rho=rho,
        )
