from dataclasses import dataclass
from src.fluid import Fluid

@dataclass
class ResProps:
    """
    Контейнер для свойств пласта.
    """
    P: float      # текущее пластовое давление [атм]
    V: float      # объём пласта [м³]
    T: float      # температура пласта [К]

class Reservoir:
    """
    Модель пласта. Отвечает только за материальный баланс.

    """
    def __init__(self, resprops: ResProps, fluid: Fluid):
        self.resprops = resprops
        self.fluid    = fluid

    def p2(self, q_total: float, dt: float = 1.0) -> float:
        if q_total <= 0:
            return self.resprops.P
        P_atm = self.resprops.P
        P_pa = P_atm * 101325.0
        rho_std = self.fluid.ro_std()          # кг/м³ при ст.усл.
        Z = self.fluid.z(P_atm)
        R = self.fluid.R                       # Дж/(моль·К) (8.314)
        T = self.resprops.T                    # К
        M = self.fluid.M                       # кг/моль
        V = self.resprops.V                    # м³
        dP_pa = (rho_std * Z * R * T) / (M * V) * q_total * dt
        dP_atm = dP_pa / 101325.0
        new_P = max(P_atm - dP_atm, 1.0)
        return new_P            