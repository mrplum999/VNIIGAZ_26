import math
import pandas as pd
from pathlib import Path
from interp import LinearInterpolator


class Fluid:

    R     = 8.314     # газовая постоянная, Дж/(моль·К)
    P_STD = 101325    # стандартное давление, Па
    T_STD = 293.15    # стандартная температура, К

    def __init__(self, M: float, rho_c: float, xa: float, xy: float, T: float):
        """
        M     — молярная масса газа, кг/моль
        rho_c — критическая плотность, кг/м³
        xa    — мольная доля N₂
        xy    — мольная доля CO₂
        T     — температура системы, К (const везде)
        """
        self.M     = M
        self.rho_c = rho_c
        self.xa    = xa
        self.xy    = xy
        self.T     = T

        csv_path = Path(__file__).parent.parent / "interp_data.csv"
        df = pd.read_csv(csv_path, sep=';')
        df.columns = df.columns.str.strip()
        self._mu_interp = LinearInterpolator(
            xs=df["pressure"].tolist(),
            ys=df["viscosity"].tolist()
        )

    def z(self, P_atm: float) -> float:
        """
        Z-фактор по GERG-91 mod. (ГОСТ 30319.2-96).
        P_atm — давление в атм.
        """
        T = self.T

        xe = 1.0 - self.xa - self.xy          # мольная доля углеводородов
        P_mpa = P_atm * 0.101325              # атм в МПа

        Zc = 1.0 - (0.074 * self.rho_c
                    - 0.006
                    - 0.063 * self.xa
                    - 0.0575 * self.xy) ** 2

        Me = (24.05525 * self.rho_c * Zc
              - 28.0135 * self.xa
              - 44.01  * self.xy) / xe

        H = 128.64 + 47.479 * Me

        B1 = (-0.425468 + 2.865e-3 * T - 4.62073e-6 * T**2
              + (8.77118e-4 - 5.56281e-6 * T + 8.8151e-9  * T**2) * H
              + (-8.24747e-7 + 4.31436e-9 * T - 6.08319e-12 * T**2) * H**2)

        B2   = -0.1446 + 7.4091e-4 * T - 9.1195e-7 * T**2
        B3   = -0.86834 + 4.0376e-3 * T - 5.1657e-6 * T**2
        B23  = -0.339693 + 1.61176e-3 * T - 2.04429e-6 * T**2
        Bst  = 0.72 + 1.875e-5 * (320 - T)**2        # B_star

        Bm = (xe**2    * B1
              + xe * self.xa * Bst * (B1 + B2)
              - 1.73 * xe * self.xy * (B1 * B3)**0.5
              + self.xa**2 * B2
              + 2 * self.xa * self.xy * B23
              + self.xy**2 * B3)

        C1 = (-0.302488 + 1.95861e-3 * T - 3.16302e-6 * T**2
              + (6.46422e-4 - 4.22876e-6 * T + 6.88157e-9  * T**2) * H
              + (-3.32805e-7 + 2.2316e-9  * T - 3.67713e-12 * T**2) * H**2)

        C2   = 7.8498e-3 - 3.9895e-5 * T + 6.1187e-8 * T**2
        C3   = 2.0513e-3 + 3.4888e-5 * T - 8.3703e-8 * T**2
        C223 = 5.52066e-3 - 1.68609e-5 * T + 1.57169e-8 * T**2
        C233 = 3.58783e-3 + 8.06674e-6 * T - 3.25798e-8 * T**2
        Cst  = 0.92 + 0.0013 * (T - 270)             

        Cm = (xe**3 * C1
              + 3 * xe**2 * self.xa * Cst * (C1**2 * C2)**(1/3)
              + 2.76 * xe**2 * self.xy * (C1**2 * C3)**(1/3)
              + 3 * xe * self.xa**2 * Cst * (C1 * C2**2)**(1/3)
              + 6.6 * xe * self.xa * self.xy * (C1 * C2 * C3)**(1/3)
              + 2.76 * xe * self.xy**2 * (C1 * C3**2)**(1/3)
              + self.xa**3 * C2
              + 3 * self.xa**2 * self.xy * C223
              + 3 * self.xa * self.xy**2 * C233
              + self.xy**3 * C3)

        b  = 1e3 * P_mpa / (2.7715 * T)    
        B0 = b * Bm
        C0 = b**2 * Cm

        A0 = 1 + 1.5 * (B0 + C0)
        A1 = 1 + B0
        A2 = (A0 - (A0**2 - A1**3)**0.5)**(1/3)
        Z = (1 + A2 + A1 / A2) / 3

        return Z/ Zc    

    def ro(self, P_atm: float) -> float:
        """Плотность газа [кг/м³] при давлении P_atm [атм]."""
        P_pa = P_atm * 101325
        Z    = self.z(P_atm)
        return P_pa * self.M / (Z * self.R * self.T)

    def ro_std(self) -> float:
        """Плотность при стандартных условиях [кг/м³]. Z_std ≈ 1."""
        return self.P_STD * self.M / (1.0 * self.R * self.T_STD)

    def bg(self, P_atm: float) -> float:
        """
        Фактор объёма газа [м³/м³].
        Bg = (P_std_atm * Z * T) / (P * T_std)
        """
        Z = self.z(P_atm)
        return (1.0 * Z * self.T) / (P_atm * self.T_STD)

    def mu(self, P_atm: float) -> float:
        """Динамическая вязкость [сП] при давлении P_atm [атм]."""
        return self._mu_interp.predict(P_atm)