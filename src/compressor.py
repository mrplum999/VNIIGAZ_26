from src.state import NodeState
class DCS:
    """
    Дожимная компрессорная станция (ДКС).

    Повышает давление газа перед подачей
    в магистральный газопровод.
    """

    def __init__(
        self,
        CR: float = 1.5,
        P_line: float = 5.0,
        q_ext: float = 0.0
    ):
        """
        Parameters
        ----------
        CR : float
            Степень сжатия (>= 1).

        P_line : float
            Давление в магистральном газопроводе, атм.

        q_ext : float
            Расход стороннего газа, ст.м³/сут.
        """

        self.CR = max(CR, 1.0)
        self.P_line = P_line
        self.q_ext = q_ext

    def P_in(self) -> float:
        """
        Минимальное давление на входе в ДКС,
        необходимое для обеспечения давления
        в магистрали.
        """

        return self.P_line / self.CR

    def set_compression_ratio(self, new_CR: float):
        """
        Изменение степени сжатия.
        Используется при анализе чувствительности.
        """

        self.CR = max(new_CR, 1.0)