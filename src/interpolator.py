class LinearInterpolator:
    """
    Линейная интерполяция.

    """

    def __init__(self, xs: list, ys: list):

        if len(xs) != len(ys):
            raise ValueError("Длины xs и ys должны совпадать.")

        if len(xs) < 2:
            raise ValueError("Необходимо минимум две точки.")

        # Проверяем сортировку
        for i in range(len(xs) - 1):
            if xs[i] >= xs[i + 1]:
                raise ValueError("Список xs должен быть отсортирован по возрастанию.")

        self.xs = xs
        self.ys = ys

    def predict(self, xp: float) -> float:
        """
        Возвращает значение функции в точке xp
        методом линейной интерполяции.
        """

        if xp < self.xs[0] or xp > self.xs[-1]:
            raise ValueError(
                f"Значение xp={xp} выходит за диапазон "
                f"[{self.xs[0]}, {self.xs[-1]}]"
            )

        # Если попали точно в узел таблицы
        for i in range(len(self.xs)):
            if xp == self.xs[i]:
                return self.ys[i]

        # Ищем нужный интервал
        for i in range(len(self.xs) - 1):

            x1 = self.xs[i]
            x2 = self.xs[i + 1]

            if x1 <= xp <= x2:

                y1 = self.ys[i]
                y2 = self.ys[i + 1]

                return y1 + (y2 - y1) * (xp - x1) / (x2 - x1)


raise RuntimeError("Ошибка интерполяции.")