from typing import List, Tuple, Union

from kily.math.ffield import GaloisField, GFElement

class LinearFeedbackShiftRegister:
    """
    Class emulates LFSR and provides useful methods to find various properties of a given LFSR

    Feedback function's coefficients are written in the following order:
    f(x) = c_L + c_(L-1)*x + c_(L-2)*x^2
    """

    def __init__(self, field: GaloisField, coefficients: Union[int, GFElement], initial_fill: Union[int, GFElement]= None):
        self._time = 0
        self._field = field
        if isinstance(coefficients, int):
            self._coefficients = GFElement.from_int(field, coefficients, strict=True)
        else:
            self._coefficients = coefficients
        if isinstance(initial_fill, int):
            self._state = GFElement.from_int(field, initial_fill, strict=True)
        else:
            self._state = initial_fill
        self.check()

    def check(self, raise_error: bool = True) -> bool:
            raise NotImplementedError

    @property
    def time(self) -> int:
        return self._time

    @property
    def coefficients(self) -> GFElement:
        return self._coefficients

    @property
    def state(self) -> GFElement:
        return self._state

    def step(self) -> int:
        new_last_state = 0
        for i in range(len(self.state.poly)):
            new_last_state = (self.coefficients.poly[-1-i] * self.state.poly[i])  % self._field.mod

        out = self.shift(1)
        self.state.poly[0] = new_last_state
        self._time += 1
        return out

    def shift(self, offset) -> int:
        if offset == 0:
            return
        self.state.poly.shift(1)
        out = self.state.poly.pop()
        return out

    def find_minimal(self):
        """
        Berlekamp-Massey algorithm
        """
        raise NotImplementedError

    @classmethod
    def parse_coefs(cls, field: GaloisField, coefficients: Union[int, GFElement], strict: bool = False) -> GFElement:
        return
