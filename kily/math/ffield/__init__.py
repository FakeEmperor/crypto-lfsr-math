from typing import List, Union, Tuple, Callable
import math

from kily.math.poly import Polynomial
from kily.math import prime_test, modulo_inverse



class GFElement:

    def __init__(self, field: 'GaloisField', poly: 'Polynomial'):
        self._field = field
        self._poly = poly

    @property
    def poly(self) -> Polynomial:
        return self._poly

    @property
    def field(self):
        return self._field

    @classmethod
    def from_int(cls, field: 'GaloisField', value: int, strict: bool = False) -> 'GFElement':
        if field.mod <= value and strict:
            raise ValueError(f'[strict mode] Value {value} is greater or equal than mod of the field: {field}')
        return GFElement(field, Polynomial.from_binary(value, field.mod_provider))


class GaloisField:
    def __init__(self, mod: int, primary_element: GFElement,  ):
        self._mod = mod
        self._primary = primary_element

    @property
    def mod(self):
        return self._mod

    @property
    def mod_provider(self) -> Callable[[], int]:
        return lambda : self.mod

    @classmethod
    def from_size(cls, field_size: Union[int, Tuple[int, int]]):
        if field_size is int:
            if not prime_test(field_size):
                raise ValueError(f"Field Size '{field_size}' is not prime")
            num = field_size
        else:
            if

        




