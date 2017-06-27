from typing import Tuple, List, Optional, Union, Callable
from copy import copy
import queue

from math import modulo_inverse, gcd, prime_test


class Polynomial(object):
    """
    A filter value on length of the polynomial,
    when linear-time optimizations kick in:
        power shifts for x^n polynomials
        poly-to-int optimizations during arithmetic ops (only for GF(2^m))
    """
    LINEAR_TIME_OPTIMIZATIONS_FROM_N = 64
    """
    A filter value on length of a polynomial
    when O(n*log(n)) should kick in:
        DFT multiplication
    """
    NLOGN_TIME_OPTIMIZATIONS_FROM_N = 128


    @property
    def ONE(self):
        return self.from_values([1])

    @property
    def ZERO(self):
        return self.from_values([0])

    @property
    def X(self):
        return self.from_values([0,1])

    def __init__(self, values: List[int], normalize: bool = False,
                 mod_or_provider: Union[int,Callable[[], int]] = None,
                 laziness: bool = False):
        self._values = list(values)
        # laziness - an operation (multiplication, division, addition by a constant), which
        # performs uniformly on all polynomial's coefficients, computes only in time, when the exact
        # value is necessary by the code.
        # this can slow down executions during intensive polynomial full-scans
        # but it can be mitigated by 'collapsing' lazy queue - compute all values and empty lazy queue
        # one advantage of this - it can use this additional knowledge on operation to
        # make optimizations on cross-polynomial operations (poly sum, poly mult, poly mod)
        # example:
        # P(x) = 4*3*(a0 + a1*x + a2*x**2 + ....) mod 6
        # then, in naive approach we will make O(n) operations
        # but using lazy queue we can make O(1) operations (3 exactly)
        # to determine, that the result will be 0
        self._laziness = laziness
        # stores operations to be performed in a format (value, opcode)
        # opcodes:
        # 0 - summation of all coefficients by a constant
        # 1 - multiplication of all coefficients by a constant (and modulo inverse)
        self._lazy_queue = queue.deque() if laziness else None
        # check if mod is provider in any way (exclusively)
        self._mod = mod_or_provider
        if normalize:
            self.normalize()

    @property
    def lazy(self):
        return self._laziness

    @lazy.setter
    def lazy(self, val: bool):
        if not val:
            # disable lazy queue - collapse all operations
            self.collapse_laziness(levels=-1)
            self._lazy_queue = None  # remove queue to free up space
        else:
            # enable lazy queue
            self._lazy_queue = queue.deque()

        self._laziness = val

    @property
    def lazed(self):
        return self._lazy_queue is not None and len(self._lazy_queue) > 0

    @property
    def mod(self):
        return self._mod if self._mod is int else self._mod()

    @property
    def deg(self):
        self.compact()
        return len(self)-1

    @property
    def is_x_pow(self):
        """
        :return: True if this polynomial has this form: A*x^n (or an A), A != 0
        """
        self.compact()
        if self[-1] == 0:
            return False
        # case for A
        if len(self) == 1:
            return True
        # if a non-zero coef a_i exists in
        # P(x) = Ax^n + a_n-1*x^n-1 + ... + a_1*x + a_0
        # then - it isn't x_pow
        for i in range(len(self)-1):
            if self[i] != 0:
                return False
        return True

    @classmethod
    def mod_check(cls, mod: int, greater: Optional[int] = None,
                  coprime: Optional[int] = None, prime: bool = False,
                  raise_error=True):
        if greater is not None and mod <= greater:
            if raise_error: raise ValueError(f"Modulo {mod} should be greater than {greater}.")
            return False
        if prime is not None and not prime_test(mod):
            if raise_error: raise ValueError(f"Modulo {mod} is not prime.")
            return False
        if coprime is not None and gcd(coprime, mod) !=1:
            if raise_error: raise ValueError(f"Modulo {mod} and {coprime} are not co-prime")
            return False

    def from_values(self, values: List[int]) -> 'Polynomial':
        return Polynomial(values, laziness=self.lazy,
                          mod_or_provider=self._mod,
                          normalize=False)

    def normalize(self) -> int:
        """
        Performs two operations: compacts polynomial and normalizes it
        (makes, so the leading coefficient is 1)
        :return: Leading coefficient value
        """
        self.compact()
        ret = self._values[-1]
        if ret and self != self.ZERO:
            self.mod_check(self.mod, coprime=ret)
            inv = modulo_inverse(self._values[-1], self.mod)
            self.mult_coefs(inv)
        return ret

    def mult_coefs(self, coef: int):
        coef %= self.mod
        if not self.lazy:
            self._compute_mult(coef)
        else:
            self._lazy_queue.append((coef, 1))

    def sum_coefs(self, coef: int):
        coef %= self.mod
        if not self.lazy:
            self._compute_sum(coef)
        else:
            self._lazy_queue.append((coef, 0))

    def compute(self, x: int):
        y = 0
        m = self.mod
        for i in range(len(self._values)):
            y = (y + (pow(x, i, m) * self[i] % m)) % m
        return y

    def compact(self):
        while self._values[-1] == 0:
            self._values.pop()

    def collapse_laziness(self, levels=-1):
        if levels < 0:
            levels = len(self._lazy_queue)
        levels = min(levels, len(self._lazy_queue))
        if levels == 0:
            return
        # if number of operations exceeds square of length - make
        # dumb full cycle
        if levels > len(self)**2:
            while len(self._lazy_queue) != 0:
                val, op = self._lazy_queue.popleft()
                if op == 0:
                    self._compute_sum(val)
                if op == 1:
                    self._compute_mult(val)

            # exit
            return
        # optimized cycle
        compacted, compacted_op = self._lazy_queue.popleft()
        i = 1
        m = self.mod
        act = lambda v, op: self._compute_mult(v) if op == 1 else self._compute_sum(v)

        while len(self._lazy_queue) and i < levels:
            ncomp, ncomp_op = self._lazy_queue.popleft()
            if ncomp_op == compacted_op:
                if compacted_op == 1:
                    compacted = (compacted * ncomp) % m
                else:
                    compacted = (compacted + ncomp) % m
            else:
                act(compacted, compacted_op)
        act(compacted, compacted_op)

    def gcd(self, other: 'Polynomial'):

        if self == self.ONE or other == self.ONE:
            return self.ONE
        if self == other:
            return copy(other)

        if self < other:
            return self.gcd(other % self)

    def __len__(self):
        return len(self._values)

    def __getitem__(self, i: int):
        assert i is int
        if self.lazy and self.lazed:
            return self.__lazy_queue_intercept(i)
        else:
            return self._values[i]

    def __setitem__(self, i: int, value: int):
        assert i is int
        assert value is int
        self._values[i] = value

    def __copy__(self):
        cp = Polynomial(self._values, self.mod, self._laziness)
        cp._lazy_queue = copy(self._lazy_queue)
        return cp

    def __sub__(self, other: 'Polynomial'):
        assert isinstance(other, Polynomial)

    def __neg__(self):
        return self.mult_coefs(-1)

    def __add__(self, other: Union['Polynomial', int]):
        assert isinstance(other, Polynomial) or isinstance(other, int)
        if isinstance(other, int):
            cp = copy(self)
            return cp._compute_sum(other)
        elif isinstance(other, Polynomial):
            mn = min(self, other)
            cp = copy(max(self, other))  # type: 'Polynomial'
            if self.mod != other.mod:
                raise ValueError("Mods are not equal")
            else:
                mod = min(self.mod, other.mod)
            limit = len(mn)
            if len(mn) == len(cp):
                while (cp[-1] + mn[-1]) % mod == 0:
                    cp._values.pop()
            else:
                for i in range(limit):
                    cp[i] = (cp[i] + mn[i]) % mod
            return cp
        else:
            raise TypeError(f"Type of other {type(other)} is not 'int' nor 'Polynomial'")

    def __mul__(self, other: Union['Polynomial', int]):
        assert isinstance(other, Polynomial) or isinstance(other, int)
        if isinstance(other, int):
            cp = copy(self)
            cp.mult_coefs(other)
            return cp
        elif isinstance(other, Polynomial):
            # naive
            m = self.mod
            values = list(0 for _ in range(len(self)+len(other)-1))
            for i in range(len(values)):
                val = 0
                for j in range(i+1):
                    val = (val + self[j] + other[i-j-1]) % m
                values[i] = val
            return Polynomial(values, False, self.mod, self.lazy)

    def __eq__(self, other: 'Polynomial'):
        if self.mod != other.mod:
            raise ValueError(f"Mods are not equal (self={self.mod}, other={other.mod})")
        self.compact()
        other.compact()

        if len(self) != len(other):
            return False
        m = self.mod
        for i in range(len(self)):
            if (self[i] - other[i]) % m != 0:
                return False
        return True


    def __pow__(self, power: int):
        # naive with binary
        if power == 0:
            return self.ONE

        val = copy(self)
        if power == 1:
            return val
        res = val
        power >>= 1
        # double and multiply
        while power != 0:
            if power & 1:
                res *= val
            power >>= 1
            val *= val
        return res


    def __mod__(self, other: 'Polynomial'):
        # naive table method
        _, mod = self._compute_poly_divide(other)
        return mod



    def __floordiv__(self, other: 'Polynomial'):
        raise NotImplementedError

    def _compute_mult(self, coef: int):
        # pop all the greater coefs
        done = False
        while not done:
            c = (self._values[-1] * coef) % self.mod
            if c == 0:
                self._values.pop()
            else:
                done = True
                self._values[-1] = c
        # compute all that remains
        for i in range(1, len(self._values)):
            self._values[-1-i] = (self._values[-1-i]*coef) % self.mod

    def _compute_sum(self, coef: int):
        if coef == 0:
            return
        for i in range(len(self._values)):
            self._values[i] += coef

    def _compute_poly_divide(self, other: 'Polynomial') -> Tuple['Polynomial', 'Polynomial']:
        self.poly_sanity(other)
        if self.is_zero:
            return self.ZERO, self.ZERO
        elif other.is_zero:
            raise ZeroDivisionError("Cannot divide by zero")
        elif other.is_one:
            return copy(self), self.ZERO
        elif other.is_coef:
            cp = copy(self)
            cp._compute_mult(modulo_inverse(other[0], self.mod))
            return cp, self.ZERO
        elif other.is_x_pow:
            power = other.deg
            a = other[-1]
            div = self.from_values(self._values[power:])
            div._compute_mult(modulo_inverse(a, self.mod))
            mod = self.from_values(self._values[:power])
            return div, mod
        else:
            cp = copy(self)
            div_vals = []
            m = self.mod
            while cp.deg > other.deg:
                shift = cp.deg - other.deg
                coef = cp[-1]*modulo_inverse(other[-1], m)
                deleter = other.shift()*coef
                



    def __lazy_queue_intercept(self, i: int):
        raise NotImplementedError("Skip lists are not implemented yet, so no"
                                  "targeted collapsing can be made")

