from typing import List

import numpy as np


class Quantity:
    def __init__(self, value, err_sys=.0, err_stat=.0):
        """
        :param data: array-like. First element is the value, following elements are errors (e.g. stat, sys)
        """
        self.value = value
        self.err_sys = err_sys
        self.err_stat = err_stat

    @property
    def errors(self):
        all_errors = [self.err_sys, self.err_stat]
        return np.array([error for error in all_errors if error > 0])

    @property
    def errors_rel(self):
        return 100*self.errors/self.value

    @property
    def error_rel(self):
        return self.errors_rel[0]

    def str(self, decimals=3, rel_errors=False):
        append_percent = lambda x: str(x) + '%'
        errors = self.errors_rel if rel_errors else self.errors
        errors = errors.round(decimals)
        # remove non zero errors (after rounding)
        errors = [error for error in errors if error > 0]
        errors = np.array(list(map(append_percent, errors))) if rel_errors else errors
        # value = self.value.round(decimals)
        value = round(self.value, decimals)
        if len(errors) > 0:
            return str(value) + " ± " + " ± ".join(map(str, errors))
        else:
            return str(value)

    def str_rel(self, decimals=2):
        return self.str(decimals=decimals, rel_errors=True)

    def scale(self, factor: float):
        return Quantity(self.value*factor, self.err_sys*factor, self.err_stat*factor)

    def to_arr(self):
        return np.array([self.value, self.err_sys, self.err_stat])

    def __try_to_convert_to_quantity__(self, obj):
        if isinstance(obj, Quantity):
            return obj
        if is_numeric(obj):
            return Quantity(obj, 0, 0)
        if isinstance(obj, np.ndarray):
            return np.array([self.__try_to_convert_to_quantity__(x) for x in obj])
        raise Exception(f"Can't convert {obj} to quantity")

    def sqrt(self):
        return self**0.5

    def __mul__(self, other):
        other = self.__try_to_convert_to_quantity__(other)
        # this represents a*[b,c,d,...] which should result in [ab, ac, ad, ...]
        if isinstance(other, np.ndarray):
            return np.array([self*x for x in other])
        else:
            # applying error propagation yields dc/da = b and dc/db = a
            # hence c_sys = |b*a_sys| + |a*b_sys| (similar for stat error)
            (a, b) = (self, other)
            c = lambda a,b: a*b
            c_a = lambda a, b: b
            c_b = lambda a, b: a
            return DerivedQuantity(c, [c_a, c_b], [a, b])

    __rmul__ = __mul__

    def __pow__(self, other):
        other = self.__try_to_convert_to_quantity__(other)
        (a, b) = (self, other)
        c = lambda a,b: a**b
        c_a = lambda a, b: b*a**(b-1)
        c_b = lambda a, b: np.log(a)*a**b
        return DerivedQuantity(c, [c_a, c_b], [a, b])

    def __rpow__(self, other):
        other = self.__try_to_convert_to_quantity__(other)
        return other**self

    def __truediv__(self, other):
        other = self.__try_to_convert_to_quantity__(other)
        # applying error propagation yields dc/da = b and dc/db = a
        # hence c_sys = |b*a_sys| + |a*b_sys| (similar for stat error)
        (a, b) = (self, other)
        c = lambda a,b: a/b
        c_a = lambda a, b: 1/b
        c_b = lambda a, b: -a/(b**2)
        return DerivedQuantity(c, [c_a, c_b], [a, b])

    def __rtruediv__(self, other):
        other = self.__try_to_convert_to_quantity__(other)
        return other/self

    def __add__(self, other):
        other = self.__try_to_convert_to_quantity__(other)
        # applying error propagation yields dc/da = 1 and dc/db = 1
        # hence c_sys = |a_sys| + |b_sys| (similar for stat error)
        (a, b) = (self, other)
        c = lambda a,b: a+b
        c_a = lambda a, b: 1
        c_b = lambda a, b: 1
        return DerivedQuantity(c, [c_a, c_b], [a, b])

    __radd__ = __add__

    def __sub__(self, other):
        return self.__add__(other*-1)

    def __rsub__(self, other):
        other = self.__try_to_convert_to_quantity__(other)
        return other-self

    def __str__(self):
        return self.str()

    def __repr__(self):
        return self.str()


    # def _ipython_display_(self):
    #     return

    def exp(self):
        # this function is called when you pass q=Quantity() to np.exp(q)
        # simple because exponential has itself as its derivative
        e = lambda x: np.exp(x)
        return DerivedQuantity(e, [e], [self])


class MeasuredQuantity(Quantity):
    def __init__(self, measurements, err_sys=0, sigma_interval=2):
        """
        :param measurements:
        :param err_sys:
        :param sigma_interval: sigma interval to use for statistical error (1=68.6%, 2=95.4%, 3=99.7%).
                               See IF0 for reference.
        References
        ----------
        https://ap.physik.unibas.ch/PDF/Manuals/German/IF0.pdf
        """
        self.measurements = np.array(measurements)
        self.err_sys = err_sys
        self.err_stat = 0
        N = len(measurements)
        errors = []
        if (self.err_sys > 0):
            errors.append(self.err_sys)
        if N > 1:
            std = self.measurements.std()
            # 2 sigma confidence interval
            self.err_stat = sigma_interval*std/np.sqrt(N)
            errors.append(self.err_stat)
        super().__init__(measurements.mean(), self.err_sys, self.err_stat)


class DerivedQuantity(Quantity):
    def __init__(self, value_function, partial_derivatives, args):
        err_sys = 0
        err_stat = 0
        quantity_values = [q_value(arg) for arg in args]
        for i in range(len(args)):
            quantity = args[i]
            partial_derivative = partial_derivatives[i]
            err_sys += np.abs(partial_derivative(*quantity_values)*q_sys_error(quantity))
            err_stat += (partial_derivative(*quantity_values)*q_stat_error(quantity))**2
        err_stat = np.sqrt(err_stat)
        super().__init__(value_function(*quantity_values), err_sys=err_sys, err_stat=err_stat)


class WeightedAverage(Quantity):
    """
    Weighted average over n error prone measurements x_i.
    Mean:    x = Σ(w_i*x_i)/w_i
    Weights: w_i = 1/Δx_i^2
    Error:   Δx = √1/(Σw_i)
    """
    def __init__(self, quantities: List[Quantity]):
        x_i = q_value(quantities)
        delta_x_i = q_error(quantities)
        w_i = 1/delta_x_i**2
        x_mean = (x_i*w_i).sum()/w_i.sum()
        delta_x = np.sqrt(1/w_i.sum())
        super().__init__(x_mean, err_stat=delta_x)


def create_list_of_quantities(data, err_sys=0, err_stat=0):
    if is_numeric(err_sys):
        err_sys = [err_sys]*len(data)
    if is_numeric(err_stat):
        err_stat = [err_stat]*len(data)
    output = []
    for i in range(len(data)):
        output.append(Quantity(data[i], err_sys=err_sys[i], err_stat=err_stat[i]))
    return np.array(output)


def apply_function(obj_or_list, func):
    """
    Apply a function either on a single object or on a list of objects
    :return: single object or array of objects, applying func on the object(s)
    """
    if isinstance(obj_or_list, list) or isinstance(obj_or_list, np.ndarray):
        return np.array([func(obj) for obj in obj_or_list])
    else:
        return func(obj_or_list)


def q_value(quantity):
    """
    Get value of quantity
    :param quantity: Quantity or List[Quantity]
    :return: float or np.array[float]
    """
    return apply_function(quantity, lambda q: q if is_numeric(q) else q.value)


def q_sys_error(quantity):
    """
    Get sys error (i.e. quantity.err_sys)
    :param quantity: Quantity or List[Quantity]
    :return: float or np.array[float]
    """
    return apply_function(quantity, lambda q: q.err_sys)


def q_stat_error(quantity):
    """
    Get stat errors (i.e. quantity.err_stat)
    :param quantity: Quantity or List[Quantity]
    :return: float or np.array[float]
    """
    return apply_function(quantity, lambda q: q.err_stat)


def q_error(quantity):
    """
    Get total error (stat and sys)
    :param quantity: Quantity or List[Quantity]
    :return: float or np.array[float]
    """
    return q_stat_error(quantity) + q_sys_error(quantity)


def is_numeric(x):
    """
    Checks if argument is numeric (int, float or any kind or numpy number type)
    :return: bool
    """
    return isinstance(x, float) or isinstance(x, int) or isinstance(x, np.number)
