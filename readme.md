# About

Python code for doing error propagation pretty quick.

**Features**
* Automatic error propagation
* Compute weighted averages

# Quick start

The following is an example which calculates the gravitational acceleration from pendulum period and length measurements.

$$
g = \frac{4\pi^2 l}{T^2}
$$

Error propagation is done automatically. 

```python
>>> from errorprop import *
>>> measurements = np.array([2.77,2.81,2.93,2.95,2.49,2.81,2.95,2.76])
>>> T = MeasuredQuantity(measurements, err_sys=0.3)
>>> print(T)
2.809 ± 0.3 ± 0.1
>>> l = Quantity(1.84, err_sys=0.7E-2, err_stat=0)
>>> print(l)
1.84 ± 0.007
>>> g = 4*np.pi**2*l/T**2
>>> print(g)
9.208 ± 2.002 ± 0.656
```

See tutorial.ipynb for a thorough walkthrough.

# Ho to use

1. Download the `errorprop.py` file and put it into the same directory as your jupyter notebook.
2. Add the following import to your notebook:
```python
from errorprop import *
```
3. Start doing error propagation... ;-)

