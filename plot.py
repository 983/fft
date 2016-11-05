import numpy as np
import numpy.linalg
import scipy.signal
from matplotlib import pyplot as plt
import time

values = open("fft.txt", "rb").read().decode("ascii").strip().split("\n")
values = list(map(float, values))
plt.plot(values)
plt.show()
