import numpy as np
import os
import sys


def best_fit(x,y):

    n = len(x)

    xm = x.mean(0)
    ym = y.mean(0)

    xa = x - xm
    ya = y - ym

    if len(ya.shape) == 3:
        xa = np.expand_dims(np.expand_dims(xa, -1), -1)

    xss = (xa ** 2).sum(0) / (n - 1)    # variance of x
    xys = (xa * ya).sum(0) / (n - 1)    # covariance

    slope = xys / xss
    intercept = ym - (slope * xm)

    return slope, intercept

