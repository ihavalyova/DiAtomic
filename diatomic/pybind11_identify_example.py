#!/usr/bin/env python

import identify
import numpy as np

if __name__ == "__main__":

    nch = 1
    identify_option = 0
    clevels = np.array([[1, 3400.0, 2, 0, 1, 1, 1, 1, 1, 1],
                        [2, 1210.0, 1, 0, 1, 1, 1, 0, 1, 1],
                        [3, 3150.0, 2, 0, 1, 1, 1, 1, 1, 1],
                        [4, 4000.0, 4, 0, 1, 1, 1, 3, 1, 1]])

    elevels = np.array([[1, 0, 1, 1000.0, 0, 0.1, 3, 1],
                        [2, 1, 2, 2000.0, 0, 0.1, 3, 1]])

    answer = identify.identify_levels(clevels, elevels, nch, identify_option)
    print(answer)
    print(answer.shape)
