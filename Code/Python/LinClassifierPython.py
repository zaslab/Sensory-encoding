# -*- coding: utf-8 -*-
"""
Created on Sun Mar  5 14:03:59 2023

@author: ZasLab
"""


import itertools
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.io
import sys

# Define conditions

from sklearn.ensemble import RandomForestClassifier


classifiers = {
    'Random_Forest': RandomForestClassifier(max_depth=20, n_estimators=300, max_features='log2'),
}
