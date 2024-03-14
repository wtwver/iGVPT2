import numpy as np

def A(v, a1, a2, a3, a4):
  return a1 * np.exp(-a2 * (v - a3)**2) / (1 + a4 * (v - a3))