# -*- coding: utf-8 -*-
"""
Created on Tue Apr 16 17:00:20 2024

@author: khajmoh18975
"""

import pvlib
from pvlib import pvsystem
import os
import pandas as pd

index_A = pd.date_range('2024-01-01', periods=10, freq='D')
index_B = pd.date_range('2024-01-01', periods=10, freq='2D')  # Different frequency for demonstration
A = pd.Series([1, 2, None, 4, None, 6, 7, None, 9, None], index=index_A)
B = pd.Series([None, 2, 3, None, 5, None, 7, 8, None, 10], index=index_B)

# Interpolate missing values in each data series
A_interpolated = A.interpolate(method='time')
B_interpolated = B.interpolate(method='time')

# Create a common index by merging the indices of A and B
common_index = A_interpolated.index.union(B_interpolated.index)

# Reindex both interpolated series to the common index
A_interpolated = A_interpolated.reindex(common_index)
B_interpolated = B_interpolated.reindex(common_index)

# Merge the interpolated series together
merged = pd.concat([A_interpolated, B_interpolated], axis=1)
merged.columns = ['A', 'B']  # Rename columns if needed

print(merged)