import numpy as np

a = [[1,2,3,np.nan,5,6],
     [1,2,3,4,5,np.nan]]
nan_mask = np.isnan(a)
b = np.array(a)[~nan_mask]