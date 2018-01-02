import matplotlib.pyplot as plt
import numpy as np


grid = np.array([[1,2,3,4],[2,3,4,5],[6,7,8,9],[10,11,12,13]])
x = [1,2,3,4]
y = [10,20,30,40]

plt.contour(x,y,grid)
plt.show()
