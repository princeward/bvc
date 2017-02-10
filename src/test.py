import numpy as np
import geo_helper
import bvc
import matplotlib.pyplot as plt

'''
# test with single robot
mypos = np.array([[2.5], 
				  [2.5]])
neighbors = np.array([[4.2, 3, 1, 2.5, 4.5],
	                  [3, 4, 1.5, 0.5, 4.5]])

corners = np.array([[0.0, 0.0, 5.0, 5.0], 
					[0.0, 5.0, 5.0, 0.0]])

cell = bvc.BVC(mypos, 0.5, corners)
cell.update_bvc(mypos, neighbors)
cell.plot_bvc()
cell.plot_show()
'''


'''
# test with multiple robots
safe_rad = 0.2
robots = np.array([[4.2, 3, 1, 2.5, 4.5],
	                  [3, 4, 1.5, 0.5, 4.5]])

corners = np.array([[0.0, 0.0, 5.0, 5.0], 
					[0.0, 5.0, 5.0, 0.0]])

cell = bvc.BVC(robots[:, 0], safe_rad, corners)
cell.update_bvc(robots[:, 0], robots[:, [1,2,3,4]])
cell.plot_bvc()

cell = bvc.BVC(robots[:, 1], safe_rad, corners)
cell.update_bvc(robots[:, 1], robots[:, [0,2,3,4]])
cell.plot_bvc()

cell = bvc.BVC(robots[:, 2], safe_rad, corners)
cell.update_bvc(robots[:, 2], robots[:, [0,1,3,4]])
cell.plot_bvc()

cell = bvc.BVC(robots[:, 3], safe_rad, corners)
cell.update_bvc(robots[:, 3], robots[:, [0,1,2,4]])
cell.plot_bvc()

cell = bvc.BVC(robots[:, 4], safe_rad, corners)
cell.update_bvc(robots[:, 4], robots[:, [0,1,2,3]])
cell.plot_bvc()

cell.plot_show()
'''

a = np.empty((2,1))
a = np.hstack( (a, np.array( [[1],[2]])) )
print a.shape