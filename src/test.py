import numpy as np

class BVC:
	def __init__(self, world_corners = None):
		'''
		world_corners must be arranged in order
		'''
		ver = np.empty([2,20]) # pre-allocate 20 vertices for the efficiency of pending
		m, n = world_corners.shape
		world_edges = np.zeros([n,3]) # every row is a line equation

		# find the boundaries of the world
		for i in range(n-1):
			world_edges[:,i] = get_line(world_corners[:,i], world_corners[:,i+1]) # TODO: direction?
		world_edges[:,n-1] = get_line(world_corners[:,n-1], world_corners[:,0])



	def get_bvc(self):
		ver = world_corner # initial BVC cell

	def get_line(self, x1, x2, inner_pt = None):
		'''
		From two points, return normal and constant
		n'[x ;y] + a = 0
		INPUT: x1, x2 points as array([x,y]), inner_pt make sure n'*inner_pt + a <= 0
		OUTPUT: line, array([[n , a]])/norm(n)
		'''
		n0 = x2 - x1
		n = np.array( [ -n0[1] , n0[0] ] )
		a = np.dot(-n, x1)
		line = np.zeros([1,3])
		line[0,0:2] = n;
		line[0,2]  = a;
		line = line / np.linalg.norm(n)

		# decide which half plane is n'[x ;y] + a <= 0
		if(inner_pt is not None):
			if(np.dot(n, inner_pt) > 0):
				line = -line
			if(np.dot(n, inner_pt) == 0):
				print("[Warn]: cannot decide half plane in get_line()")

		return line




'''
bvc = BVC()

x1 = np.array([1,0])
x2 = np.array([2,1])
inner_pt = np.array([-1,1])
line = bvc.get_line(x1, x2)
print line
'''

a = np.ones([2,4])
m, n = a.shape
print n