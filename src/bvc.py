import numpy as np
import geo_helper

class BVC:
	# static variables
	NUM_VER_PRE_ALLOCATE = 10

	def __init__(self, own_pos, world_corners = None):
		'''
		world_corners must be arranged in order: 2 x n (n is the # of corners)
		own_pos is used to determine the interior of the world
		'''

		m, self.n_world_corners = world_corners.shape
		self.world_corners = world_corners
		self.world_edges = np.zeros([self.n_world_corners,3]) # every row is a line equation

		# find the boundaries of the world
		n = self.n_world_corners
		for i in range(n-1):
			self.world_edges[i,:] = geo_helper.get_line(self.world_corners[:,i], self.world_corners[:,i+1], own_pos)
		self.world_edges[n-1:,:] = geo_helper.get_line(self.world_corners[:,n-1], self.world_corners[:,0], own_pos)

		self.ver = np.empty([2, self.NUM_VER_PRE_ALLOCATE]) # pre-allocate 20 vertices for the efficiency of pending
		self.lines = self.world_edges # BVC boundaries are initialized to be the world edges


	def update_bvc(self, own_pos, nbr_pos):
		# initial BVC cell
		self.ver[:, 0:self.n_world_corners] = self.world_corners
		self.lines = self.world_edges
		
		m, n_nbr = nbr_pos.shape
		dist_sq = [] # distance square
		for i in range(n_nbr):
			temp = own_pos - nbr_pos[:,i]
			dist_sq.append(np.linalg.norm(temp))

		








mypos = np.array([[1.0], 
				  [1.0]])
neighbors = np.array([[4, 3],
	                  [3, 4]])

corners = np.array([[0.0, 0.0, 5.0, 5.0], 
					[0.0, 5.0, 5.0, 0.0]])

bvc = BVC(mypos, corners)
bvc.update_bvc(mypos, neighbors)