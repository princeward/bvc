import numpy as np
import pdb # for debug
import geo_helper
import matplotlib.pyplot as plt

class BVC:
	# static variables
	NUM_VER_PRE_ALLOCATE = 10

	def __init__(self, own_pos, safe_rad, world_corners):
		'''
		world_corners must be arranged in order: 2 x n (n is the # of corners)
		own_pos is used to determine the interior of the world
		'''

		m, self.n_world_corners = world_corners.shape
		self.world_corners = world_corners
		self.world_edges = np.zeros([self.n_world_corners,3]) # every row is a line equation

		self.safe_rad = safe_rad

		# find the boundaries of the world
		n = self.n_world_corners
		for i in range(n-1):
			self.world_edges[i,:] = geo_helper.get_line(self.world_corners[:,i], self.world_corners[:,i+1], own_pos)
		self.world_edges[n-1:,:] = geo_helper.get_line(self.world_corners[:,n-1], self.world_corners[:,0], own_pos)

		self.ver = np.empty([2, self.NUM_VER_PRE_ALLOCATE]) # pre-allocate 20 vertices for the efficiency of pending
		self.lines = np.empty([self.NUM_VER_PRE_ALLOCATE, 3]) # pre-allocate


	def update_bvc(self, own_pos, nbr_pos):
		# FACT: voronoi cell always has the same number of vertices and edges
		self.own_pos = own_pos.reshape([2,1])
		self.nbr_pos = nbr_pos
		# initial BVC cell
		self.ver = np.empty([2, self.NUM_VER_PRE_ALLOCATE]) # pre-allocate 20 vertices for the efficiency of pending
		self.lines = np.empty([self.NUM_VER_PRE_ALLOCATE, 3]) # pre-allocate
		self.ver[:, 0:self.n_world_corners] = self.world_corners
		self.lines[0:self.n_world_corners, :] = self.world_edges
		
		num_line = self.n_world_corners  # number of lines
		num_ver = num_line # number of vertices
		
		m, n_nbr = nbr_pos.shape
		dist_sq = np.zeros(n_nbr) # distance square
		for i in range(n_nbr):
			temp = own_pos - nbr_pos[:,i]
			dist_sq[i] = np.linalg.norm(temp)
		idx = np.argsort(dist_sq) # ascending order
		nbr_pos = nbr_pos[:,idx] # sort nbr_pos in ascending order

		for i in range(n_nbr):
			buf_edge = geo_helper.get_one_bvc_edge(own_pos, nbr_pos[:,i], self.safe_rad)
			is_valid_edge = False	
			# find new intersection points with existing lines
			for j in range(num_line):
				int_pt = geo_helper.get_line_intersect(buf_edge, self.lines[j,:])
				if int_pt is not None:
					if geo_helper.is_in_hull(int_pt, self.lines[0:num_line, :]): # valid intersection point, append it as new vertice
						###################
						self.ver[:,num_ver] = int_pt.reshape(2) # append the vertex
						num_ver += 1
						is_valid_edge = True

			# process if the new buf edge is indeed a BVC edge
			if is_valid_edge:
				# append the new edge
				self.lines[num_line, :] = buf_edge
				num_line += 1
				# remove redundant vertices, get the new set of vertices
				old_ver = self.ver
				old_num_ver = num_ver
				self.ver = np.empty([2, self.NUM_VER_PRE_ALLOCATE]) # reset the vertices
				num_ver = 0
				for vv in range(old_num_ver):
					if geo_helper.is_in_hull(old_ver[:,vv], self.lines[0:num_line, :]):
						self.ver[:, num_ver] = old_ver[:,vv]
						num_ver += 1
				# remove redandent edges
				# criteria: if none of the vertices is on this edge -> redandent!
				old_lines = self.lines
				old_num_line = num_line
				self.lines = np.empty([self.NUM_VER_PRE_ALLOCATE, 3])
				num_line = 0
				for ll in range(old_num_line):
					is_useful = False
					for vv in range(num_ver):
						if geo_helper.is_on_line(self.ver[:,vv], old_lines[ll,:]):
							is_useful = True
					if is_useful:
						self.lines[num_line, :] = old_lines[ll,:]
						num_line += 1

			if num_line != num_ver:
				print "[Err]: something wrong in update_bvc()."

		self.N = num_ver # number of edges and vertices

		# put the vertices in order
		angles = np.zeros(self.N)
		for i in range(self.N):
			angles[i] = np.arctan2(self.ver[1,i]-self.own_pos[1,0],  self.ver[0,i]-self.own_pos[0,0])
		angIdx = np.argsort(angles)
		self.ver = self.ver[:,angIdx]

	def find_closest_to_goal(self, goal):
		'''
		find the closest point on the cell to a goal point
		'''
		goal = np.reshape(goal,(2,1))
		# 1. check if the goal is already inside BVC
		if geo_helper.is_in_hull(goal, self.lines[0:self.N, :]): # self.lines has redandent pre-allocated elements, need to disregard
			return goal

		# 2. check all vertices
		dist_min = float('inf')
		pt_min = np.zeros((2,1))
		for i in range(self.N):
			dist = np.linalg.norm(self.ver[:,i] - goal)
			if dist < dist_min:
				dist_min = dist
				pt_min = self.ver[:,i]

		# 3. check projections to edges
		for i in range(self.N):
			proj_pt = geo_helper.project_pt_to_line(goal, self.lines[i,:])
			if proj_pt is None: # it's possible there is no intersection point
				continue
			if geo_helper.is_in_hull(proj_pt, self.lines[0:self.N, :]):
				dist = np.linalg.norm(proj_pt - goal)
				if dist < dist_min:
					dist_min = dist
					pt_min = proj_pt

		return pt_min


	def plot_bvc(self):
		# draw own position
		plt.plot(self.own_pos[0], self.own_pos[1],'b^')
		# draw neighbors
		plt.plot(self.nbr_pos[0,:], self.nbr_pos[1,:],'b^')
		# draw BVC
		plt.plot(self.ver[0, 0:self.N], self.ver[1, 0:self.N], 'k')
		plt.plot([self.ver[0,self.N-1], self.ver[0,0]], [self.ver[1,self.N-1], self.ver[1,0]], 'k')

		plt.axis([0,5,0,5])

	def plot_point(self, pt, mkr):
		pt = pt.reshape(2)
		plt.plot(pt[0], pt[1], mkr)

	def plot_show(self):
		plt.draw()

	def plot_reset(self):
		plt.clf()

	def plot_pause(self):
		plt.pause(0.05)
