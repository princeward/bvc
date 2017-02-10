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

		# ver and lines are BVC vertices and edges
		self.ver = np.empty([2, self.NUM_VER_PRE_ALLOCATE]) # pre-allocate 20 vertices for the efficiency of pending
		self.lines = np.empty([self.NUM_VER_PRE_ALLOCATE, 3]) # pre-allocate

	def get_valid_nbr(self):
		'''
		Valid neighbor is the neighbor who share an Voronoi edge 
		return is ordered (in terms of arctan2 angle)
		'''
		# sort based on arctan2: [-pi,pi]
		n_nbr = len(self.valid_nbr_id)
		angles = np.zeros(n_nbr)
		for i in range(n_nbr):
			angles[i] = np.arctan2(self.valid_nbr_pos[1,i]-self.own_pos[1,0],  self.valid_nbr_pos[0,i]-self.own_pos[0,0])
		angIdx = np.argsort(angles)

		# sort position
		self.valid_nbr_pos = self.valid_nbr_pos[:, angIdx]

		# sort id, since self.valid_nbr_id is not a numpy array, have to do the for loop
		temp = self.valid_nbr_id
		self.valid_nbr_id = []
		for i in range(n_nbr):
			self.valid_nbr_id.append( temp[angIdx[i]] )

		return (self.valid_nbr_pos, self.valid_nbr_id) # the columns in valid_nbr_pos and valid_nbr_id correspond to the same rbt

	def update_bvc(self, own_pos, nbr_pos, nbr_id):
		# FACT: voronoi cell always has the same number of vertices and edges
		self.valid_nbr_id = [] # neighbors that create a BVC edge 
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
		nbr_id = nbr_id[idx]

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
				self.valid_nbr_id.append(nbr_id[i]) # memorize the id of this valid neighbor
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

		# update the valid neighbor index/position
		self.valid_nbr_pos = np.zeros( (2,len(self.valid_nbr_id)) )
		for itr, valid_id in enumerate(self.valid_nbr_id):
			col = np.where(nbr_id == valid_id)
			self.valid_nbr_pos[:, itr] = nbr_pos[:, col[0][0]]

		# truncate self.lines and self.ver to get rid of the useless pre-allocated empty elements
		self.lines = self.lines[0:self.N, :]
		self.ver = self.ver[:, 0:self.N]

	def find_closest_to_goal(self, goal, pre_pt_min=None, dist_penalty_ratio = 0.05):
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
			if (pre_pt_min is not None) and (pre_pt_min.size != 0): # penalize on jumping too far from previous
				dist += dist_penalty_ratio * geo_helper.get_two_point_dist(self.ver[:,i], pre_pt_min)
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
				if (pre_pt_min is not None) and (pre_pt_min.size != 0): # penalize on jumping too far from previous
					dist += dist_penalty_ratio * geo_helper.get_two_point_dist(proj_pt, pre_pt_min)
				if dist < dist_min:
					dist_min = dist
					pt_min = proj_pt

		return pt_min

	def get_line_intersect_bvc(self, target_line):
		'''
		Given a line, find the intersection point with the BVC
		There may be no solution, one point, or two solutions
		INPUT:
			target_line: 3x1 line equation: ax + by + c = 0
		'''
		result_pt = None

		for i in range(self.lines.shape[0]): # consider all voronoi edges
			line = self.lines[i,:].reshape((1,3))
			inter_pt = geo_helper.get_line_intersect(line, target_line)
			if (inter_pt is not None) and geo_helper.is_in_hull(inter_pt, self.lines): # the intersection point is in the BVC
				if result_pt is None:
					result_pt = inter_pt.reshape((2,1))
				else:
					result_pt = np.hstack( (result_pt, inter_pt.reshape((2,1))) )

		return result_pt


	def plot_bvc(self):
		# draw own position
		#plt.plot(self.own_pos[0], self.own_pos[1],'bo', markersize = int(140*self.safe_rad))
		# draw neighbors
		#plt.plot(self.nbr_pos[0,:], self.nbr_pos[1,:], 'bo', markersize = int(140*self.safe_rad))
		# draw BVC
		plt.plot(self.ver[0, 0:self.N], self.ver[1, 0:self.N], 'k')
		plt.plot([self.ver[0,self.N-1], self.ver[0,0]], [self.ver[1,self.N-1], self.ver[1,0]], 'k')

		plt.axis([-1.1,6.1,-1.1,6.1])

	def plot_self(self, marker = 'o', clr = 'b'):
		plt.plot(self.own_pos[0], self.own_pos[1], marker, markersize = int(120*self.safe_rad), alpha = 0.8)

	def plot_point(self, pt, mkr = 'o', clr = 'b', ms = 5):
		if pt is not None:
			pt = pt.reshape(2)
			plt.plot(pt[0], pt[1], mkr, color = clr, markersize = ms)

	def plot_show(self):
		plt.draw()

	def plot_reset(self):
		plt.clf()

	def plot_pause(self, nsec):
		plt.pause(nsec)
