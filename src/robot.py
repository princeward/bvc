import numpy as np
import geo_helper
import bvc
import random
import pdb
import matplotlib.pyplot as plt

class Robot:
	# static variable
	MAX_SPD = 0.05
	N_ROBOT_MAX = 10
	DEADLOCK_THRES_RATIO = 5.0 # threshold for checking deadlock (3-agent case): DEADLOCK_THRES_RATIO * safety_radius
								# this parameter must be larger than 2

	def __init__(self, id = None, pos = None, goal = None):
		if pos is not None:
			self.pos = np.reshape(pos,(2,1))
		
		if id is not None:
			self.id = id
		
		if goal is not None:
			self.goal = goal

		self.old_nbr_dist = -np.ones(self.N_ROBOT_MAX) # for tracking the distance with neighbors, pre-allocate space, -1 means not tracked
		self.new_nbr_dist = -np.ones(self.N_ROBOT_MAX)
		self.action_nbr_idx = None

	def set_bvc(self, own_pos, safe_rad, world_corners):
		self.pos = np.reshape(own_pos,(2,1))
		self.cell = bvc.BVC(self.pos, safe_rad, world_corners)

	def set_goal(self, goal):
		self.goal = np.reshape(goal, (2,1))

	def get_pos(self):
		return self.pos

	def mem_nbr_dist(self, nbr_pos, nbr_id):
		'''
		memorize the distance to neighbors
		assumed the ordering of the nbr_pos is preserved every time
		'''

		self.nbr_pos = nbr_pos

		if self.action_nbr_idx is not None:
			return

		for m in range(self.N_ROBOT_MAX):
			if m not in self.cell.valid_nbr_id: # m is not a valid neighbor, don't track
				self.old_nbr_dist[m] = -1
				self.new_nbr_dist[m] = -1
			else:
				if self.new_nbr_dist[m] != -1:
					self.old_nbr_dist[m] = self.new_nbr_dist[m]
					idd = np.where(nbr_id == m)
					dist = np.linalg.norm(nbr_pos[:,idd].reshape(2,1) - self.pos)
					self.new_nbr_dist[m] = dist
					#pdb.set_trace()
					if self.new_nbr_dist[m] - self.old_nbr_dist[m] < np.minimum(- 1.3 * self.MAX_SPD, -0.1*self.new_nbr_dist[m]):
						self.action_nbr_idx = idd # TODO: right now only take care of one such neighbor
						self.action_nbr_range = self.old_nbr_dist[m]
						# print "dectected"
				else:
					idd = np.where(nbr_id == m)
					dist = np.linalg.norm(nbr_pos[:,idd].reshape(2,1) - self.pos)
					self.new_nbr_dist[m] = dist


	def bvc_find_closest_to_goal(self, goal = None):
		if goal is not None:
			self.goal = goal

		if self.check_potential_deadlock(): # robot is close to potential deadlock
			return self.find_critical_unstable_pt() # Note: may return None
		else: # no potential deadlock, return regular closest point
			return self.cell.find_closest_to_goal(self.goal)
		
		'''
		# precaution action, old code
		if self.action_nbr_idx is not None:
			if geo_helper.is_in_circle(self.closest, self.nbr_pos[:,self.action_nbr_idx], self.action_nbr_range):
				# choose a vertex to detour along
				dist_min = float('inf')
				for i in range(self.cell.N):
					if not geo_helper.is_in_circle(self.cell.ver[:,i], self.nbr_pos[:,self.action_nbr_idx], self.action_nbr_range):
						#self.closest = self.cell.ver[:,i]
						candidate = self.cell.ver[:,i]
						dist = np.linalg.norm(candidate-self.goal)
						if dist < dist_min:
							self.closest = candidate
			else:
				self.action_nbr_idx = None
		'''

	def check_potential_deadlock(self):
		'''
		Check whether a potential deadlock exists (three-robot case)
		Input:
			None
		Return:
			True: has potential deadlock
			False: no deadlock
		'''
		self.close_nbr_pos = None # a list of close neighbors that may cause potential deadlock
										   # two consecutive columns are pairs

		# get valid neighbors (those have a shared Voronoi edge)
		valid_nbr_pos, valid_nbr_id = self.cell.get_valid_nbr()
		self.valid_nbr_pos = valid_nbr_pos
		self.valid_nbr_id = valid_nbr_id
		
		check_result = False
		n_nbr = len(valid_nbr_id)
		if n_nbr < 2:
			return False
		else:
			# calculate distances
			dist1 = np.zeros(n_nbr+1) # dist between two neighbors
			dist2 = np.zeros(n_nbr) # dist between self and neighbor 1
			dist3 = np.zeros(n_nbr) # dist between self and neighbor 2

			for i in range(n_nbr-1): # check pairs, note neighbors are sorted by angle
				dist1[i] = geo_helper.get_two_point_dist(valid_nbr_pos[:,i], valid_nbr_pos[:, i+1]) 
				dist2[i] = geo_helper.get_two_point_dist(valid_nbr_pos[:,i], self.pos) # dist between 
				dist3[i] = geo_helper.get_two_point_dist(valid_nbr_pos[:,i+1], self.pos)
			# distance between last neighbor and first neighbor
			# when there are only two neighbors, this is also true, but redundant computation
			dist1[n_nbr-1] = geo_helper.get_two_point_dist(valid_nbr_pos[:,n_nbr-1], valid_nbr_pos[:, 0])

			# check potential deadlock
			for i in range(n_nbr-1):
				if (dist1[i] <= self.DEADLOCK_THRES_RATIO * self.cell.safe_rad) and \
				(dist2[i] <= self.DEADLOCK_THRES_RATIO * self.cell.safe_rad) and \
				(dist3[i] <= self.DEADLOCK_THRES_RATIO * self.cell.safe_rad):
					#print "potential deadlock for robot {0}".format(self.id)
					# TODO: Check whether the closest point is on the vertice of the BVC
					if self.close_nbr_pos is None:
						self.close_nbr_pos = valid_nbr_pos[:, i].reshape((2,1))
					else:
						self.close_nbr_pos = np.hstack( (self.close_nbr_pos, valid_nbr_pos[:, i].reshape((2,1))) )
					self.close_nbr_pos = np.hstack( (self.close_nbr_pos, valid_nbr_pos[:, i+1].reshape( (2,1) )) )
					check_result = True
			# check tail and head case
			if (dist1[n_nbr-1] <= self.DEADLOCK_THRES_RATIO * self.cell.safe_rad) and \
			(dist2[n_nbr-1] <= self.DEADLOCK_THRES_RATIO * self.cell.safe_rad) and \
			(dist3[n_nbr-1] <= self.DEADLOCK_THRES_RATIO * self.cell.safe_rad):
				#print "potential deadlock for robot {0}".format(self.id)
				# TODO: Check whether the closest point is on the vertice of the BVC
				if self.close_nbr_pos is None:
					self.close_nbr_pos = valid_nbr_pos[:, n_nbr-1].reshape( (2,1) )
				else:
					self.close_nbr_pos = np.hstack( (self.close_nbr_pos, valid_nbr_pos[:, n_nbr-1].reshape( (2,1) ) ) )
				self.close_nbr_pos = np.hstack( (self.close_nbr_pos, valid_nbr_pos[:, 0].reshape( (2,1) ) ) )
				check_result = True

		# no potential deadlock detected
		return check_result

	def find_critical_unstable_pt(self):
		'''
		When check_potential_deadlock() returns true, find a better point to go (a unstable position in two-robot case)
		The unstable point makes sure that we do not end up in a stable deadlock situation
		Important:
			Must call check_potential_deadlock() before to update close neighbors
		'''

		dist_min = float('Inf')
		result_pt = None

		# iterate through all valid neighbors
		for i in range(len(self.valid_nbr_id)): # assume self.valid_nbr_id has been updated in check_potential_deadlock()
			# first, make sure the neighbor robot "blocks" the way to goal
			nbr = self.valid_nbr_pos[:, i].reshape((2,1))
			v1 = self.goal - self.pos
			v2 =  nbr - self.pos
			if np.dot(v1.T, v2) > 0: # this neighbor is indeed blocking the way
				line = geo_helper.get_line(self.goal, nbr) # get a line between goal and this neighbor
				inter_pts = self.cell.get_line_intersect_bvc(line) # may have 0, 1, 2 intersection point(s)
				# first, check if there is no intersection point
				if inter_pts is None:
					#break
					inter_pts = self.cell.ver
				# consider all the intersection points
				for i in range(inter_pts.shape[1]):
					# second, check if the point is too close to deadlock situation
					for j in range(self.close_nbr_pos.shape[1]/2):
						if (geo_helper.get_two_point_dist(inter_pts[:,i], self.close_nbr_pos[:, 2*j]) < self.DEADLOCK_THRES_RATIO * self.cell.safe_rad) or \
						(geo_helper.get_two_point_dist(inter_pts[:,i], self.close_nbr_pos[:, 2*j+1]) < self.DEADLOCK_THRES_RATIO * self.cell.safe_rad):
							break
						# lastly, compare which one is a better point
						dist = geo_helper.get_two_point_dist(inter_pts[:, i], self.goal)
						if (dist < dist_min) and (dist > self.DEADLOCK_THRES_RATIO * self.cell.safe_rad):
							dist_min = dist
							result_pt = inter_pts[:, i].reshape((2,1))
					# TODO: there is a chance that intersection point is outside of BVC, but there are good points to find

		return result_pt



	def move(self, set_pt):
		'''
		move to a setpoint under speed limit
		update its own position
		'''
		if set_pt is None: # if set_pt is none, do not move
			print "set_pt is None for robot {0}".format(self.id)
			return

		set_pt = np.reshape(set_pt, (2,1))
		dist = np.linalg.norm(set_pt - self.pos)
		unit_dir = (set_pt - self.pos) / dist
		if dist < self.MAX_SPD:
			self.pos = set_pt - 0.02 * unit_dir # for stability reason
		else:
			self.pos += self.MAX_SPD * unit_dir

	def plot_id(self):
		plt.text(self.pos[0], self.pos[1], str(self.id))



