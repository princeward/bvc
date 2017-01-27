import numpy as np
import geo_helper
import bvc
import random
import pdb

class Robot:
	# static variable
	MAX_SPD = 0.05
	N_ROBOT_MAX = 10

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
		self.closest = self.cell.find_closest_to_goal(self.goal)
		
		
		# precaution action
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


	def move(self, set_pt):
		'''
		move to a setpoint under speed limit
		update its own position
		'''
		set_pt = np.reshape(set_pt, (2,1))
		dist = np.linalg.norm(set_pt - self.pos)
		unit_dir = (set_pt - self.pos) / dist
		if dist < self.MAX_SPD:
			self.pos = set_pt - 0.02 * unit_dir # for stability reason
		else:
			self.pos += self.MAX_SPD * unit_dir



