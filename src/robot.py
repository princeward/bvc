import numpy as np
import geo_helper
import bvc
import random

class Robot:
	# static variable
	MAX_SPD = 0.05

	def __init__(self, id = None, pos = None, goal = None):
		if pos is not None:
			self.pos = np.reshape(pos,(2,1))
		
		if id is not None:
			self.id = id
		
		if goal is not None:
			self.goal = goal

	def set_bvc(self, own_pos, safe_rad, world_corners):
		self.pos = np.reshape(own_pos,(2,1))
		self.cell = bvc.BVC(self.pos, safe_rad, world_corners)

	def set_goal(self, goal):
		self.goal = np.reshape(goal, (2,1))

	def get_pos(self):
		return self.pos

	def bvc_find_closest_to_goal(self, goal = None):
		if goal is not None:
			self.goal = goal
		self.closest = self.cell.find_closest_to_goal(self.goal)

	def plot_point(self, pt, mkr):
		self.cell.plot_point(pt, mkr)

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



