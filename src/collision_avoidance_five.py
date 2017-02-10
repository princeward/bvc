import numpy as np
import random
import robot

COLORS = ['b', 'g', 'r', 'c', 'm', 'y']

# collison avoidance example with five robots
#corners = np.array([[0.0, 0.0, 5.0, 5.0], 
#					[0.0, 5.0, 5.0, 0.0]])
corners = np.array([[-1.0, -1.0, 6.0, 6.0], 
					[-1.0, 6.0, 6.0, -1.0]])
safe_rad = 0.2

robots = [robot.Robot(i) for i in range(5)]
robots[0].set_bvc([2.5 + 0.1*random.random(), 0.5 + 0.1*random.random()], safe_rad, corners)
robots[1].set_bvc([4.5 + 0.1*random.random(), 2.5 + 0.1*random.random()], safe_rad, corners)
robots[2].set_bvc([4.0 + 0.1*random.random(), 4.5 + 0.1*random.random()], safe_rad, corners)
robots[3].set_bvc([1.0 + 0.1*random.random(), 4.5 + 0.1*random.random()], safe_rad, corners)
robots[4].set_bvc([0.5 + 0.1*random.random(), 2.5 + 0.1*random.random()], safe_rad, corners)

robots[0].set_goal([2.5, 4.5])
robots[1].set_goal([0.5, 2.5])
robots[2].set_goal([0.5, 0.5])
robots[3].set_goal([4.5, 0.5])
robots[4].set_goal([4.5, 2.5])


aa = [0,1,2,3,4]
for loop in range(300):
	# print loop
	# all robots re-compute BVC
	for i in range(5):
		all_pos = np.zeros((2,5))
		for r in range(5):
			all_pos[:,r,None] = robots[r].get_pos() # none is used to preserve dimensionality

		bb = [x for x in aa if x != i]
		bb_array = np.array(bb) # this is the IDs of neighboring robots
		
		# extract own pos for robot i, as well as all other neighbor robots
		own_pos = all_pos[:,i]
		nbr_pos = all_pos[:,bb_array]
		
		robots[i].cell.update_bvc(own_pos, nbr_pos, bb_array)
		robots[i].mem_nbr_dist(nbr_pos, bb_array)
		robots[i].cell.plot_bvc()

	# all robots do geometric solution
	for rbt in robots:
		closest = rbt.bvc_find_closest_to_goal()
		rbt.move(closest)

		rbt.cell.plot_self('o', COLORS[rbt.id])
		rbt.plot_id()
		rbt.cell.plot_point(closest, 'o', 'y')
		rbt.cell.plot_point(rbt.goal, '*', COLORS[rbt.id], 10)

	robots[0].cell.plot_show()
	robots[0].cell.plot_pause(nsec = 0.01)
	robots[0].cell.plot_reset()

	# debug
	for rbt in robots:
		'''
		if rbt.id == 0:
			# check potential deadlock
			if rbt.check_potential_deadlock():
				rbt.find_critical_unstable_pt()
		'''
		pass

	#if loop % 10 == 0:
	#	raw_input("Press Enter to continue...")