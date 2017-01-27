import numpy as np
import random
import robot

# collison avoidance example with five robots
corners = np.array([[0.0, 0.0, 5.0, 5.0], 
					[0.0, 5.0, 5.0, 0.0]])
safe_rad = 0.2

robots = [robot.Robot(i) for i in range(5)]
robots[0].set_bvc([2.5 + 0.2*random.random(), 0.5 + 0.2*random.random()], safe_rad, corners)
robots[1].set_bvc([4.5 + 0.2*random.random(), 2.5 + 0.2*random.random()], safe_rad, corners)
robots[2].set_bvc([4.0 + 0.2*random.random(), 4.5 + 0.2*random.random()], safe_rad, corners)
robots[3].set_bvc([1.0 + 0.2*random.random(), 4.5 + 0.2*random.random()], safe_rad, corners)
robots[4].set_bvc([0.5 + 0.2*random.random(), 2.5 + 0.2*random.random()], safe_rad, corners)

robots[0].set_goal([2.5, 4.5])
robots[1].set_goal([0.5, 2.5 - 0.3])
robots[2].set_goal([0.5, 0.5])
#robots[2].set_goal([3.0, 2.0])
robots[3].set_goal([4.5, 0.5])
robots[4].set_goal([4.5, 2.5+0.2])


aa = [0,1,2,3,4]
for loop in range(200):
	print loop
	# all robots re-compute BVC
	for i in range(5):
		all_pos = np.zeros((2,5))
		for r in range(5):
			all_pos[:,r,None] = robots[r].get_pos() # none is used to preserve dimensionality

		bb = [x for x in aa if x != i]
		bb_array = np.array(bb)
		
		# extract own pos for robot i, as well as all other neighbor robots
		own_pos = all_pos[:,i]
		nbr_pos = all_pos[:,bb_array]
		
		robots[i].cell.update_bvc(own_pos, nbr_pos)
		robots[i].cell.plot_bvc()

	# all robots do geometric solution
	for rbt in robots:
		rbt.bvc_find_closest_to_goal()
		rbt.plot_point(rbt.closest, 'ro')
		rbt.plot_point(rbt.goal, 'go')

	robots[0].cell.plot_show()
	robots[0].cell.plot_pause()
	robots[0].cell.plot_reset()

	# all robots move
	for rbt in robots:
		rbt.move(rbt.closest)