import numpy as np

def get_line(x1, x2, inner_pt = None):
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