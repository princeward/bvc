import numpy as np
import pdb

def is_on_line(pt, line, tol = 1e-3):
	'''
	check whether a point is on the line
	'''
	line = line.reshape(1,3)
	if np.absolute( np.dot(line[0,0:2], pt) + line[0,2] ) <= tol:
		return True
	else:
		return False

def is_in_half_plane(pt, line, tol = 1e-3):
	'''
	line: 1x3 array
	tol: tolerance for numerical reason
	'''
	line = line.reshape(1,3)
	if np.dot(line[0,0:2], pt) + line[0,2] <= tol:
		return True
	else:
		return False

def is_in_hull(pt, lines, tol = 1e-3):
	'''
	determine if a point is inside a convex hull formed by lines
	lines: nx3 matrix
	'''
	nline, garbage = lines.shape
	for i in range(nline):
		if not is_in_half_plane(pt, lines[i,:]):
			return False
	return True

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
		if(np.dot(n, inner_pt) + a > 0):
			line = -line
		if(np.dot(n, inner_pt) + a == 0):
			print("[Warn]: cannot decide half plane in get_line()")

	return line

def get_one_bvc_edge(pi, po, safe_rad):
	'''
	pi - inner point, po - outer point, safe_rad - safety radius
	Returns parameters for line [a b]*q + c = 0 (buffered voronoi edge)
	'''
	pi = pi.reshape(2,1)
	po = po.reshape(2,1)

	line = np.zeros([1,3])
	dist = np.linalg.norm(po-pi)
	if(dist/2.0 < safe_rad):
		print("[Err]: Distance is smaller than 2xSafetyRadius.")
		return line

	unit_dir = (po-pi)/dist
	break_pt = pi + unit_dir*(dist/2-safe_rad) # the retracted mid-point
	c = - np.dot(unit_dir.T, break_pt)

	line[0,0:2] = unit_dir.T
	line[0,2] = c
	
	if is_in_half_plane(pi, line):
		return line
	else:
		return -line

def get_line_intersect(line1, line2):
	# Solve intersection point: Ax + c = 0, x = - A^-1 c
	line1 = line1.reshape(1,3)
	line2 = line2.reshape(1,3)

	A = np.zeros([2,2])
	c = np.zeros([2,1])

	A[0,:] = line1[0,0:2]
	c[0,0] = line1[0,2]
	A[1,:] = line2[0,0:2]
	c[1,0] = line2[0,2]
	Ainv = np.linalg.inv(A) # TODO: A not invertible?

	return (-np.dot(Ainv, c) ) # note there should be a neg sign




