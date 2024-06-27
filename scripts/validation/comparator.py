import sys
import math
import logging

def absDiff( one, two , ignoreSign = False):
	if ignoreSign:
		one, two = math.fabs(one), math.fabs(two)
	return max(one, two) - min(one, two)
	
def compare( fExpected, fNew, epsilon, scalar=1.0, ignoreSign=False):
	fone = open( fExpected, 'r' )
	ftwo = open( fNew, 'r' )

	epsilon = float( epsilon )

	diffs = 0
	
	i = 0
	for lineone in fone:
		linetwo = ftwo.readline()

		elementsone = lineone.split()
		elementstwo = linetwo.split()

		i = i + 1

		lenone = len( elementsone )
		lentwo = len( elementstwo )

		if lenone != lentwo:
			diffs = diffs + 1
			logging.debug( "Line: %d differs in size." % ( i ) )
			logging.debug( "Should be %d elements but there are %d." % ( lenone, lentwo ) )
		else:
			for j in range( lenone ):
				# used to ignore XYZ meta data such as # of frames and # of atoms
				try:
					feone = int(elementsone[j])
					fetwo = int(elementstwo[j])
					# these are ints -- skip this pair
					continue
				except ValueError:
					pass
				try:
					feone = float( elementsone[j] ) * scalar
					fetwo = float( elementstwo[j] )

					fediff = absDiff( feone, fetwo, ignoreSign)

					if fediff > epsilon:
						diffs = diffs + 1
						logging.debug( "Line %d, Element %d Differs" % ( i, j ) )
						logging.debug( "Expected: %f, Actual: %f, Difference: %f" % ( feone, fetwo, fediff ) )
				except ValueError:
					pass
	return diffs == 0


