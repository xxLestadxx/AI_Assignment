import matplotlib.pyplot as plt
import numpy as np
import random
import os
from math import sqrt
from random import shuffle

# Reads the file  of colours
# Returns the number of colours in the file and a list with the colours (RGB) values

def read_file(fname):
    with open(fname, 'r') as afile:
        lines = afile.readlines()
    n = int(lines[3])    # number of colours  in the file
    col = []
    lines = lines[4:]    # colors as rgb values
    for l in lines:
        rgb = l.split()
        col.append(rgb)
    return n, col

# Display the colours in the order of the permutation in a pyplot window
# Input, list of colours, and ordering  of colours.
# They need to be of the same length

def plot_colours(col, perm):

	assert len(col) == len(perm)

	ratio = 10 # ratio of line height/width, e.g. colour lines will have height 10 and width 1
	img = np.zeros((ratio, len(col), 3))
	for i in range(0, len(col)):
		img[:, i, :] = colours[perm[i]]

	fig, axes = plt.subplots(1, figsize=(8,4)) # figsize=(width,height) handles window dimensions
	axes.imshow(img, interpolation='nearest')
	axes.axis('off')
	plt.show()


#####_______main_____######

# Get the directory where the file is located
dir_path = os.path.dirname(os.path.realpath(__file__))
os.chdir(dir_path) # Change the working directory so we can read the file


ncolors, colours = read_file('colours.txt')  # Total number of colours and list of colours

test_size = 100  # Size of the subset of colours for testing
test_colours = colours[0:test_size]  # list of colours for testing

#permutation = random.sample(range(test_size), test_size) # produces random pemutation of lenght test_size, from the numbers 0 to test_size -1
#plot_colours(test_colours, permutation)

def constructive(numberOfColours,randomColour):
    constructive_order = []
    res = []
    sumOfSolutions = 0
    colours_copy = colours[0:numberOfColours]
    n = randomColour
    # where n is the first random color and then is the next color
    constructive_order.append(colours.index(colours_copy[n]))
    del colours_copy[n]
    while(len(constructive_order)<numberOfColours):
#        print "n e ", n
        for p in range(len(colours_copy)):
            result = sqrt((float(colours[n][0]) - float(colours_copy[p][0]))**2 + (float(colours[n][1]) - float(colours_copy[p][1]))**2 + (float(colours[n][2]) - float(colours_copy[p][2]))**2)
            res.append(result)
            best = min(res)
            i = res.index(best)

#        print "best is ", best
#        print "res ot ",i, " e to ", res[i]
#        print "colors copy ot i e ", colours_copy[i]
#        print "color e ", colours.index(colours_copy[i])
        sumOfSolutions += best
        n = colours.index(colours_copy[i])
#        print "n e ", n
#        print "colors original ot i e " ,colours[i]
        res = []
        constructive_order.append(colours.index(colours_copy[i]))
        del colours_copy[i]

    return constructive_order, sumOfSolutions

n = random.randint(0,test_size)
greedyColoursList, sumOfGreedySol = constructive(test_size,n)

#print "the list of the indexes is: " ,ko
print "the sum of the distances is", sumOfGreedySol
#plot_colours(test_colours, greedyColoursList)

#shuffling the list of colours so that I can different solutions
shuffledColours = random.sample(test_colours, len(test_colours))
#print "The shuffled is ", shuffledColours

def evaluate(sol):
    total = 0
    for i in range(len(sol)-1):
        result = sqrt((float(sol[i][0]) - float(sol[i+1][0]))**2 + (float(sol[i][1]) - float(sol[i+1][1]))**2 + (float(sol[i][2]) - float(sol[i+1][2]))**2)
        total += result

    return total

randomSol = evaluate(shuffledColours)
#print "random e tukaaa ", randomSol
