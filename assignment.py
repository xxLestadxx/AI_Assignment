import matplotlib.pyplot as plt
import numpy as np
import random
import os
from math import sqrt
from random import shuffle
from copy import deepcopy


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

def constructive(setOfColours,randomColour):
    constructive_order = []
    res = []
    sumOfSolutions = 0
    colours_copy = setOfColours[0:len(setOfColours)]
    n = randomColour
    # where n is the first random color and then is the next color
    constructive_order.append(setOfColours.index(colours_copy[n]))
    #print "index ", setOfColours.index(colours_copy[n])
    #print "first color ", colours_copy[n]
    del colours_copy[n]

    while(len(constructive_order)<len(setOfColours)):
#        print "n e ", n
        for p in range(len(colours_copy)):
            #print "setOfColours color", setOfColours[n]
            #print "colours_copy color ", colours_copy[p]
            result = sqrt((float(setOfColours[n][0]) - float(colours_copy[p][0]))**2 + (float(setOfColours[n][1]) - float(colours_copy[p][1]))**2 + (float(setOfColours[n][2]) - float(colours_copy[p][2]))**2)
            res.append(result)
            best = min(res)
            i = res.index(best)

        #print "best is ", best
    #    print "res ot ",i, " e to ", res[i]
    #    print "colors copy ot i e ", colours_copy[i]
#        print "color e ", colours.index(colours_copy[i])
        sumOfSolutions += best
        n = setOfColours.index(colours_copy[i])
#        print "n e ", n
#        print "colors original ot i e " ,colours[i]
        res = []
        constructive_order.append(setOfColours.index(colours_copy[i]))
        del colours_copy[i]

    return constructive_order, sumOfSolutions



def listOfIndexesIntoColours(listOfIndexes):
    list = []
    for i in listOfIndexes:
        list.append(colours[i])

    return list

def swap_random(seq):
         idx = range(len(seq))
        # print "lengtha e tolkoz", len(seq)
         i1, i2 = random.sample(idx, 2)
         seq[i1], seq[i2] = seq[i2], seq[i1]


def evaluate(sol):
    total = 0
    for i in range(len(sol)-1):
        result = sqrt((float(sol[i][0]) - float(sol[i+1][0]))**2 + (float(sol[i][1]) - float(sol[i+1][1]))**2 + (float(sol[i][2]) - float(sol[i+1][2]))**2)
        total += result

    return total


def hill_climbing(setOfColours, numberOfIterations):
    #shuffling the list of colours so that i can have a random start
    sol = random.sample(setOfColours, len(setOfColours)) # list of colours for testing
    bestSolList = [[]]
    solListIndexes = []
    eval_sol = evaluate(sol)
    while (numberOfIterations>0):
        #This function swaps two of the indexes does not return anything, as the change is done to the list itself
        swap_random(sol)
        eval_neighbour_sol = evaluate(sol)
        if (eval_neighbour_sol < eval_sol):
            #print "eval_neighbour_sol e ", eval_neighbour_sol
            bestSolList = deepcopy(sol)
            eval_sol = eval_neighbour_sol

        numberOfIterations -=1

    for i in range(len(bestSolList)):
        solListIndexes.append(colours.index(bestSolList[i]))
#    print " tui to " ,solListIndencies

    return solListIndexes, eval_sol




#to see if two lists are with the same color
def equal_ignore_order(a, b):
    unmatched = list(b)
    for element in a:
        try:
            unmatched.remove(element)
        except ValueError:
            return False
    return not unmatched


def multi_hc(tries):
    bestSolIndexes = []
    bestSolEval = 0

    values = []

    for i in range(0, tries):
        solIndexes,solSum = hill_climbing(test_colours,hillclimbingIterations)

        values.append(solSum)

        if len(bestSolIndexes) < 1:
            bestSolIndexes = solIndexes[:]
            bestSolEval = solSum
            continue;

        if (bestSolEval > solSum):
            bestSolIndexes = solIndexes
            bestSolEval = solSum

    return values, bestSolIndexes, bestSolEval


#####_______main_____######

# Get the directory where the file is located
dir_path = os.path.dirname(os.path.realpath(__file__))
os.chdir(dir_path) # Change the working directory so we can read the file


ncolors, colours = read_file('colours.txt')  # Total number of colours and list of colours

test_size = 100 # Size of the subset of colours for testing
test_colours = colours[0:test_size]  # list of colours for testing

permutation = random.sample(range(test_size), test_size) # produces random pemutation of lenght test_size, from the numbers 0 to test_size -1
#plot_colours(test_colours, permutation)

firstRandomColor = random.randint(0,test_size)
greedyColoursList, sumOfGreedySol = constructive(test_colours,firstRandomColor)
#print "the sum of the distances in greedy is : ", sumOfGreedySol
#plot_colours(test_colours, greedyColoursList)
#print "greedy is: ", greedyColoursList
permutationColours = listOfIndexesIntoColours(greedyColoursList)

b = equal_ignore_order(test_colours, permutationColours)
print b
#greedyColoursList, sumOfGreedySol = constructive(permutationColours,firstRandomColor)
#print "the sum of the distances in greedy permutation is: ", sumOfGreedySol
#plot_colours(test_colours,greedyColoursList)

hillclimbingIterations = 200
multiHCIterations = 20

hillclimbingIndexes, sumofHillclimbingSol = hill_climbing(test_colours, hillclimbingIterations)
#print " the sum in hill_climbing is : ", sumofHillclimbingSol
#plot_colours(test_colours, hillclimbingIndexes)



#return values, bestSolIndexes, bestSolEval
v, b, s = multi_hc(multiHCIterations)
print "list of values", v
print "best value", s
plot_colours(test_colours, b)
