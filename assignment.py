import matplotlib.pyplot as plt
import random
import numpy as np
import os
import time
import math

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
         i1, i2 = random.sample(idx, 2)
         seq[i1], seq[i2] = seq[i2], seq[i1]


def inversion(listToSwap, first, second):

    #listToSwap[idxStart], listToSwap[idxEnd] = listToSwap[idxEnd], listToSwap[idxStart]
    newList = listToSwap[first : second]
    newList.reverse()

    return newList



def evaluate(sol):
    total = 0
    for i in range(len(sol)-1):
        #it takes the real values of the list
        first, second = test_colours[sol[i]], test_colours[sol[i + 1]]
        result = sqrt((float(first[0]) - float(second[0]))**2 + (float(first[1]) - float(second[1]))**2 + (float(first[2]) - float(second[2]))**2)
        total += result

    return total


def random_neighbour(s):
        neighbour = deepcopy(s)
        list = [0,0]
        list[0] = random.randint(0,len(neighbour)-1)
        list[1] = random.randint(0,len(neighbour)-1)
        #this is to make sure that the first is going to be bigger then the second index
        list.sort()
        neighbour[list[0] : list[1] + 1] = inversion(neighbour, list[0], list[1] + 1)

        return neighbour

def hill_climbing(setOfColours, numberOfIterations):
    #shuffling the list of colours so that i can have a random start
    sol = random.sample(range(len(setOfColours)), len(setOfColours)) # list of colours for testing
    #list that will save all improvements from the eval_neighbour_sol to be then used for graph
    solListEval = []
    eval_sol = evaluate(sol)
    solListEval.append(eval_sol)
    while (numberOfIterations>0):
        neighbour_sol = random_neighbour(sol)
        eval_neighbour_sol = evaluate(neighbour_sol)
        if (eval_neighbour_sol < eval_sol):
            #print "eval_neighbour_sol e ", eval_neighbour_sol
            sol = deepcopy(neighbour_sol)
            eval_sol = eval_neighbour_sol

        numberOfIterations -=1
        solListEval.append(eval_sol)

    return sol, eval_sol, solListEval



def multi_hc(tries, setOfColours,hillclimbingIterations):

    bestSolIndexes = []
    bestSolEval = 0

    worstSolIndexes =[]
    worstSolEval = 0

    values = []

    for i in range(0, tries):
        solIndexes,solSum, val = hill_climbing(setOfColours,hillclimbingIterations)

        #values.append(solSum)

        if len(bestSolIndexes) < 1:
            bestSolIndexes = solIndexes[:]
            bestSolEval = solSum

            worstSolIndexes = solIndexes[:]
            worstSolEval = solSum
            continue;

        if (bestSolEval > solSum):
            bestSolIndexes = solIndexes[:]
            bestSolEval = solSum

        if (worstSolEval < solSum):
            worstSolIndexes = solIndexes[:]
            worstSolEval = solSum

    return bestSolIndexes, bestSolEval, worstSolIndexes, worstSolEval



def annealing(setOfColours,temperature,numberOfIterations):
    #shuffling the list of colours so that i can have a random start
    sol = random.sample(range(len(setOfColours)), len(setOfColours)) # list of colours for testing

    solListEval = []
    solListIndices = []
    eval_sol = evaluate(sol)

    tempmin = 0.00001
    alpha = 0.7
    solListEval.append(eval_sol)

    while (tempmin < temperature):
        for i in range(numberOfIterations):
            eval_neighbour_sol = evaluate(sol)
            deltaEval = (eval_neighbour_sol - eval_sol)
            if (eval_neighbour_sol < eval_sol):
                solListEval.append(eval_neighbour_sol)
                #print "eval_neighbour_sol e ", eval_neighbour_sol
                bestSolList = deepcopy(sol)
                eval_sol = eval_neighbour_sol
            else:
                if((1 - math.exp(-deltaEval/temperature)) > random.uniform(0,1)):
                    solListEval.append(eval_neighbour_sol)
                    #print "eval_neighbour_sol e ", eval_neighbour_sol
                    bestSolList = deepcopy(sol)
                    eval_sol = eval_neighbour_sol
        temperature = temperature * alpha



    return bestSolList, eval_sol, solListEval


#####_______main_____######

# Get the directory where the file is located
dir_path = os.path.dirname(os.path.realpath(__file__))
os.chdir(dir_path) # Change the working directory so we can read the file


ncolors, colours = read_file('colours.txt')  # Total number of colours and list of colours

test_size = 100 # Size of the subset of colours for testing
test_colours = colours[0:test_size]  # list of colours for testing

test_size500 = 500 # Size of the subset of colours for testing
test_colours500 = colours[0:test_size500]  # list of colours for testing

#permutation = random.sample(range(test_size), test_size) # produces random pemutation of lenght test_size, from the numbers 0 to test_size -1
#plot_colours(test_colours, permutation)


#Greedy algorithm starting with random colour for 100
'''
print "------ Greedy algorithm for 100 -----"
firstRandomColor = random.randint(0,test_size)
start = time.time()
greedyColoursList, sumOfGreedySol = constructive(test_colours,firstRandomColor)
end = time.time()
print "The time of execution for the greedy algotrithm with 100 colours is: ", end-start, ". The sum of the distances is: ", sumOfGreedySol
#plot_colours(test_colours, greedyColoursList)
'''

#this is to get the range for the sum of distance of Greedy algorithm solution for 100 colours
'''
resultCollectorGreedySol100 = []

while(len(resultCollectorGreedySol100)<100):
    firstRandomColor = random.randint(0,test_size)
    greedyColoursList, sumOfGreedySol = constructive(test_colours,firstRandomColor)
    resultCollectorGreedySol100.append(sumOfGreedySol)

minSolution = min(resultCollectorGreedySol100)
maxSolution = max(resultCollectorGreedySol100)

print "The approximate range of solution for distance with the greedy algorithm for 100 colours is between ", minSolution ," and ",maxSolution
'''

#Greedy algotrithm starting with random colour for 500
'''
print "------ Greedy algorithm for 500 -----"
firstRandomColor = random.randint(0,test_size500)
start = time.time()
greedyColoursList, sumOfGreedySol = constructive(test_colours500,firstRandomColor)
end = time.time()
print "The time of execution for the greedy algotrithm with 500 colours is: ", end-start, ". The sum of the distances is: ", sumOfGreedySol
#plot_colours(test_colours500, greedyColoursList)
'''

#this is to get the range for the sum of distance of Greedy algorithm solution for 500 colours
'''
resultCollectorGreedySol500 = []

while(len(resultCollectorGreedySol500)<100):
    firstRandomColor = random.randint(0,test_size500)
    greedyColoursList, sumOfGreedySol = constructive(test_colours500,firstRandomColor)
    resultCollectorGreedySol500.append(sumOfGreedySol)

print resultCollectorGreedySol500
minSolution = min(resultCollectorGreedySol500)
maxSolution = max(resultCollectorGreedySol500)

print "The approximate range of solution for distance with the greedy algorithm for 500 colours is between ", minSolution ," and ",maxSolution
'''

hillclimbingIterations = 10000
multiHCIterations = 20

#Starting with 100
start = time.time()
hillclimbingIndexes, sumofHillclimbingSol, valuesHillClimb = hill_climbing(test_colours, hillclimbingIterations)
end = time.time()
print "The time of execution for the hill-climbing with 200 iterations for 100 sample list is: ", end - start
print "The sum of distances for hill_climbing is: ", sumofHillclimbingSol
plot_colours(test_colours, hillclimbingIndexes)

#populating a list from 1 until the the lenght of valuesHillClimb
improvements = []
for iter in range(len(valuesHillClimb)):
    improvements.append(iter)

plt.figure()
plt.plot(improvements,valuesHillClimb)
plt.ylabel('values')
plt.xlabel('number of times improved')
plt.title("Improvements in Hill Climbing algotrithm for 100")
plt.show()
'''
#starting with 500
start = time.time()
hillclimbingIndexes, sumofHillclimbingSol, valuesHillClimb = hill_climbing(test_colours500, hillclimbingIterations)
end = time.time()
print "The time of execution for the hill-climbing with 200 iterations for 500 sample list is: ", end - start
print "The sum of distances for hill_climbing is: ", sumofHillclimbingSol
plot_colours(test_colours500, hillclimbingIndexes)

#populating a list from 1 until the the lenght of valuesHillClimb
improvements = []
for iter in range(len(valuesHillClimb)):
    improvements.append(iter)


plt.figure()
plt.plot(improvements,valuesHillClimb)
plt.ylabel('values')
plt.xlabel('number of times improved')
plt.title("Improvements in Hill Climbing algotrithm for 500")
plt.show()
'''
'''
#return values bestSolIndexes, bestSolEval,  (worstSolIndexes, worstSolEval, are only for the range of solutions)
#starting with 100
start = time.time()
bestSolIndexes, bestSolEval, worstSolIndexes,worstSolEval = multi_hc(multiHCIterations, test_colours, hillclimbingIterations)
end = time.time()

print "The time of the multi_hc is: ", end - start
print "The best solution's sum of distances: ", bestSolEval
print "The worst solution is: ", worstSolEval
#print "list of values", bestSolIndexes
#print "list of worst values", worstSolIndexes
#plot_colours(test_colours, bestSolIndexes)

#starting with 500

start = time.time()
bestSolIndexes, bestSolEval, worstSolIndexes, worstSolEval = multi_hc(multiHCIterations, test_colours500, hillclimbingIterations)
end = time.time()
print "The time of the multi_hc is: ", end - start
print "The best solution's sum of distances: ", bestSolEval
print "The worst solution is: ", worstSolEval
#print "list of values", bestSolIndexes
#plot_colours(test_colours500, bestSolIndexes)
'''
'''
#annealing for 100
temperatureAnnealing = 1
start = time.time()
annealingIndecies, sumOfAnnealingEval, valuesAnnealing = annealing(test_colours,temperatureAnnealing, hillclimbingIterations)
end = time.time()
print "The time of execution for the annealing with 200 iterations and temp of 100 for 100 sample list is: ", end - start
print "The sum of distances for simulated_annealing is: ", sumOfAnnealingEval
#print "Values gathered", valuesAnnealing
#plot_colours(test_colours, annealingIndecies)
'''
