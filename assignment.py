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

'''
This is the Greedy algorithm
@param setOfColours takes the set of colours
@param randomColour takes a random number to indicate the first random colour
'''
def constructive(setOfColours,randomColour):
    #initializing of variables and lists
    constructive_order = []
    res = []
    sumOfSolutions = 0
    # where n is the first random color and then is the next color
    n = randomColour

    #creating a copy of the original list, to keep track of used colours
    colours_copy = setOfColours[0:len(setOfColours)]

    #add the first colour to the list and deletes it from the copy list
    constructive_order.append(setOfColours.index(colours_copy[n]))
    del colours_copy[n]

    while(len(constructive_order)<len(setOfColours)):
        for p in range(len(colours_copy)):
            result = sqrt((float(setOfColours[n][0]) - float(colours_copy[p][0]))**2 + (float(setOfColours[n][1]) - float(colours_copy[p][1]))**2 + (float(setOfColours[n][2]) - float(colours_copy[p][2]))**2)
            res.append(result)
            best = min(res)
            i = res.index(best)

        sumOfSolutions += best
        n = setOfColours.index(colours_copy[i])
        #resetting the result list so that there would be updated correct values
        res = []
        #adding the colour with minimum distance and deleting it from the copy list
        constructive_order.append(setOfColours.index(colours_copy[i]))
        del colours_copy[i]

    return constructive_order, sumOfSolutions


'''
Inverts the two colours in the list
taking as parameters the listToSwap, and the two colours
'''
def inversion(listToSwap, first, second):

    newList = listToSwap[first : second]
    newList.reverse()

    return newList


'''
Evaluates sum of distances of the provided list of colours
Returns the total sum of distances
'''
def evaluate(sol):
    total = 0
    for i in range(len(sol)-1):
        #it takes the real values of the list
        first, second = test_colours[sol[i]], test_colours[sol[i + 1]]
        result = sqrt((float(first[0]) - float(second[0]))**2 + (float(first[1]) - float(second[1]))**2 + (float(first[2]) - float(second[2]))**2)
        total += result

    return total

'''
A method to create a new random neighbour, by taking the list and then
inverting randomly two indexes
'''
def random_neighbour(s):
        neighbour = deepcopy(s)
        list = [0,0]
        list[0] = random.randint(0,len(neighbour)-1)
        list[1] = random.randint(0,len(neighbour)-1)
        #sorting is to make sure that the first is going to be bigger then the second index and that everything is going to be alright
        list.sort()
        neighbour[list[0] : list[1] + 1] = inversion(neighbour, list[0], list[1] + 1)

        return neighbour
'''
Simple implementation of the Hill Climbing algorithm
Takes setOfColours and numberOfIterations
Returns the best solution, the best solution's sum of distance and the list with all evaluations from the execution
'''
def hill_climbing(setOfColours, numberOfIterations):
    #shuffling the list of colours so that i can have a random start
    sol = random.sample(range(len(setOfColours)), len(setOfColours)) # list of colours for testing
    #list that will save all improvements from the eval_neighbour_sol to be then used for graph
    solListEval = []
    #evaluating the first solution and adding it to the list with solutions
    eval_sol = evaluate(sol)
    solListEval.append(eval_sol)

    while (numberOfIterations>0):
        #taking the first/best solution and inverting 2 indexes to create new neighbour then evaluating it
        neighbour_sol = random_neighbour(sol)
        eval_neighbour_sol = evaluate(neighbour_sol)
        #if the new solution is better, then the new solution becomes the best one
        if (eval_neighbour_sol < eval_sol):
            #print "eval_neighbour_sol e ", eval_neighbour_sol
            sol = deepcopy(neighbour_sol)
            eval_sol = eval_neighbour_sol

        numberOfIterations -=1
        #at the end appending a value, to keep track for the graph
        solListEval.append(eval_sol)

    return sol, eval_sol, solListEval

'''
Multi hill climbing algotrithm, runs the hill climbing for multiple times, takes the best solution at the end
Takes: tries - the number of iterations of the HC algorithm, setOfColours and hillclimbingIterations
Returns: best solution, the sum of distances of the best, list with evaluations of all the solutions
'''

def multi_hc(tries, setOfColours,hillclimbingIterations):

    bestSolIndecies = []
    bestSolEval = 0
    valuesMulti = []

    for i in range(0, tries):
        solIndecies,solSum, val = hill_climbing(setOfColours,hillclimbingIterations)
        valuesMulti.append(solSum)
        if len(bestSolIndecies) < 1:
            bestSolIndecies = solIndecies[:]
            bestSolEval = solSum
            continue;
        if (bestSolEval > solSum):
            bestSolIndecies = solIndecies[:]
            bestSolEval = solSum

    return bestSolIndecies, bestSolEval, valuesMulti

'''
Multi Simulated annealing algotrithm, runs the SA for multiple times, takes the best solution at the end
Takes: tries - the number of iterations of the SA algorithm, setOfColours, temperature and simulatedAnnealingIterations
Returns: best solution, the sum of distances of the best, list with evaluations of all the solutions
'''

def multi_sa(tries, setOfColours,temperatureAnnealing,simulatedAnnealingIterations):

    bestSolIndecies = []
    bestSolEval = 0
    valuesMulti = []

    for i in range(0, tries):
        annealingIndecies, sumOfEval, valuesAnnealing = annealing(test_colours,temperatureAnnealing, simulatedAnnealingIterations)
        valuesMulti.append(sumOfEval)
        if len(bestSolIndecies) < 1:
            bestSolIndecies = annealingIndecies[:]
            bestSolEval = sumOfEval
            continue;
        if (bestSolEval > sumOfEval):
            bestSolIndecies = annealingIndecies[:]
            bestSolEval = sumOfEval

    return bestSolIndecies, bestSolEval, valuesMulti

'''
Simple implementation of simulated annealing algorithm
Takes setOfColours, temperature and number of Iterations, for every time the temperature drops
Returns best solution, the sum distance of the best solution and the list of all solutions
'''
def annealing(setOfColours,temperature,numberOfIterations):
    #shuffling the list of colours so that i can have a random start
    sol = random.sample(range(len(setOfColours)), len(setOfColours)) # list of colours for testing
    solListEval = []

    tempmin = 0.00001
    alpha = 0.5

    eval_sol = evaluate(sol)
    solListEval.append(eval_sol)

    while (tempmin < temperature):
        for i in range(numberOfIterations):
            #taking the first/best solution and inverting 2 indexes to create new neighbour then evaluating it
            neighbour_sol = random_neighbour(sol)
            eval_neighbour_sol = evaluate(neighbour_sol)
            #if the new solution is better, then the new solution becomes the best one
            deltaEval = eval_neighbour_sol - eval_sol

            if (eval_neighbour_sol < eval_sol):
                #print "eval_neighbour_sol e ", eval_neighbour_sol
                sol = deepcopy(neighbour_sol)
                eval_sol = eval_neighbour_sol
            else:
                if((math.exp(-deltaEval/temperature)) > random.uniform(0,1)):
                    sol = deepcopy(neighbour_sol)
                    eval_sol = eval_neighbour_sol

            solListEval.append(eval_sol)
        temperature = temperature * alpha

    return sol, eval_sol, solListEval


#####_______main_____########################################################

# Get the directory where the file is located
dir_path = os.path.dirname(os.path.realpath(__file__))
os.chdir(dir_path) # Change the working directory so we can read the file

ncolors, colours = read_file('colours.txt')  # Total number of colours and list of colours

test_size = 100 # Size of the subset of colours for testing
test_colours = colours[0:test_size]  # list of colours for testing

#permutation = random.sample(range(test_size), test_size) # produces random pemutation of lenght test_size, from the numbers 0 to test_size -1
#plot_colours(test_colours, permutation)


#####Greedy algorithm #####

print "------ Greedy algorithm for ",test_size," sample-----"
firstRandomColor = random.randint(0,test_size)
start = time.time()
greedyColoursList, sumOfGreedySol = constructive(test_colours,firstRandomColor)
end = time.time()
print "The time of execution for the greedy algotrithm with ",test_size," colours is: ", end-start, ". The sum of the distances is: ", sumOfGreedySol
plot_colours(test_colours, greedyColoursList)


#this is to get the range for the sum of distance of Greedy algorithm solution
'''
resultCollectorGreedySol = []

while(len(resultCollectorGreedySol)<test_size):
    firstRandomColor = random.randint(0,test_size)
    greedyColoursList, sumOfGreedySol = constructive(test_colours,firstRandomColor)
    resultCollectorGreedySol100.append(sumOfGreedySol)

minSolution = min(resultCollectorGreedySol)
maxSolution = max(resultCollectorGreedySol)

print "The approximate range of solution for distance with the greedy algorithm for ",test_size," colours is between ", minSolution ," and ",maxSolution
'''

hillclimbingIterations = 3000
multiIterations = 30

#### Hill Climbing algotrithm ####
print "------ Hill Climbing algorithm for ",test_size," sample-----"
start = time.time()
hillclimbingIndexes, sumofHillclimbingSol, valuesHillClimb = hill_climbing(test_colours, hillclimbingIterations)
end = time.time()
print "The time of execution for the hill-climbing with",hillclimbingIterations,"iterations for ",test_size," sample list is: ", end - start
print "The sum of distances for hill_climbing is: ", sumofHillclimbingSol
plot_colours(test_colours, hillclimbingIndexes)

#populating a list from 1 until the the lenght of valuesHillClimb to be used in the graph
improvements = []
for iter in range(len(valuesHillClimb)):
    improvements.append(iter)

plt.figure()
plt.plot(improvements,valuesHillClimb)
plt.ylabel('values')
plt.xlabel('number of times iterated')
plt.title("Improvements in Hill Climbing algotrithm")
plt.show()

'''
##### Multi Hill Climbing #####
print "------ Multi Hill Climbing algorithm for ", multiIterations ," times with Hill Climbing at ",hillclimbingIterations," iterations -----"
start = time.time()
bestSolIndexes, bestSolEval, valuesMulti = multi_hc(multiIterations, test_colours, hillclimbingIterations)
end = time.time()
print "The time of the multi_hc is: ", end - start
print "The best solution's sum of distances: ", bestSolEval
print "The Standard Deviation is: ", np.std(valuesMulti)
print "The mean is: ", np.mean(valuesMulti)
print "The median is: ", np.median(valuesMulti)
fig4, ax4 = plt.subplots()
ax4.set_title('Multi Hill Climbing boxplot')
ax4.boxplot(valuesMulti, showfliers=False)
plt.show()
plot_colours(test_colours, bestSolIndexes)


temperatureAnnealing = 3
simulatedAnnealingIterations = 400


#### Simulated Annealing ####
print "------ Simulated Annealing algorithm for ",test_size," sample-----"
start = time.time()
annealingIndecies, sumOfAnnealingEval, valuesAnnealing = annealing(test_colours,temperatureAnnealing, simulatedAnnealingIterations)
end = time.time()
print "The time of execution for the annealing with ",simulatedAnnealingIterations," iterations and temperature of ",temperatureAnnealing," for ",test_size," sample list is: ", end - start
print "The sum of distances for simulated_annealing is: ", sumOfAnnealingEval
plot_colours(test_colours, annealingIndecies)

print "------ Multi Simulated Annealing algorithm for ", multiIterations ," times with SA at ",simulatedAnnealingIterations," iterations -----"
start = time.time()
bestSolIndecies, bestSolEval, valuesMulti = multi_sa(multiIterations, test_colours,temperatureAnnealing, simulatedAnnealingIterations)
end = time.time()
print "The time of the multi_hc is: ", end - start
print "The best solution's sum of distances: ", bestSolEval
print "The Standard Deviation is: ", np.std(valuesMulti)
print "The mean is: ", np.mean(valuesMulti)
print "The median is: ", np.median(valuesMulti)
plot_colours(test_colours, bestSolIndecies)
fig4, ax4 = plt.subplots()
ax4.set_title('Multi Simulated annealing boxplot')
ax4.boxplot(valuesMulti, showfliers=False)
plt.show()
'''
