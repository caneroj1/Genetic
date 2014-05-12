genetic
=======

I'm experimenting with genetic algorithms because they're pretty cool.


I found this interesting resource at http://www.ai-junkie.com/ga/intro/gat2.html that talks a bit about genetic algorithms, serving as a nice overview, but they also posed an interesting problem to exercise understanding of genetic algorithms.

The problem reads: Given the digits 0 through 9 and the operators +, -, * and /,  find a sequence that will represent a given target number. The operators will be applied sequentially from left to right as you read.

So pretty much a chromosome in this "population" is a series of bits, where groupings of 4 bits represent a gene, encoding either a digit in 0-9 or one of the arithmetic operators. 

A genetic algorithm will apply the principles from the Theory of Evolution in order to take a population of chromosomes and evolve them until one of them approaches the solution. One of the main challenges encountered so far was parsing the chromosome in order to derive an expression that could be evaluated arithmetically. The value of that expression is used to compute the fitness of each chromosome.

There is still a lot to be implemented in this project, especially crossovers, mutations, and the actual evolution parts.