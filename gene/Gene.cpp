//
//
//  Gene.cpp
//  Created by Joseph Canero on 5/12/14.
//
//

#include <iostream>
#include <bitset>
#include <string>
#include <sstream>
#include <cmath>
#include <queue>
#include <vector>
#include <limits.h>

//  each chromosome will be made up of 9 genes, each made up of 4 bits
//  we will have 10 chromosomes
#define NGENES 9                        //  genes per chromosome
#define NBITS 4                         //  bits per gene
#define NCHROME 30                      //  number of chromosomes
#define TARGET 23                       //  target number for evolution
#define GENERATIONS 10                  //  number of generations to go through
const int CBITS = NGENES * NBITS;       //  bits per chromosome
enum OperatorType {ADD, SUB, MUL, DIV}; //  enum for the valid operators

using namespace std;

//  create a struct to store a bitset (chromosome) and the associated fitness
struct gene_fit {
    float fitness = INT_MIN;
    bitset<CBITS> gene;
};

//  function for sorting the vector of gene_fits
bool compare(gene_fit g1, gene_fit g2) { return (g1.fitness < g2.fitness); }

//  functions to calculate the fitness, decode genes in chromosome, evaluate corresponding expression,
//  perform crossovers, run tournaments, mutate chromosome, toString, and advance to the next generation
float getFitness(bitset<CBITS>);
float parseChromosome(bitset<CBITS>);
float evaluateExpression(queue<int>, queue<OperatorType>);
void crossover(bitset<CBITS>, bitset<CBITS>, bitset<CBITS>[]);
gene_fit tournament(vector<gene_fit>);
vector<gene_fit> generation(vector<gene_fit>);
bitset<CBITS> mutate(bitset<CBITS>);
string toString(gene_fit);

int main() {
    //  create a vector of these gene and fitness pair structs
    vector<gene_fit> genepool(NCHROME);
    
    //  initialize the gene pool
    //  seed RNG
    srand(time(NULL));
    
    cout << "Generation 0: " << endl;
    
    //  initialize the chromosomes in the vector
    //  generate a random number in the range of 2^n, where n is the number of bits per chromosome
    //  compute the initial fitness
    for (vector<gene_fit>::iterator iter = genepool.begin(); iter != genepool.end(); iter++) {
        iter->gene = rand() % (int)(pow(2, CBITS));
        iter->gene[32] = rand() % 2;
        iter->gene[33] = rand() % 2;
        iter->gene[34] = rand() % 2;
        iter->gene[35] = rand() % 2;
        iter->fitness = getFitness(iter->gene);
        cout << toString(*iter);
    }
    
    //  start the evolution
    for (int gen = 0; gen < GENERATIONS; gen++) {
        cout << "Generation " << gen + 1 << ": " << endl;
        //  advance the generation
        generation(genepool);
        
        //  compute fitness
        for (vector<gene_fit>::iterator iter = genepool.begin(); iter != genepool.end(); iter++) {
            iter->fitness = getFitness(iter->gene);
            cout << toString(*iter);
        }
    }
    
    return 0;
}

//  function to compute the fitness of a given chromosome
//  passes the chromosome to the parse function which will pass the proper queues to the evaluator function
//  the result will be returned as the fitness
float getFitness(bitset<CBITS> chromosome) {
    return parseChromosome(chromosome);
}

//  function to decode a chromosome
//  almost every gene has a mapping except for two
//  genes map to a digit in the range of 0-9 or one of the arithmetic operators: +, -, *, /
//  a chromosome must also be of the form operand, operator, operand, operator, etc.
float parseChromosome(bitset<CBITS> chromosome) {
    //  this will store the operands of the expression
    queue<int> operands;
    
    //  this will store the operators of the expression
    queue<OperatorType> expression;
    
    //  convert each chromosome to its string representation in order to initialize a new gene
    //  from groupings of 4 in the string
    string conversion = chromosome.to_string<char, std::string::traits_type,std::string::allocator_type>();
    
    //  this boolean value indicates when we should add an operand to the queue.
    //  if it is true, we are looking for an operand, if false we push the next
    //  operator we see
    bool operand = true;
    
    //  we need to parse the bits of the chromosome 4 at a time and get their corresponding value
    for (int g = 0; g < NGENES; g++) {
        bitset<NBITS> gene(conversion.substr(g * 4, 4));
        int value = (int)gene.to_ulong();
        
        //  we want to add an operand
        //  if value is less than 10, it corresponds to 0-9
        if(operand) {
            if (value < 10) {
                operands.push(value);
                operand = false;
            }
        }
        //  we want to add an operator
        //  10 = +, 11 = -, 12 = *, 13 = div
        else {
            if (value >= 10) {
                switch (value) {
                    case 10:
                        expression.push(ADD);
                        operand = true;
                        break;
                        
                    case 11:
                        expression.push(SUB);
                        operand = true;
                        break;
                        
                    case 12:
                        expression.push(MUL);
                        operand = true;
                        break;
                        
                    case 13:
                        expression.push(DIV);
                        operand = true;
                        break;
                }
            }
        }
    }
    return evaluateExpression(operands, expression);
}

//  function to evaluate a chromosome based upon the expression it encodes
//  a fitness value will be returned that indicates how close the expression comes to TARGET
//  the two queues are used in unison in order to provide the data for the expression
float evaluateExpression(queue<int> operands, queue<OperatorType> expression) {
    int op1, op2;           //  operands for the expression
    float fitness = 0;      //  the fitness of the equation as determined by how close it gets to TARGET
    bool first = true;      //  the algorithm proceeds differently on the first iteration
    bool cont = true;       //  bool to indicate if the expression has been formed correctly up to current
                            //  point. true means we should proceed
    
    do {
        //  first iteration
        if (first) {
            if (operands.size() >= 2) {             //  pop 2 operands from the stack if they are available
                op1 = operands.front();
                operands.pop();
                op2 = operands.front();
                operands.pop();
            }
            else cont = false;                      //  if not, we cannot proceed
            if (expression.size() > 0 && cont) {    //  pop an operator, if available
                switch (expression.front()) {
                    case ADD:
                        fitness = op1 + op2;
                        break;
                        
                    case SUB:
                        fitness = op1 - op2;
                        break;
                        
                    case MUL:
                        fitness = op1 * op2;
                        break;
                        
                    case DIV:
                        if (op2 != 0) {
                            fitness = op1 / op2;
                        }
                        else cont = false;
                        break;
                }
                expression.pop();
            }
            else cont = false;                      //  if we couldn't pop one, the expression is invalid
        }
        //  other iterations after the first
        else {
            if (operands.size() > 0) {              //  pop one operand at a time since we have the
                op1 = fitness;                      //  previous result available to us
                op2 = operands.front();
                operands.pop();
            }
            else cont = false;
            if (expression.size() > 0 && cont) {    //  pop an operator, if available
                switch (expression.front()) {
                    case ADD:
                        fitness = op1 + op2;
                        break;
                        
                    case SUB:
                        fitness = op1 - op2;
                        break;
                        
                    case MUL:
                        fitness = op1 * op2;
                        break;
                        
                    case DIV:
                        if (op2 != 0) {
                            fitness = op1 / op2;
                        }
                        else cont = false;
                        break;
                }
                expression.pop();
            }
            else cont = false;                      //  else, the expression is malformed
        }
        first = false;                              //  we have moved on from the first round
    } while (cont);
    
    //  return the result of the expression
    //  fitness is how close it comes to calculating the value of the target
    
    //  distance of solution from target
    int dist = abs(fitness - TARGET);
    
    return (float)(TARGET - dist)/TARGET;
}

//  this function will perform a crossover operation between two parents in order to create two child
//  chromosomes from them. a random point along the chromosome will be generated and all bits after that
//  points on both chromosomes will be swapped in order to create the child chromosome
void crossover(bitset<CBITS> parent1, bitset<CBITS> parent2, bitset<CBITS> child[]) {
    //  generate the crossover point
    int crossover_point = rand() % CBITS;
    
    //  the bits on chromosome 1 and chromosome 2 up to the crossover point should be identical on child 1
    //  and child 2 respectively
    for (int i = 0; i < crossover_point; i++) {
        child[0][i] = parent1[i];
        child[1][i] = parent2[i];
    }
    
    //  now we swap the bits past the crossover point
    for (int i = crossover_point; i < CBITS; i++) {
        child[0][i] = parent2[i];
        child[0][i] = parent1[i];
    }
}

//  this function will take the entire genepool and run tournament selection amongst them. individuals from
//  the population will be randomly selected to compete in the tournament. the individual with the highest
//  fitness is the winner of that tournament and will be returned. this individual will serve as a parent.
gene_fit tournament(vector<gene_fit> population) {
    //  the number of individuals in each tournament is a third of the total population size
    int participants = NCHROME / 3;
    
    //  number of individuals selected so far
    int selected = 0;
    
    //  the winner of the tournament
    gene_fit winner;

    while (selected < participants) {
        //  compute a random number in the range of the population size to choose an individual
        //  create the iterator and advance it the necessary number of places from the beginning
        int chosen = rand() % population.size();
        vector<gene_fit>::iterator iter = population.begin();
        advance(iter, chosen-1);
        
        if (winner.fitness < iter->fitness) {
            winner.gene = iter->gene;
            winner.fitness = iter->fitness;
        }
        
        //  erase the chosen element from the population
        population.erase(iter);
        selected++;
    }
    
    //  return the winner of the tournament
    return winner;
}

//  this function will advance the population to the next generation
//  it will run the appropriate number of tournaments, use the winners to create the children, and
//  create a new population for the next generation. it will return a vector of the next population
vector<gene_fit> generation(vector<gene_fit> population) {
    //  we will run 16 tournaments
    int num_tournaments = 16;
    int chosen;
    
    //  copy of population for sorting purposes
    vector<gene_fit> copy = population;
    
    //  parent1, parent2
    gene_fit parent1, parent2;
    
    //  array of children to be passed to the crossover function
    bitset<CBITS> children[2];
    
    //  vectors that will hold the parents to be used for crossover and the next gen of individuals
    vector<gene_fit> parents, next_gen;
    
    //  run the tournaments and add the winners to the next generation
    for (int i = 0; i < num_tournaments; i++) {
        parents.push_back(tournament(population));
    }
    
    //  select parents randomly from the parents vector to perform crossover
    for (int i = 0; i < num_tournaments/2; i++) {
        //  choose the first parent
        chosen = rand() % parents.size();
        vector<gene_fit>::iterator iter = parents.begin();
        advance(iter, chosen-1);
        
        //  assign values, remove the parent from the vector
        parent1 = *iter;
        parents.erase(iter);
        
        //  choose the second parent
        chosen = rand() % parents.size();
        iter = parents.begin();
        advance(iter, chosen-1);
        
        //  assign values, remove the parent from the vector
        parent2 = *iter;
        parents.erase(iter);
        
        //  crossover between the two parents
        crossover(parent1.gene, parent2.gene, children);
        
        //  copy the bitset from the crossover to the new child and add it to the next generation
        //  do it for both children
        gene_fit *c = new gene_fit;
        c->gene = children[0];
        next_gen.push_back(*c);
        delete c;
        
        gene_fit *d = new gene_fit;
        d->gene = children[1];
        next_gen.push_back(*d);
        delete d;
    }
    
    //  at this point, we only have 16 children in the next generation, we need to choose the 14 most fit
    //  individuals from the last generation
    sort(copy.begin(), copy.end(), compare);
    for (vector<gene_fit>::iterator iter = copy.begin(); iter < copy.begin()+14; iter++) {
        next_gen.push_back(*iter);
    }
    
    //  now that we have the full set of individuals for the next generation, we should do mutation on them
    for (vector<gene_fit>::iterator iter = next_gen.begin(); iter != next_gen.end(); iter++) {
        iter->gene = mutate(iter->gene);
    }
    return next_gen;
}

//  this function will mutate a chromosome. the kind of mutation implemented here is a point mutation in
//  single bits. the probability that a single bit will be mutated is 1 / (CBITS), so generally 1 bit per
//  chromosome is mutated. a mutation in a bit will make 0 become 1 and a 1 become 0.
bitset<CBITS> mutate(bitset<CBITS> chromosome) {
    int m;
    for (int i = 0; i < CBITS; i++) {
        //  generate a number in the range of 0-CBITS. A 0 indicates a mutation will occur.
        m = rand() % CBITS;
        if (m == 0) {
            chromosome.flip(i);
        }
    }
    return chromosome;
}

//  this function accepts a single gene_fit item and returns a neatly formatted string of that chromosome
string toString(gene_fit gene) {
    //  string stream for the results
    ostringstream result;
    result << "Fitness: " << gene.fitness << endl << "Expression: ";
    
    //  convert the gene's chromosome to a string so it can be parsed
    string conversion = gene.gene.to_string<char, std::string::traits_type,std::string::allocator_type>();
    
    //  this boolean value indicates when we should add an operand to the queue.
    //  if it is true, we are looking for an operand, if false we push the next
    //  operator we see
    bool operand = true;
    
    //  we need to parse the bits of the chromosome 4 at a time and get their corresponding value
    for (int g = 0; g < NGENES; g++) {
        bitset<NBITS> gene(conversion.substr(g * 4, 4));
        int value = (int)gene.to_ulong();
        
        //  we want to add an operand
        //  if value is less than 10, it corresponds to 0-9
        if(operand) {
            if (value < 10) {
                result << value << " ";
                operand = false;
            }
        }
        //  we want to add an operator
        //  10 = +, 11 = -, 12 = *, 13 = div
        else {
            if (value >= 10) {
                switch (value) {
                    case 10:
                        result << " + ";
                        operand = true;
                        break;
                        
                    case 11:
                        result << " - ";
                        operand = true;
                        break;
                        
                    case 12:
                        result << " * ";
                        operand = true;
                        break;
                        
                    case 13:
                        result << " / ";
                        operand = true;
                        break;
                }
            }
        }
    }
    result << "\n\n";
    return result.str();
}
