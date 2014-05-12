//
//
//  Gene.cpp
//  Created by Joseph Canero on 5/12/14.
//
//

#include <iostream>
#include <bitset>
#include <string>
#include <cmath>
#include <queue>

//  each chromosome will be made up of 9 genes, each made up of 4 bits
//  we will have 10 chromosomes
#define NGENES 9                    //  genes per chromosome
#define NBITS 4                     //  bits per gene
#define NCHROME 10                  //  number of chromosomes
#define TARGET 23                   //  target number for evolution
const int CBITS = NGENES * NBITS;   //  bits per chromosome

using namespace std;

//  function to calculate the fitness
float getFitness(bitset<CBITS>);

int main() {
    //  need an array of bitsets to hold the population
    bitset<CBITS> genepool[NCHROME];
    
    //  initialize the gene pool
    //  seed RNG
    srand(time(NULL));
    for (int i = 0; i < NCHROME; i++) {
        //  generate a random number in the range of 2^n, where n is the number of bits per chromosome
        genepool[i] = rand() % (int)(pow(2, CBITS));
    }
    
    //  compute the fitness for each chromosome in the population
    for (int i = 0; i < NCHROME; i++) {
        cout << getFitness(genepool[i]) << endl;
    }
    
    return 0;
}

//  function to compute the fitness of a given chromosome
//  almost every gene has a mapping but two
//  genes map to a digit in the range of 0-9 and the arithmetic operators: +, -, *, /
//  a chromosome must also be of the form operand, operator, operand, operator, etc.
float getFitness(bitset<CBITS> chromosome) {
    
    //  operands for the expression,
    queue<int> operands;
    
    //  enum for the valid operators
    enum OperatorType {ADD, SUB, MUL, DIV};

    //  this will store the operators of the expression
    queue<OperatorType> expression;
    
    //  convert each chromosome to its string representation in order to initialize a new gene
    //  from groupings of 4 in the string
    string conversion = chromosome.to_string<char, std::string::traits_type,std::string::allocator_type>();
    
    //  this boolean value indicates when we should add an operand to the queue
    //  if it is true, we are looking for an operand, if false we push the next
    //  operator we see
    bool operand = true;
    
    //  we need to parse the bits of the chromosome 4 at a time
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
    
    //  use the two queues to create the expression and evaluate the fitness
    int op1, op2;           //  operands for the expression
    float fitness = 0;      //  the fitness of the equation as determined by how close it gets to TARGET
    int round = 1;          //  the algorithm proceeds differently on the first iteration, this keeps track
    bool cont = true;       //  bool to indicate if the expression has been formed correctly up to current
                            //  point and means we should proceed
    
    do {
        //  first iteration
        if (round == 1) {
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
        round++;                                    //  increment the round counter
    } while (cont);

    //  return the result of the expression
    //  fitness is how close it comes to calculating the value of the target
    return fitness/TARGET;
}