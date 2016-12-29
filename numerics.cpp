#include <iostream>
#include <ctype.h>
#include <cstdlib>
#include <cstdio>
#include <map>
#include <string.h>
#include <stdio.h>
#include <sstream>
#include <string>
#include <utility>
#include <vector>
#include "funky.h"

/* HELPER FUNCTIONS */
void skipSpace(const char ** strp);

int predefined(std::string fname, const std::vector<Function *> fs);

bool valuesAreOK(std::vector<double> values, std::string fname, const std::vector<Function *> * fs, unsigned * linenumber, char type='e');

Function * parseDefine(const char ** strp, unsigned * linenumber, const std::vector<Function *> * fdefs); 

std::vector<double> parseEvaluate(const char ** strp, std::string * funcname, const std::vector<Function *> * fdefs);

bool commandIs(const char * command, const char ** line);

void printThisDefinition(Function * f, std::vector<Function *> fs);

void printEvaluation(const char * str, std::vector<double> params, std::string fname, const std::vector<Function *> * fdefs);

void checkToPrintWarning(const char ** strp, unsigned * linenumber);

std::string parseVar(const char ** strp);

numint_t * parseNumInt(const char ** strp, Function * f, size_t count);

bool numIntStructOK(numint_t * ni, unsigned * linenumber);

double numericalIntegration(numint_t * ni, size_t n_iter);

void manipulateNumIntStruct(numint_t * ni);

void printIntegration(const char * str, double * answer, char type='n');

size_t determineLoopCount(numint_t * ni);

mcint_t * parseMCInt(const char ** strp, Function * f, size_t count);

bool mcIntStructOK(mcint_t * mci, unsigned * linenumber);

double monteCarloIntegration(mcint_t * mci);

gradient_t * parseGradient(const char ** strp, Function * f, size_t count);



/* MAIN */
int main(int argc, char ** argv) {
  if (argc != 2) {
    std::cerr << "Usage: ./numerics inputfile" << std::endl;
    return EXIT_FAILURE;
  }

  FILE * commandfile;
  if ((commandfile = fopen(argv[1], "r")) == NULL) {
    std::cerr << "Error. Could not open " << argv[1];
  }

  std::vector<Function *> functions;
  Function * func;
  numint_t * numint; // declare struct for numerical integral
  mcint_t * mcint; // struct for mc integral
  gradient_t * gradient; // struct for gradient asc/desc

  char * line = NULL;
  size_t sz;
  unsigned linenum = 1;
  while (getline(&line, &sz, commandfile) != -1) {
    const char * temp = line;
    if (commandIs("define", &temp)) {
      if ((func = parseDefine(&temp, &linenum, &functions)) != NULL) {
        if (predefined(func->getName(), functions) != -1) {
          std::cerr << "Error on line " << linenum << ". Already defined function ";
          func->printName();
          std::cerr << " as: ";
          printThisDefinition(func, functions);
          delete func;
          linenum++;
          continue;
        }
        functions.push_back(func);
        functions.back()->printDefinition();
        checkToPrintWarning(&temp, &linenum);
      }
    }
    else if (commandIs("evaluate", &temp)) {
      skipSpace(&temp);
      const char * evalstr = temp;
      std::string funcname;
      std::vector<double> input = parseEvaluate(&temp, &funcname, &functions);
      if (valuesAreOK(input, funcname, &functions, &linenum)) {
        printEvaluation(evalstr, input, funcname, &functions);
        checkToPrintWarning(&temp, &linenum);
      }
    }
    else if (commandIs("numint", &temp)) {
      skipSpace(&temp);
      const char * intstr = temp;
      std::string funcname = parseVar(&temp);
      int ind;
      if ((ind = predefined(funcname, functions)) == -1) {
        std::cerr << "Error on line " << linenum << ". Parsed function '" << funcname << "' has not been defined yet." << std::endl;
        linenum++;
        continue;
      }
      numint = parseNumInt(&temp, functions[ind], functions[ind]->getParamCount());
      if (numint == NULL) {
        std::cerr << "Error on line " << linenum << ". Cannot integrate function with no parameters." << std::endl;
        linenum++;
        continue;
      }
      if (valuesAreOK(numint->los, funcname, &functions, &linenum, 'n')) {
        if (valuesAreOK(numint->his, funcname, &functions, &linenum, 'n')) {
          if (!numIntStructOK(numint, &linenum)) {
            delete numint;
            linenum++;
            continue;
          }
          manipulateNumIntStruct(numint);
          size_t loops = determineLoopCount(numint);
          double answer = numericalIntegration(numint, loops);
          printIntegration(intstr, &answer);
          checkToPrintWarning(&temp, &linenum);
        }
      }
      delete numint;
    }
    else if (commandIs("mcint", &temp)) {
      skipSpace(&temp);
      const char * intstr = temp;
      std::string funcname = parseVar(&temp);
      int ind;
      if ((ind = predefined(funcname, functions)) == -1) {
        std::cerr << "Error on line " << linenum << ". Parsed function '" << funcname << "' has not been defined yet." << std::endl;
        linenum++;
        continue;
      }
      mcint = parseMCInt(&temp, functions[ind], functions[ind]->getParamCount());
      if (mcint == NULL) {
        std::cerr << "Error on line " << linenum << ". Cannot integrate function with no parameters." << std::endl;
        linenum++;
        continue;
      }
      if (valuesAreOK(mcint->los, funcname, &functions, &linenum, 'm')) {
        if (valuesAreOK(mcint->his, funcname, &functions, &linenum, 'm')) {
          if (!mcIntStructOK(mcint, &linenum)) {
            delete mcint;
            linenum++;
            continue;
          }
          double answer = monteCarloIntegration(mcint);
          printIntegration(intstr, &answer, 'm');
          checkToPrintWarning(&temp, &linenum);
        }
      }      
      delete mcint;
    }
    else if (commandIs("max", &temp)) {
      skipSpace(&temp);
      // const char * gradstr = temp;
      std::string funcname = parseVar(&temp);
      int ind;
      if ((ind = predefined(funcname, functions)) == -1) {
        std::cerr << "Error on line " << linenum << ". Parsed function '" << funcname << "' has not been defined yet." << std::endl;
        linenum++;
        continue;
      }
      gradient = parseGradient(&temp, functions[ind], functions[ind]->getParamCount());
      if (gradient == NULL) {
        std::cerr << "Error on line " << linenum << ". Cannot integrate function with no parameters." << std::endl;
        linenum++;
        continue;
      }
      std::cout << "Actual code to find maximum not reached." << std::endl;
      delete gradient;
    }
    else if (commandIs("min", &temp)) {
      skipSpace(&temp);
      // const char * gradstr = temp;
      std::string funcname = parseVar(&temp);
      int ind;
      if ((ind = predefined(funcname, functions)) == -1) {
        std::cerr << "Error on line " << linenum << ". Parsed function '" << funcname << "' has not been defined yet." << std::endl;
        linenum++;
        continue;
      }
      gradient = parseGradient(&temp, functions[ind], functions[ind]->getParamCount());
      if (gradient == NULL) {
        std::cerr << "Error on line " << linenum << ". Cannot integrate function with no parameters." << std::endl;
        linenum++;
        continue;
      }
      std::cout << "Actual code to find minimum not reached." << std::endl;
      delete gradient;
    }
    else {
      std::cerr << "Error. Invalid command on line " << linenum << "." << std::endl;
    }

    linenum++;
  }

  fclose(commandfile);

  free(line);

  for (size_t i = 0; i < functions.size(); i++) {
    delete functions[i];
  }
  return EXIT_SUCCESS;
}
