#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include <utility>
#include <vector>
#include <assert.h>
#include "funky.h"

Function::Function(std::string fname, std::vector<std::pair<std::string, double *> > vars) : name(std::string(fname)), params(vars), rawExpression(NULL), expression(NULL) {}

Function::~Function() {
  for (size_t i = 0; i < getParamCount(); i++) {  // pointers to param values
    if (params[i].second != NULL) delete params[i].second;
  }
  if (rawExpression != NULL) delete[] rawExpression;  // raw expression string
  if (expression != NULL) delete expression;  // expression
}

void Function::setRawExpression(char * re) {
  rawExpression = re;
  return;
}

void Function::setExpression(Expression * e) {
  expression = e;
  return;
}

size_t Function::getParamCount() const {
  return params.size();
}

void Function::setParams(std::vector<double> vals) {
  assert(vals.size() == params.size());
  for (size_t i = 0; i < params.size(); i++) {
    if (params[i].second == NULL) {
      params[i].second = new double;
      *(params[i].second) = vals[i];
    }
    else {
      *(params[i].second) = vals[i];
    }
  }
  return;
}

double Function::evaluate() const {
  return expression->evaluate();
}

void Function::printDefinition() const {
  std::cout << "defined (" << name;
  for (size_t i = 0; i < params.size(); i++) {
    std::cout << " " << params[i].first;
  }
  std::cout << ")" << std::endl;
}

bool Function::isParam(std::string & var) const {
  for (size_t i = 0; i < params.size(); i++) {
    if (var == params[i].first) {
      return true;
    }
  }

  return false;
}

double Function::getParamValue(const std::string & var) const {
  for (size_t i = 0; i < params.size(); i++) {
    if (var == params[i].first) {
      return *(params[i].second);
    }
  }
  std::cout << "Error. Could not get value for parameter " << var << ". Has not been set yet." << std::endl;
  abort();
}

void Function::printName() const {
  std::cout << name;
  return;
}

std::string Function::getName() const {
  return name;
}

const char * Function::getRawExpr() const {
  return rawExpression;
}

std::vector<std::string> Function::getParams() const {
  std::vector<std::string> ps;

  for (size_t i = 0; i < params.size(); i++) {
    ps.push_back(params[i].first);
  }

  return ps;
}

