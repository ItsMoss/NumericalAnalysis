#include <iostream>
#include <iomanip>
#include <ctype.h>
#include <cstdlib>
#include <cstdio>
#include <map>
#include <cmath>
#include <string.h>
#include <stdio.h>
#include <sstream>
#include <string>
#include <utility>
#include <vector>
#include "funky.h"


Expression * parseExpr(const char ** strp, Function * parent, const std::vector<Function *> * fdefs);

Expression * parseOp(const char ** strp, Function * parent, const std::vector<Function *> * fdefs);

std::string parseVar(const char ** strp);

std::vector<double> parseEvaluate(const char ** strp, std::string * funcname, const std::vector<Function *> * fdefs);

Expression * parseFuncExpr(const char ** strp, std::map<std::string, char *> params, Function * parent, Function * ogparent, const std::vector<Function *> * fdefs);

Expression * parseFuncOp(std::string op, const char ** rawexpr, const char ** strp, std::vector<std::string> params, Function * opparent, Function * parent, const std::vector<Function *> * fdefs);


void skipSpace(const char ** strp) {
  while(isspace(**strp)) {
    *strp = *strp + 1;
  }
}


double parseNumber(const char ** strp) {
  skipSpace(strp);
  char * endp;
  double num = strtod(*strp, &endp);
  if (endp == *strp) {
    std::cerr << "Error. Expected a number, but found '" << **strp << "'." << std::endl;
    return 0.49394;
  }
  *strp = endp;

  return num;
}


unsigned long parseUnsignedNumber(const char ** strp) {
  skipSpace(strp);

  if (**strp == '-') {
    std::cerr << "Error. Cannot have a negative number." << std::endl;
    return 0;
  }

  char * endp;
  unsigned long num = strtoul(*strp, &endp, 10);

  if (endp == *strp) {
    std::cerr << "Error. Expected a number, but found '" << **strp << "'." << std::endl;
    return 0;
  }
  else if (*endp == '.') {
    std::cerr << "Error. Expected an integer but found a float." << std::endl;
    return 0;
  }
  *strp = endp;
 
  return num; 
}


Expression * makeExpr(std::string op, Expression * lhs, Expression * rhs) {
  if (op == "+") {
    return new PlusExpression(lhs,rhs);
  }
  else if (op == "-") {
    return new MinusExpression(lhs,rhs);
  }
  else if (op == "*") {
    return new TimesExpression(lhs,rhs);
  }
  else if (op == "/") {
    return new DivideExpression(lhs,rhs);
  }
  else if (op == "pow") {
    return new PowExpression(lhs, rhs);
  }
  else if (op == "sqrt") {
    return new SqrtExpression(lhs);
  }
  else {
    std::cerr << "Impossible op char: " << op << std::endl;
    abort();
  }
}


bool isValidOp(std::string op) {
  std::string allOps("+-*/");
  return allOps.find(op) != allOps.npos || op == "pow" || op == "sqrt";
}


std::string parseOpSymbol(const char ** strp) {
  std::stringstream ss;
  while (!isspace(**strp)) {
    ss << **strp;
    *strp = *strp + 1;
  }

  return ss.str();
}


int predefined(std::string fname, const std::vector<Function *> fs) {
  for (size_t i = 0; i < fs.size(); i++) {
    if (fname == fs[i]->getName()) {
      return (int)i;
    }
  }
  return -1;
}


Expression * parseFuncExprOp(const char ** strp, std::map<std::string, char *> params, Function * parent, Function * ogparent, const std::vector<Function *> * fdefs) {
  skipSpace(strp);
  std::string op = parseOpSymbol(strp);
  if (!isValidOp(op)) {
    int ind;
    if ((ind = predefined(op, *fdefs)) != -1) {
      const char * predef = (*fdefs)[ind]->getRawExpr();
      std::vector<std::string> fparams = (*fdefs)[ind]->getParams();
      *strp = *strp + 1;
      return parseFuncOp(op, &predef, strp, fparams, (*fdefs)[ind], parent, fdefs);
    }
    else {
      std::cerr << "Invalid op: "<< op<< std::endl;
      return NULL;
    }
  }
  *strp = *strp + 1;
  Expression * lhs = parseFuncExpr(strp, params, parent, ogparent, fdefs);
  if (lhs == NULL) {
    return NULL;
  }
  Expression * rhs = NULL;
  if (op != "sqrt") {
    rhs = parseFuncExpr(strp, params, parent, ogparent, fdefs);
    if (rhs == NULL) {
      delete lhs;
      return NULL;
    }
  }
  skipSpace(strp);
  if (**strp == ')') {
    *strp = *strp + 1;
    return makeExpr(op,lhs,rhs);
  }
  std::cerr <<"Expected ) but found '" << **strp << "'." << std::endl;
  delete lhs;
  if (rhs != NULL) delete rhs;
  return NULL;
}


Expression * parseFuncExpr(const char ** strp, std::map<std::string, char *> params, Function * parent, Function * ogparent, const std::vector<Function *> * fdefs) {
  skipSpace(strp);
  if (**strp == '(') {
   // operation
    *strp = *strp + 1;
    return parseFuncExprOp(strp, params, parent, ogparent, fdefs);
  }
  else if (isdigit(**strp) || ((**strp == '.' || **strp == '-') && isdigit((*strp)[1])) || (**strp == '-' && (*strp)[1] == '.' && isdigit((*strp)[2]))) {
    // number
    char * endp;
    double num = strtod(*strp, &endp);
    if (endp == *strp) {
      std::cerr << "Error. Expected a number, but found '" << **strp << "'." << std::endl;
      return NULL;
    }
    *strp = endp;
    return new NumExpression(num);
  }
  else if (**strp != '\0') {
    // variable
    std::string var = parseVar(strp);
    if (parent->isParam(var)) {
      return parseExpr((const char **)&(params[var]), ogparent, fdefs);
    }
    else {
      std::cerr << "Error. Variable '" << var << "' not included in function: ";
      parent->printDefinition();
      return NULL;
    }
  }
  else {
    // something wrong occurred
    std::cerr << "Error. Null terminator occurred unexpectedly while parsing." << std::endl;
    abort();
  }
}


void parseFOEhelper(const char ** strp, size_t * sz) {
  // either its a nested operation
  if (**strp == '(') {
    while (**strp != ')') {
      *strp = *strp + 1;
      *sz = *sz + 1;
      if (**strp == '(') {
        parseFOEhelper(strp, sz);
      }
    }
    *strp = *strp + 1;
    *sz = *sz + 1;
  }
  // or its a number
  else {
    char * endp;
    strtod(*strp, &endp);  // do NOT actually need return value here!
    if (endp == *strp) {
      std::cerr << "Error. Expected a number, but found '" << **strp << "'." << std::endl;
    }
    else { // SUCCESS
      *sz = *sz + (endp - *strp);
      *strp = endp;
    }
  }

  return;
} 


char * parseFuncOpExpr(const char ** strp, size_t retsz=1) {
  skipSpace(strp);
  parseFOEhelper(strp, &retsz);
  if (retsz == 1) {
    return NULL;
  }
  
  char * retstr = (char *)malloc(retsz * sizeof(*retstr));
  const char * retbegin = (*strp) - retsz + 1;

  strncpy(retstr, retbegin, retsz - 1);
  retstr[retsz-1] = '\0';  // don't forget to NULL TERMINATE

  return retstr;
}


Expression * parseFuncOp(std::string op, const char ** rawexpr, const char ** strp, std::vector<std::string> params, Function * opparent, Function * parent, const std::vector<Function *> * fdefs) {
  std::map<std::string, char *> parmap;
  size_t parcount = params.size();
  char * parex;

  int i = 0;
  while (**strp != ')') {
    parex = parseFuncOpExpr(strp);
    if (parex == NULL) {
      for (size_t j = 0; j < parmap.size(); j++) {
        free(parmap[params[j]]);
      }
      return NULL;
    }
    parmap.insert(std::pair<std::string, char *>(params[i], parex));
    if (parmap.size() > parcount) {
      std::cerr << "Error. Too many parameters passed in function operator (" << op << ") for definiton of " << parent->getName() << "." << std::endl;
      for (size_t j = 0; j < parmap.size(); j++) {
        free(parmap[params[j]]);
      }
      return NULL;
    }

    skipSpace(strp);
    i++;
  }

  if (parmap.size() < parcount) {
    std::cerr << "Error. Too few parameters passed in function operator (" << op << ") for definition of " << parent->getName() << "." << std::endl;
    for (size_t j = 0; j < parmap.size(); j++) {
      free(parmap[params[j]]);
    }
    return NULL;
  }

  *strp = *strp + 1;
  Expression * prize = parseFuncExpr(rawexpr, parmap, opparent, parent, fdefs);
  for (size_t j = 0; j < parmap.size(); j++) {
    free(parmap[params[j]]);
  }

  return prize;
}


Expression * parseOp(const char ** strp, Function * parent, const std::vector<Function *> * fdefs) {
  skipSpace(strp);
  std::string op = parseOpSymbol(strp);
  if (!isValidOp(op)) {
    int ind;
    if ((ind = predefined(op, *fdefs)) != -1) {
      const char * predef = (*fdefs)[ind]->getRawExpr();
      std::vector<std::string> fparams = (*fdefs)[ind]->getParams();
      *strp = *strp + 1;
      return parseFuncOp(op, &predef, strp, fparams, (*fdefs)[ind], parent, fdefs);
    }
    else {
      std::cerr << "Invalid op: "<< op<< std::endl;
      return NULL;
    }
  }
  *strp = *strp + 1;
  Expression * lhs = parseExpr(strp, parent, fdefs);
  if (lhs == NULL) {
    return NULL;
  }
  Expression * rhs = NULL;
  if (op != "sqrt") {
    rhs = parseExpr(strp, parent, fdefs);
    if (rhs == NULL) {
      delete lhs;
      return NULL;
    }
  }
  skipSpace(strp);
  if (**strp == ')') {
    *strp = *strp + 1;
    return makeExpr(op,lhs,rhs);
  }
  std::cerr <<"Expected ) but found '" << **strp << "'." << std::endl;
  delete lhs;
  if (rhs != NULL) delete rhs;
  return NULL;
}


std::string parseVar(const char ** strp) {
  skipSpace(strp);
  if (**strp == '\0') {
    std::cerr << "Error. Null terminator unexpectedly found." << std::endl;
    abort();
  }
  std::stringstream ss;
  while (!isdigit(**strp) && !isspace(**strp)) {
    ss << **strp;
    *strp = *strp + 1;
    if (isdigit(**strp)) {
      std::cerr << "Error. ID contains a digit." << std::endl;
      return ss.str();
    }
    if (**strp == '(' || **strp == ')') {
      break;  // parentheses cannot be contained in a variable
    }
  }

  return ss.str();
}


Expression * parseExpr(const char ** strp, Function * parent, const std::vector<Function *> * fdefs) {
  skipSpace(strp);
  if (**strp == '(') {
    // operation
    *strp = *strp + 1;
    return parseOp(strp, parent, fdefs);
  }
  else if (isdigit(**strp) || ((**strp == '.' || **strp == '-') && isdigit((*strp)[1])) || (**strp == '-' && (*strp)[1] == '.' && isdigit((*strp)[2]))) {
    // number
    char * endp;
    double num = strtod(*strp, &endp);
    if (endp == *strp) {
      std::cerr << "Error. Expected a number, but found '" << **strp << "'." << std::endl;
      return NULL;
    }
    *strp = endp;
    return new NumExpression(num);
  }
  else if (**strp != '\0') {
    // variable
    std::string var = parseVar(strp);
    if (parent->isParam(var)) {
      return new VarExpression(var, parent);
    }
    else {
      std::cerr << "Error. Variable '" << var << "' not included in function: ";
      parent->printDefinition();
      return NULL;
    }
  }
  else {
    // something wrong occurred
    std::cerr << "Error. Null terminator occurred unexpectedly while parsing." << std::endl;
    abort();
  }
}


Function * parseIDs(const char ** strp) {
  std::string name = parseVar(strp);
  if (name == "" || isdigit(**strp)) return NULL;
  std::vector<std::pair<std::string, double *> > parameters;
  skipSpace(strp);
  while (**strp != ')') {
    std::pair<std::string, double *> par;
    par.first = parseVar(strp);
    if (par.first == "" || isdigit(**strp)) return NULL;
    par.second = NULL;
    parameters.push_back(par);
    skipSpace(strp);
  }

  *strp = *strp + 1; // put string pointer just after the ')'

  // Create Function
  return new Function(name, parameters);
}


Function * parseDefine(const char ** strp, unsigned * linenumber, const std::vector<Function *> * fdefs) {
  skipSpace(strp);
  Function * f = NULL;
  if (**strp == '(') {
    // (name id id ...)
    *strp = *strp + 1;
    skipSpace(strp);
    f = parseIDs(strp);
    if (f == NULL) {
      std::cerr << "Error on line " << *linenumber << ". Invalid id - contains a digit." << std::endl;
      return NULL;
    }
  }
  else {
    std::cerr << "Error on line " << *linenumber << ". Expected a '(' to parse." << std::endl;
    return NULL;
  }
  skipSpace(strp);
  if (**strp == '=') {
    *strp = *strp + 1;
    skipSpace(strp);
  }
  else {
    std::cerr << "Error on line " << *linenumber << ". Expected a '=' to parse." << std::endl;
    delete f;
    return NULL;
  }

  // Create function
  char * raw = new char[strlen(*strp) + 1]; // need to dynamically create expression
  strcpy(raw, *strp);
  f->setRawExpression(raw);
  Expression * expr = parseExpr(strp, f, fdefs);
  f->setExpression(expr);
  if (expr == NULL) {
    delete f;
    f = NULL;
  }

  return f;
}


double ParseInput(const char ** strp, std::string * funcname, const std::vector<Function *> * fdefs) {
  skipSpace(strp);

  // number parameter
  if (isdigit(**strp) || ((**strp == '.' || **strp == '-') && isdigit((*strp)[1])) || (**strp == '-' && (*strp)[1] == '.' && isdigit((*strp)[2]))) {
    char * endp;
    double num = strtod(*strp, &endp);
    if (endp == *strp) {
      std::cerr << "Error. Expected a number, but found '" << *strp[0] << "'." << std::endl;
      return 0.49394;  // the magic error handling value !!!
    }
    *strp = endp;
    return num;
  }
  // function parameter
  else if (**strp == '(') {
    std::string funcparname;
    std::vector<double> values = parseEvaluate(strp, &funcparname, fdefs);
    for (size_t i = 0; i < (*fdefs).size(); i++) {
      if (funcparname == (*fdefs)[i]->getName()) {
        (*fdefs)[i]->setParams(values);
        return (*fdefs)[i]->evaluate();
      }
    }

    return 0.49394;
  }

  std::cerr << "Error. Expected to parse an input but found '" << **strp << "'." << std::endl;
  return 0.49394;
}


std::vector<double> parseEvaluate(const char ** strp, std::string * funcname, const std::vector<Function *> * fdefs) {
  std::vector<double> values;
  skipSpace(strp);
  if (**strp == '(') {
    *strp = *strp + 1;
  }

  // parse function name
  skipSpace(strp);
  *funcname = parseVar(strp);
  int ind;
  if ((ind = predefined(*funcname, *fdefs)) == -1) {
    std::cerr << "Error. Cannot evaluate undefined function '" << *funcname << "'." << std::endl;
    return values;
  }

  size_t nparams = (*fdefs)[ind]->getParamCount();

  for (size_t i = 0; i < nparams; i++) {
    values.push_back(ParseInput(strp, funcname, fdefs));
    if (values.back() == 0.49394) {
      return values;
    }
  }

  skipSpace(strp);
  if (**strp == ')') {
    *strp = *strp + 1;
    return values;
  }

  std::cerr << "Error. Cannot evaluate '" << *funcname << "', expected a ')' to parse." << std::endl;
  return values;
}



bool commandIs(const char * command, const char ** line) {
  skipSpace(line);
  size_t L = strlen(command);
  const char * start = *line;

  for (size_t i = 0; i < L; i++) {
    if (command[i] != (char)tolower(**line)) {
      *line = start;
      return false;
    }
    *line = *line + 1;
  }

  // Make sure next character is a space
  if (!isspace(**line)) {
    *line = start;
    return false;
  }
  return true;
}


void printThisDefinition(Function * f, std::vector<Function *> fs) {
  for (size_t i = 0; i < fs.size(); i++) {
    if (f->getName() == fs[i]->getName()) {
      fs[i]->printDefinition();
    }
  }
  return;
}


void printEvaluation(const char * str, std::vector<double> params, std::string fname, const std::vector<Function *> * fdefs) {
  double answer;
  size_t i;

  for (i = 0; i < (*fdefs).size(); i++) {
    if (fname == (*fdefs)[i]->getName()) {
      (*fdefs)[i]->setParams(params);
      answer = (*fdefs)[i]->evaluate();
      break;
    }
  }

  skipSpace(&str);
  size_t L = strlen(str);  // L - 1 is what we want to print (i.e. exclude newline character and nullterminator
  std::string evalStr(str, L-1);

  std::cout << evalStr << " = " << std::setprecision(3) << answer << std::endl;
  std::cout.unsetf(std::ios::floatfield);

  return;
}


bool valuesAreOK(std::vector<double> values, std::string fname, const std::vector<Function *> * fdefs, unsigned * linenumber, char type='e') {
  // Let's check values are not eqaul to magical error number of 0.49394
  // AND that the correct number of values are provided given fname definition
  std::string command, input;

  if (type == 'e') {
    command = "evaluate";
    input = "input parameter";
  }
  else {
    command = "integrate";
    input = "pair of bound parameters";
  }

  for (size_t i = 0; i < (*fdefs).size(); i++) {
    if (fname == (*fdefs)[i]->getName()) {
      if (values.size() != (*fdefs)[i]->getParamCount()) {
        std::cerr << "Error on line " << *linenumber << ". Too few input parameters provided to evaluate function: '" << fname << "'." << std::endl;
        return false;
      }
      break;
    }
  }

  if (!values.size()) return false;  // do not forget an empty value set

  for (size_t i = 0; i < values.size(); i++) {
    if (values[i] == 0.49394) {
      std::cerr << "Error on line " << *linenumber << ". Could not " << command << " due to problem with " << i+1;
      switch (i) {
        case 0:
          std::cerr << "st " << input << "." << std::endl;
          return false;
        case 1:
          std::cerr << "nd " << input << "." << std::endl;
          return false;
        case 2:
          std::cerr << "rd " << input << "." << std::endl;
          return false;
        default:
          std::cerr << "th " << input << "." << std::endl;
          return false;
      }
    }
  }

  return true;
}


void checkToPrintWarning(const char ** strp, unsigned * linenumber) {
  if (**strp == ')') {
    *strp = *strp + 1;
  }

  skipSpace(strp);

  if (isalnum(**strp) || ispunct(**strp)) {
    std::cerr << "Warning on line " << *linenumber << ". Ignoring the following text: " << *strp;
  }
}


void cascadeBack(std::vector<double> * curr, size_t back, numint_t * ni) {
  if ((*curr)[back] > ni->his[back]) {
    if (back == 0) return;
    (*curr)[back] = ni->los[back];  // reset
    back--;
    (*curr)[back] += ni->stepsize;  // increase new back
    cascadeBack(curr, back, ni);
  }
  return;
}


size_t determineLoopCount(numint_t * ni) {
  size_t count = 1;
  // std::cout << "Struct count=" << ni->count << std::endl;
  for (size_t i = 0; i < ni->count; i++) {
    count *= (size_t)((ni->his[i] - ni->los[i]) / ni->stepsize) + 1;
    // std::cout << "current hi=" << ni->his[i] << "and low=" << ni->los[i] << std::endl;
  }
  // std::cout << "Return count=" << count << std::endl;
  return count;
}


double numericalIntegration(numint_t * ni, size_t n_iter) {
  std::vector<double> curr = ni->los;
  double sum = 0.0;
  
  for (size_t i = 0; i < n_iter; i++) {
    ni->function->setParams(curr);
    double result = pow(ni->stepsize, ni->count) * ni->function->evaluate();
    curr.back() += ni->stepsize;
    cascadeBack(&curr, ni->count - 1, ni);
    if (!isnan(result)) {  // ignore all of the nans
      sum += result;
    }
  }

  return sum;
}


numint_t * parseNumInt(const char ** strp, Function * f, size_t count) {
  // Before anything, make sure count is > 0
  if (count == 0) return NULL;

  skipSpace(strp);
  numint_t * ni = new numint_t;

  ni->function = f;
  ni->count = count;
  ni->stepsize = parseNumber(strp);

  for (size_t i = 0; i < ni->count; i++) {
    ni->los.push_back(parseNumber(strp));
    if (ni->los[i] == 0.49394) break;
    ni->his.push_back(parseNumber(strp));
    if (ni->his[i] == 0.49394) break;
  }

  return ni;
}


bool numIntStructOK(numint_t * ni, unsigned * linenumber) {
  if (ni->stepsize <= 0) {
    std::cerr << "Error on line " << *linenumber << ". Stepsize must be greater than zero." << std::endl;
    return false;
  }

  for (size_t i = 0; i < ni->count; i++) {
    if (ni->stepsize > ni->his[i] - ni->los[i]) {
      std::cerr << "Error on line " << *linenumber << ". Stepsize is greater than a bound range." << std::endl;
      return false;
    }
    if (ni->his[i] <= ni->los[i]) {
      std::cerr << "Error on line " << *linenumber << ". Cannot have a decreasing bound range." << std::endl;
      return false;
    }
  }

  return true;
}


void manipulateNumIntStruct(numint_t * ni) {
  double halfstep = ni->stepsize / 2;

  for (size_t i = 0; i < ni->count; i++) {
    ni->los[i] += halfstep;
    ni->his[i] -= halfstep;
  }

  return;
}


void printIntegration(const char * str, double * answer, char type='n') {
  std::string inttype;
  if (type == 'm') {
    inttype = "mcint";
  }
  else {
    inttype = "numint";
  }

  std::string strclean(str, strlen(str) - 1);

  std::cout << inttype << " " << strclean << " = " << std::setprecision(3) << *answer << std::endl;
  std::cout.unsetf(std::ios::floatfield);

  return;
}


mcint_t * parseMCInt(const char ** strp, Function * f, size_t count) {
  // Before anything, make sure count is > 0
  if (count == 0) return NULL;

  skipSpace(strp);
  mcint_t * mci = new mcint_t;

  mci->function = f;
  mci->count = count;
  mci->n_iter = parseUnsignedNumber(strp);

  for (size_t i = 0; i < mci->count; i++) {
    mci->los.push_back(parseNumber(strp));
    if (mci->los[i] == 0.49394) break;
    mci->his.push_back(parseNumber(strp));
    if (mci->his[i] == 0.49394) break;
  }

  return mci;
}


bool mcIntStructOK(mcint_t * mci, unsigned * linenumber) {
  if (mci->n_iter <= 0) {
    std::cerr << "Error on line " << *linenumber << ". Number of random iterations must be greater than zero." << std::endl;
    return false;
  }

  for (size_t i = 0; i < mci->count; i++) {
    if (mci->his[i] <= mci->los[i]) {
      std::cerr << "Error on line " << *linenumber << ". Cannot have a decreasing bound range." << std::endl;
      return false;
    }
  }

  return true;
}


double computeNDArea(mcint_t * mci) {
  double area = 1.0;

  for (size_t i = 0; i < mci->count; i++) {
    area *= (mci->his[i] - mci->los[i]);
  }

  return area;
}


std::vector<double> mcRandom(mcint_t * mci) {
  std::vector<double> randpoint;

  for (size_t i = 0; i < mci->count; i++) {
    double rho = (mci->his[i] - mci->los[i] + 1);
    int randparam = rand() % (int)rho + (int)mci->los[i];
    randpoint.push_back((double)randparam);
  }

  return randpoint;
}


double monteCarloIntegration(mcint_t * mci) {
  double sum = 0.0;
  double mcarea = computeNDArea(mci);

  for (unsigned long i = 0; i < mci->n_iter; i++) {
    std::vector<double> curr = mcRandom(mci);
    mci->function->setParams(curr);
    double result = mci->function->evaluate();
    if (!isnan(result)) {
      sum += result;
    }
  }
  
  return sum / mci->n_iter * mcarea;
}


gradient_t * parseGradient(const char ** strp, Function * f, size_t count) {
  // Before anything, make sure count is > 0
  if (count == 0) return NULL;

  skipSpace(strp);
  gradient_t * gr = new gradient_t;

  gr->function = f;
  gr->count = count;
  gr->gamma = parseNumber(strp);
  gr->convergeDist = parseNumber(strp);

  for (size_t i = 0; i < gr->count; i++) {
    gr->startPoint.push_back(parseNumber(strp));
    if (gr->startPoint[i] == 0.49394) break;
  }

  return gr;
}

