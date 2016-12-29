#ifndef __FUNKY_H__
#define __FUNKY_H__

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <sstream>
#include <utility>
#include <vector>
#include <assert.h>


/* STEP 1 */

class Expression;

class Function {
  private:
    std::string name;
    std::vector<std::pair<std::string, double *> > params;
    const char * rawExpression;
    Expression * expression;
  public:
    Function(std::string fname, std::vector<std::pair<std::string, double *> > vars);
    ~Function();
    void setRawExpression(char * re);
    void setExpression(Expression * e); 
    size_t getParamCount() const;
    bool isParam(std::string & var) const;
    void setParams(std::vector<double> vals);
    double evaluate() const;
    void printDefinition() const;
    double getParamValue(const std::string & var) const;
    void printName() const;
    std::string getName() const;
    const char * getRawExpr() const;
    std::vector<std::string> getParams() const;
};

class Expression {
  public:
    virtual ~Expression() {}
    virtual double evaluate() const = 0;
};

class NumExpression : public Expression {
  private:
    double number;
  public:
    NumExpression(double num) : number(num) {}
    virtual ~NumExpression() {}
    virtual double evaluate() const {
      return number;
    }
};

class VarExpression : public Expression {
  private:
    std::string id;
    Function * parent;
  public:
    VarExpression(std::string var, Function * f) : id(var), parent(f) {}
    virtual ~VarExpression() {}
    virtual double evaluate() const {
      return parent->getParamValue(id);
    }
};

class OpExpression : public Expression {
  protected:
    Expression * lho;
    Expression * rho;
  public:
    OpExpression(Expression * operand) : lho(operand), rho(NULL) {}
    OpExpression(Expression * lhs, Expression * rhs) : lho(lhs), rho(rhs) {}
    virtual ~OpExpression() {
      delete lho;
      delete rho;
    }
};

class PlusExpression : public OpExpression {
  public:
    virtual ~PlusExpression() {}
    PlusExpression(Expression * lhs, Expression * rhs) : OpExpression(lhs, rhs) {}
    virtual double evaluate() const {
      return lho->evaluate() + rho->evaluate();
    }
};

class MinusExpression : public OpExpression {
  public:
    virtual ~MinusExpression() {}
    MinusExpression(Expression * lhs, Expression * rhs) : OpExpression(lhs, rhs) {}
    virtual double evaluate() const {
      return lho->evaluate() - rho->evaluate();
    }
};

class TimesExpression : public OpExpression {
  public:
    virtual ~TimesExpression() {}
    TimesExpression(Expression * lhs, Expression * rhs) : OpExpression(lhs, rhs) {}
    virtual double evaluate() const {
      return lho->evaluate() * rho->evaluate();
    }
};

class DivideExpression : public OpExpression {
  public:
    virtual ~DivideExpression() {}
    DivideExpression(Expression * lhs, Expression * rhs) : OpExpression(lhs, rhs) {}
    virtual double evaluate() const {
      return lho->evaluate() / rho->evaluate();
    }
};

class PowExpression : public OpExpression {
  public:
    virtual ~PowExpression() {}
    PowExpression(Expression * lhs, Expression * rhs) : OpExpression(lhs, rhs) {}
    virtual double evaluate() const {
      return pow(lho->evaluate(), rho->evaluate());
    }
};

class SqrtExpression : public OpExpression {
  public:
    virtual ~SqrtExpression() {}
    SqrtExpression(Expression * square) : OpExpression(square) {}
    virtual double evaluate() const {
      return sqrt(lho->evaluate());
    }
};


/* STEP 2 */

struct numint {
  Function * function;
  double stepsize;
  std::vector<double> los;
  std::vector<double> his;
  size_t count;
};
typedef struct numint numint_t;

/* STEP 3 */

struct mcint {
  Function * function;
  unsigned long n_iter;
  std::vector<double> los;
  std::vector<double> his;
  size_t count;
};
typedef struct mcint mcint_t;

/* STEP 4 */

struct gradient {
  Function * function;
  double gamma;
  double convergeDist;
  std::vector<double> startPoint;
  size_t count;
};
typedef struct gradient gradient_t;

#endif
