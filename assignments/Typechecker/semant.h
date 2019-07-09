#ifndef SEMANT_H_
#define SEMANT_H_

#include <assert.h>
#include <iostream>
#include <list>
#include <vector> 
#include <map> 
#include "cool-tree.h"
#include "stringtab.h"
#include "symtab.h"
#include "list.h"

#define TRUE 1
#define FALSE 0

// This is a structure that may be used to contain the semantic
// information such as the inheritance graph.  You may use it or not as
// you like: it is only here to provide a container for the supplied
// methods.

ostream& semant_error(Class_ c);
ostream& semant_error(tree_node* node);
ostream& semant_error(Symbol filename, tree_node* node);
ostream& semant_error();

bool isBasic(char* name);
void buildTable(Classes classes);
void buildMethods();
void checkMethods();
int errors();
void typecheck(Classes classes);
void checkMain(Classes classes);
Class_ lca(Class_ c1, Class_ c2);
std::vector<Class_> getInheritanceChain(Symbol name);
std::vector<Class_> getInheritanceChain(Class_ c);
bool subtype(Symbol type1, Symbol type2);

class Graph {
private:
  int V; // number of vertices
  std::vector<std::list<int> > graph;
  bool hasDirectedCycle;
  std::vector<int> cycle;

  void dfs(int v, 
	   std::vector<bool>& marked, 
	   std::vector<bool>& onStack,
	   std::vector<int>& edgeTo);

public:
  Graph(int V);
  std::list<int> adj(int v);
  int getV();
  void addEdge(int i, int j); // graph vertices as numbers starting from 0
  void checkCycle();
  bool hasCycle();
  std::vector<int> getCycle();
};

struct char_comp {
  bool operator()(char const* a, char const* b) const;
};

class GraphChecker {
private:
  // Inheritance graph
  Graph graph;
  int count;
  Classes classes;

  // error routines
  void install_basic_classes();

  // check routines
  void checkCycle(Classes classes);
  void checkInheritance(Classes classes);
  void checkDefinitions(Classes classes);
public:
  GraphChecker(Classes);
};


#endif

