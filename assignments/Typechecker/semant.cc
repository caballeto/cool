
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <string>
#include <set>
#include <algorithm>
#include "semant.h"
#include "tree.h"
#include "cool-tree.h"
#include "utilities.h"

extern int semant_debug;
extern char *curr_filename;

static int semant_errors = 0;
static ostream& error_stream = std::cerr;

static Class_ curr_class = NULL;

std::map<char*, int, char_comp> class_name_map;
std::map<int, Class_> vertice_map;

typedef std::map<Symbol, Class_> ClassMap;
ClassMap classMap;

SymbolTable<Symbol, Symbol> objects;
std::map<Class_, std::vector<method_class*> > methods;

//////////////////////////////////////////////////////////////////////
//
// Symbols
//
// For convenience, a large number of symbols are predefined here.
// These symbols include the primitive type and method names, as well
// as fixed names used by the runtime system.
//
//////////////////////////////////////////////////////////////////////
static Symbol 
    arg,
    arg2,
    Bool,
    concat,
    cool_abort,
    copy,
    Int,
    in_int,
    in_string,
    IO,
    length,
    Main,
    main_meth,
    No_class,
    No_type,
    Object,
    out_int,
    out_string,
    prim_slot,
    self,
    SELF_TYPE,
    Str,
    str_field,
    substr,
    type_name,
    val;
//
// Initializing the predefined symbols.
//
static void initialize_constants(void) {
    arg         = idtable.add_string("arg");
    arg2        = idtable.add_string("arg2");
    Bool        = idtable.add_string("Bool");
    concat      = idtable.add_string("concat");
    cool_abort  = idtable.add_string("abort");
    copy        = idtable.add_string("copy");
    Int         = idtable.add_string("Int");
    in_int      = idtable.add_string("in_int");
    in_string   = idtable.add_string("in_string");
    IO          = idtable.add_string("IO");
    length      = idtable.add_string("length");
    Main        = idtable.add_string("Main");
    main_meth   = idtable.add_string("main");
    //   _no_class is a symbol that can't be the name of any 
    //   user-defined class.
    No_class    = idtable.add_string("_no_class");
    No_type     = idtable.add_string("_no_type");
    Object      = idtable.add_string("Object");
    out_int     = idtable.add_string("out_int");
    out_string  = idtable.add_string("out_string");
    prim_slot   = idtable.add_string("_prim_slot");
    self        = idtable.add_string("self");
    SELF_TYPE   = idtable.add_string("SELF_TYPE");
    Str         = idtable.add_string("String");
    str_field   = idtable.add_string("_str_field");
    substr      = idtable.add_string("substr");
    type_name   = idtable.add_string("type_name");
    val         = idtable.add_string("_val");
}

int errors() {
  return semant_errors;
}

bool isBasic(Symbol name) {
  return name == Int || name == Bool || name == Str;
}

GraphChecker::GraphChecker(Classes classes) : graph(0), count(0) {
  install_basic_classes(); // initialize basic classes
  this->classes = classes;
  
  checkDefinitions(classes);
  if (errors() > 0) return;
  checkInheritance(classes);
  if (errors() > 0) return;
  checkCycle(classes);
}

void GraphChecker::checkDefinitions(Classes classes) {
  count = 5; // number of base classes
 
  // checking redefined classes and classes that extend basic classes
  for (int i = classes->first(); classes->more(i); i = classes->next(i)) {
    Class_ clazz = classes->nth(i);
    Symbol name = clazz->get_name();
    
    // basic class redefinition
    if (isBasic(name) || name == Object|| name == IO || name == SELF_TYPE) {
      semant_error(clazz) << "Redefinition of basic class " << name << "." << std::endl;
    } else if (class_name_map.count(name->get_string()) == 1) {
      semant_error(clazz) << "Class " << name << " was previously defined." << std::endl;
    }
    
    vertice_map[count] = clazz;
    class_name_map[name->get_string()] = count++;
  }
}

void GraphChecker::checkInheritance(Classes classes) {
  for (int i = classes->first(); classes->more(i); i = classes->next(i)) {
    Class_ clazz = classes->nth(i);
    
    // inheriting from base classes
    if (isBasic(clazz->get_parent()) || clazz->get_parent() == SELF_TYPE) {
      semant_error(clazz)<< "Class " << clazz->get_name() << " cannot inherit class " 
			 << clazz->get_parent() << "." << std::endl;
    } 

    if (class_name_map.count(clazz->get_parent()->get_string()) == 0) {
      semant_error(clazz) << "Class " << clazz->get_name() << " inherits from an undefined class "
			  << clazz->get_parent() << "." << std::endl;
    }
  }
}

void GraphChecker::checkCycle(Classes classes) {
  Graph G(count);

  // default inheritance [Int, Bool, String, IO]
  G.addEdge(1, 0);
  G.addEdge(2, 0);
  G.addEdge(3, 0);
  G.addEdge(4, 0);

  // Build graph
  for (int i = classes->first(); classes->more(i); i = classes->next(i)) {
    Class_ clazz = classes->nth(i);
    G.addEdge(class_name_map[clazz->get_name()->get_string()], class_name_map[clazz->get_parent()->get_string()]);
  }

  G.checkCycle(); // check for cycle in a graph

  if (G.hasCycle()) {
    std::vector<int> cycle_vertices = G.getCycle();
    for (std::vector<int>::iterator it = cycle_vertices.begin(); it != cycle_vertices.end(); ++it) {
      Class_ clazz = vertice_map[*it];
      Symbol name = clazz->get_name();
      semant_error(clazz) << "Class " << name << ", or an ancestor of " << name << 
	", is involved in an inheritance cycle." << std::endl;  
    }
  }

  this->graph = G; // assign built graph
}

//Class_ ClassTable::lca(Class_ a, Class_ b) {
//  return vertice_map[graph.lca(class_name_map[a->get_name()->get_string()], 
//			       class_name_map[b->get_name()->get_string()])];
//}

void GraphChecker::install_basic_classes() {

    // The tree package uses these globals to annotate the classes built below.
   // curr_lineno  = 0;
    Symbol filename = stringtable.add_string("<basic class>");
    
    // The following demonstrates how to create dummy parse trees to
    // refer to basic Cool classes.  There's no need for method
    // bodies -- these are already built into the runtime system.
    
    // IMPORTANT: The results of the following expressions are
    // stored in local variables.  You will want to do something
    // with those variables at the end of this method to make this
    // code meaningful.

    // 
    // The Object class has no parent class. Its methods are
    //        abort() : Object    aborts the program
    //        type_name() : Str   returns a string representation of class name
    //        copy() : SELF_TYPE  returns a copy of the object
    //
    // There is no need for method bodies in the basic classes---these
    // are already built in to the runtime system.

    Class_ Object_class =
	class_(Object, 
	       No_class,
	       append_Features(
			       append_Features(
					       single_Features(method(cool_abort, nil_Formals(), Object, no_expr())),
					       single_Features(method(type_name, nil_Formals(), Str, no_expr()))),
			       single_Features(method(copy, nil_Formals(), SELF_TYPE, no_expr()))),
	       filename);

    // 
    // The IO class inherits from Object. Its methods are
    //        out_string(Str) : SELF_TYPE       writes a string to the output
    //        out_int(Int) : SELF_TYPE            "    an int    "  "     "
    //        in_string() : Str                 reads a string from the input
    //        in_int() : Int                      "   an int     "  "     "
    //
    Class_ IO_class = 
	class_(IO, 
	       Object,
	       append_Features(
			       append_Features(
					       append_Features(
							       single_Features(method(out_string, single_Formals(formal(arg, Str)),
										      SELF_TYPE, no_expr())),
							       single_Features(method(out_int, single_Formals(formal(arg, Int)),
										      SELF_TYPE, no_expr()))),
					       single_Features(method(in_string, nil_Formals(), Str, no_expr()))),
			       single_Features(method(in_int, nil_Formals(), Int, no_expr()))),
	       filename);  

    //
    // The Int class has no methods and only a single attribute, the
    // "val" for the integer. 
    //
    Class_ Int_class =
	class_(Int, 
	       Object,
	       single_Features(attr(val, prim_slot, no_expr())),
	       filename);

    //
    // Bool also has only the "val" slot.
    //
    Class_ Bool_class =
	class_(Bool, Object, single_Features(attr(val, prim_slot, no_expr())),filename);

    //
    // The class Str has a number of slots and operations:
    //       val                                  the length of the string
    //       str_field                            the string itself
    //       length() : Int                       returns length of the string
    //       concat(arg: Str) : Str               performs string concatenation
    //       substr(arg: Int, arg2: Int): Str     substring selection
    //       
    Class_ Str_class =
	class_(Str, 
	       Object,
	       append_Features(
			       append_Features(
					       append_Features(
							       append_Features(
									       single_Features(attr(val, Int, no_expr())),
									       single_Features(attr(str_field, prim_slot, no_expr()))),
							       single_Features(method(length, nil_Formals(), Int, no_expr()))),
					       single_Features(method(concat, 
								      single_Formals(formal(arg, Str)),
								      Str, 
								      no_expr()))),
			       single_Features(method(substr, 
						      append_Formals(single_Formals(formal(arg, Int)), 
								     single_Formals(formal(arg2, Int))),
						      Str, 
						      no_expr()))),
	       filename);

    // init base classes mappings
    classMap[Object_class->get_name()] = Object_class;
    classMap[Int_class->get_name()] = Int_class;
    classMap[Bool_class->get_name()] = Bool_class;
    classMap[IO_class->get_name()] = IO_class;
    classMap[Str_class->get_name()] = Str_class;

    class_name_map["Object"] = 0;
    class_name_map["Int"] = 1;
    class_name_map["Bool"] = 2;
    class_name_map["IO"] = 3;
    class_name_map["String"] = 4;

    vertice_map[0] = Object_class;
    vertice_map[1] = Int_class;
    vertice_map[2] = Bool_class;
    vertice_map[3] = IO_class;
    vertice_map[4] = Str_class;
}

// Graph implementation

Graph::Graph(int V) 
  : V(V), graph(std::vector<std::list<int> >(V, std::list<int>())), hasDirectedCycle(false) 
{ }

bool Graph::hasCycle() {
  return hasDirectedCycle;
}

std::list<int> Graph::adj(int v) {
  return graph[v];
}

int Graph::getV() {
  return V;
}

void Graph::addEdge(int i, int j) {
  graph[i].push_back(j);
}

void Graph::checkCycle() {
  std::vector<bool> marked(V, false);
  std::vector<bool> onStack(V, false);
  std::vector<int> edgeTo(V, 0);
  for (int v = 0; v < V; v++) {
    if (!marked[v] && !hasDirectedCycle) {
      dfs(v, marked, onStack, edgeTo);
    }
  }
}

void Graph::dfs(int v, std::vector<bool>& marked, std::vector<bool>& onStack, std::vector<int>& edgeTo) {
  onStack[v] = true;
  marked[v] = true;
  std::list<int> adjacent = adj(v);
  for (std::list<int>::iterator it = adjacent.begin(); it != adjacent.end(); ++it) {
    int w = *it;
    if (hasDirectedCycle) return;
    else if (!marked[w]) {
      edgeTo[w] = v;
      dfs(w, marked, onStack, edgeTo);
    } else if (onStack[w]) {
      hasDirectedCycle = true;
      for (int x = v; x != w; x = edgeTo[x]) {
	cycle.push_back(x);
      }
      cycle.push_back(w);
    }
  }
  onStack[v] = false;
}

std::vector<int> Graph::getCycle() {
  return cycle;
}

// char* comparator
bool char_comp::operator()(char const* a, char const* b) const {
  return strcmp(a, b) < 0;
}

////////////////////////////////////////////////////////////////////
//
// semant_error is an overloaded function for reporting errors
// during semantic analysis.  There are three versions:
//
//    ostream& ClassTable::semant_error()                
//
//    ostream& ClassTable::semant_error(Class_ c)
//       print line number and filename for `c'
//
//    ostream& ClassTable::semant_error(Symbol filename, tree_node *t)  
//       print a line number and filename
//
///////////////////////////////////////////////////////////////////

ostream& semant_error(Class_ c) {
  return semant_error(c->get_filename(), c);
}

ostream& semant_error(tree_node* node) {                                                             
    return semant_error(curr_class->get_filename(), node);
}

ostream& semant_error(Symbol filename, tree_node *t) {
    error_stream << filename << ":" << t->get_line_number() << ": ";
    return semant_error();
}

ostream& semant_error() {                                                 
    semant_errors++;                            
    return error_stream;
} 

Class_ lca(Symbol c1, Symbol c2) {
  std::vector<Class_> chain1 = getInheritanceChain(c1);
  std::vector<Class_> chain2 = getInheritanceChain(c2);
  
  std::reverse(chain1.begin(), chain1.end());
  std::reverse(chain2.begin(), chain2.end());

  std::size_t i;
  for (i = 1; i < std::min(chain1.size(), chain2.size()); i++)
    if (chain1[i] != chain2[i])
      return chain1[i - 1];
  return chain1[i - 1];
}

std::vector<Class_> getInheritanceChain(Symbol name) {
  if (name == SELF_TYPE)
    name = curr_class->get_name();
  return getInheritanceChain(classMap[name]);
}

std::vector<Class_> getInheritanceChain(Class_ c) {
  std::vector<Class_> chain;
  while (c->get_name() != Object) {
    chain.push_back(c);
    c = classMap[c->get_parent()];
  }

  chain.push_back(classMap[Object]);
  return chain;
}

bool subtype(Symbol type1, Symbol type2) {
  if (type1 == SELF_TYPE && type2 == SELF_TYPE) return true;
  if (type1 != SELF_TYPE && type2 == SELF_TYPE) return false;
  if (type1 == SELF_TYPE) type1 = curr_class->get_name();
  if (type1 == No_type && type2 == No_type) return true;
  if (type1 == No_type) return false;

  Class_ c1 = classMap[type1], c2 = classMap[type2];
  std::vector<Class_> chain = getInheritanceChain(c1);
  
  for (std::size_t i = 0; i < chain.size(); i++)
    if (chain[i] == c2)
      return true;
  return false;
}

void checkMain(Classes classes) {
  bool hasMain = false;
  Class_ main = NULL;
  for (int i = classes->first(); classes->more(i); i = classes->next(i)) {
    if (classes->nth(i)->get_name() == Main) {
      hasMain = true;
      main = classes->nth(i);
      break;
    }
  }

  if (!hasMain) {
    semant_error() << "Class Main is not defined." << std::endl;
  } else {
    hasMain = false;
    std::vector<method_class*> _methods = methods[main];
    for (std::size_t i = 0; i < _methods.size(); i++) {
      if (_methods[i]->get_name() == main_meth) {
	hasMain = true;
	break;
      }
    }

    if (!hasMain) {
      semant_error(curr_class) << "No 'main' method in class Main. " << std::endl;
    }
  }
}

void buildTable(Classes classes) {
  for (int i = classes->first(); classes->more(i); i = classes->next(i)) {
    Class_ clazz = classes->nth(i);
    classMap[clazz->get_name()] = clazz;
  }
}

void buildMethods() {
  for (ClassMap::iterator it = classMap.begin(); it != classMap.end(); it++) {
    Features features = it->second->get_features();
    std::vector<method_class*> _methods;
    for (int i = features->first(); features->more(i); i = features->next(i)) {
      if (strcmp(features->nth(i)->constructor(), "method") == 0) {
	method_class* method = static_cast<method_class*>(features->nth(i));
	bool existed = false;
	for (std::size_t j = 0; j < _methods.size(); j++)
	  if (_methods[j]->get_name() == method->get_name())
	    existed = true;
	
	if (existed)
	  semant_error(method) << "Method " << method->get_name() << " is multiply defined." << std::endl;
	else
	  _methods.push_back(static_cast<method_class*>(features->nth(i)));
      }
    }

    methods[it->second] = _methods;
  }
}

void typecheck(Classes classes) {
  buildTable(classes);
  //std::cerr << "Built table." << std::endl;
  buildMethods();
  //std::cerr << "Built methods." << std::endl;
  checkMain(classes);
  //std::cerr << "Checked main." << std::endl;
  checkMethods();
  //std::cerr << "Checked methods." << std::endl;
}

method_class* get_method(Class_ c, Symbol method_name) {
  std::vector<method_class*> _methods = methods[c];
  for (std::size_t i = 0; i < _methods.size(); i++)
    if (_methods[i]->get_name() == method_name)
      return _methods[i];
  return NULL;
}

void checkMethods() {
  for (ClassMap::iterator it = classMap.begin(); it != classMap.end(); it++) {
    if (isBasic(it->first) || it->first == IO || it->first == Object) continue;
    Symbol class_name = it->first;
    curr_class = it->second;
    //std::cerr << "Processing class " << class_name << std::endl;
    std::vector<Class_> chain = getInheritanceChain(curr_class);
    chain.push_back(curr_class);
    for (std::size_t k = 1; k < chain.size(); k++) {
      Class_ ancestor_class = chain[k];
      Features ancestor_features = ancestor_class->get_features();
      objects.enterscope();
      for (int i = ancestor_features->first(); ancestor_features->more(i); i = ancestor_features->next(i)) {
	if (strcmp(ancestor_features->nth(i)->constructor(), "attr") != 0) continue;
	attr_class* attr = static_cast<attr_class*>(ancestor_features->nth(i));
	if (objects.lookup(attr->get_name()) != NULL)
	    semant_error(curr_class) << "Attribute " << attr->get_name() << " is an attribute of an inherited class." << std::endl;
	objects.addid(attr->get_name(), new Symbol(attr->get_type()));
      }
    }

    //std::cerr << "Processed ancestor classes. " << std::endl;
    Features features = curr_class->get_features();
    //std::cerr << "Got features." << std::endl;
    for (int i = features->first(); features->more(i); i = features->next(i)) {
      // std::cerr << "Processing features: " << std::endl;
      if (strcmp(features->nth(i)->constructor(), "method") == 0) {
	method_class* curr_method = static_cast<method_class*>(features->nth(i));
	curr_method->typecheck();
	//std::cerr << "Processed method. " << std::endl;
	//std::cerr << "Processed method: " << curr_method->get_name() << std::endl;
	for (std::size_t k = 1; k < chain.size(); k++) {
	  method_class* ancestor_method = get_method(chain[k], curr_method->get_name());
	  if (!ancestor_method) continue;
	  
	  if (curr_method->get_return_type() != ancestor_method->get_return_type())
	    semant_error(curr_method) << "In redefined method " << curr_method->get_name() << ", return type " << curr_method->get_return_type() 
				      << " is different from original return type " << ancestor_method->get_return_type() << "." << std::endl;

	  Formals curr_formals = curr_method->get_formals();
	  Formals ancestor_formals = ancestor_method->get_formals();
	  
	  int k1 = curr_formals->first(), k2 = ancestor_formals->first();
	  while (curr_formals->more(k1) && ancestor_formals->more(k2)) {
	    if (curr_formals->nth(k1)->get_type() != ancestor_formals->nth(k2)->get_type())
	      semant_error(curr_formals->nth(k1)) << "In redefined method " << curr_method->get_name() 
						  << ", type " << curr_formals->nth(k1)->get_type()
						  << " is different from original type " << ancestor_formals->nth(k2)->get_type() 
						  << "." << std::endl;
	    k1 = curr_formals->next(k1);
	    k2 = ancestor_formals->next(k2);
	    if (curr_formals->more(k1) ^ ancestor_formals->more(k2))
	      semant_error(curr_method) << "Incompatible number of formal parameters in redefined method " << curr_method->get_name() << "." << std::endl;
	  }
	}
      } else {
	attr_class* curr_attr = static_cast<attr_class*>(features->nth(i));
      //std::cerr << "Processing attribute: " << curr_attr->get_name() << std::endl;
	Symbol expr_type = curr_attr->get_init_expr()->typecheck();
	if (classMap.find(curr_attr->get_type()) == classMap.end())
	  semant_error(curr_attr) << "Class " << curr_attr->get_type() << " of attribute " << curr_attr->get_name() << " is undefined." << std::endl;
	else if (classMap.find(expr_type) != classMap.end() && !subtype(expr_type, curr_attr->get_type()))
	  semant_error(curr_attr) << "Inferred type " << expr_type << " of initialization attribute " << curr_attr->get_name() 
				  << " does not conform to declared type " << curr_attr->get_type() << "." << std::endl;
	else if (curr_attr->get_name() == self)
	  semant_error(curr_attr) << "Attribute with illegal name 'self'." << std::endl;
      }
    }

    for (std::size_t k = 1; k < chain.size(); k++)
      objects.exitscope();
  }
}

void method_class::typecheck() {
  objects.enterscope();
  for (int i = formals->first(); formals->more(i); i = formals->next(i)) {
    if (formals->nth(i)->get_name() == self)
      semant_error(formals->nth(i)) << "'self' cannot be the name of formal parameter." << std::endl;
    else if (objects.probe(formals->nth(i)->get_name()))
      semant_error(formals->nth(i)) << "Formal parameter " << formals->nth(i)->get_name() << " is multiply defined." << std::endl;
    else if (classMap.find(formals->nth(i)->get_type()) == classMap.end())
      semant_error(formals->nth(i)) << "Class " << formals->nth(i)->get_type() << " of formal parameter " 
				    << formals->nth(i)->get_name() << " is undefined." << std::endl;
    else 
      objects.addid(formals->nth(i)->get_name(), new Symbol(formals->nth(i)->get_type()));
  }

  Symbol expr_type = expr->typecheck();
  if (return_type != SELF_TYPE && classMap.find(return_type) == classMap.end())
    semant_error(this) << "Undefined return type " << return_type << " in method " << name << "." << std::endl;
  else if (!subtype(expr_type, return_type))
    semant_error(this) << "Inferred return type " << expr_type << " of method " << name 
		       << " does not conform to declared return type " << return_type << "." << std::endl;
  objects.exitscope();
}

Symbol branch_class::typecheck() {
  objects.enterscope();
  objects.addid(name, new Symbol(type_decl));
  type = expr->typecheck();
  objects.exitscope();
  return type;
}

Symbol assign_class::typecheck() {
  Symbol rtype = expr->typecheck();

  if (objects.lookup(name) == NULL) {
    semant_error(this) << "Assignment to undeclared variable " << name << "." << std::endl;
    type = rtype;
    return type;
  }
  
  Symbol ltype = *objects.lookup(name);
  
  if (!subtype(rtype, ltype)) {
    semant_error(this) << "Type " << rtype << " of assigned expression does not conform to declared type " 
		       << ltype << " of identifier " << name << "." << std::endl;
    type = ltype;
    return type;
  }

  type = rtype;
  return type;
}

Symbol static_dispatch_class::typecheck() {
  bool error = false;
  Symbol expr_type = expr->typecheck();

  if (type_name != SELF_TYPE && classMap.find(type_name) == classMap.end()) {
    semant_error(this) << "Static dispatch to undefined class " << type_name << "." << std::endl;
    type = Object;
    return type;
  }

  if (expr_type != SELF_TYPE && classMap.find(expr_type) == classMap.end()) {
    type = Object;
    return type;
  }

  if (!subtype(expr_type, type_name)) {
    error = true;
    semant_error(this) << "Expression type " << expr_type << " does not conform to declared static type " 
			 << type_name << "." << std::endl;
  }

  method_class* method = NULL;
  std::vector<Class_> chain = getInheritanceChain(type_name);
  for (std::size_t i = 0; i < chain.size(); i++) {
    std::vector<method_class*> _methods = methods[chain[i]];
    for (std::size_t j = 0; j < _methods.size(); j++) {
      if (_methods[j]->get_name() == name) {
	method = _methods[j];
	break;
      }
    }
  }

  if (method == NULL) {
    error = true;
    semant_error(this) << "Static dispatch to undefined method " << name << "." << std::endl;
  } else {
    Formals formals = method->get_formals();
    int k1 = actual->first(), k2 = formals->first();
    while (actual->more(k1) && formals->more(k2)) {
      Symbol actual_type = actual->nth(k1)->typecheck();
      Symbol formal_type = formals->nth(k2)->get_type();
      if (!subtype(actual_type, formal_type)) {
	error = true;
	semant_error(this) << "In call method " << name << ", type " << actual_type 
			   << " of parameter " << formals->nth(k2)->get_name() 
			   << " does not conform to declared type " << formal_type << "." << std::endl;
      }
      k1 = actual->next(k1);
      k2 = formals->next(k2);
      if (actual->more(k1) ^ formals->more(k2)) {
	error = true;
	semant_error(this) << "Method " << name << " called with wrong number of arguments." << std::endl;
      }
    }
  }

  if (error) {
    type = Object;
  } else {
    type = method->get_return_type();
    if (type == SELF_TYPE)
      type = expr_type;
  }

  return type;
}

Symbol dispatch_class::typecheck() {
  bool error = false;
  Symbol expr_type = expr->typecheck();

  if (expr_type != SELF_TYPE && classMap.find(expr_type) == classMap.end()) {
    semant_error(this) << "Dispatch on undefined class " << expr_type << "." << std::endl;
    type = Object;
    return type;
  }

  method_class* method = NULL;
  std::vector<Class_> chain = getInheritanceChain(expr_type);
  for (std::size_t i = 0; i < chain.size(); i++) {
    std::vector<method_class*> _methods = methods[chain[i]];
    for (std::size_t j = 0; j < _methods.size(); j++) {
      if (_methods[j]->get_name() == name) {
	method = _methods[j];
	break;
      }
    }
  }

  if (method == NULL) {
    semant_error(this) << "Dispatch to undefined method " << name << "." << std::endl;
    error = true;
  } else {
    Formals formals = method->get_formals();
    int k1 = actual->first(), k2 = formals->first();
    while (actual->more(k1) && formals->more(k2)) {
      Symbol actual_type = actual->nth(k1)->typecheck();
      Symbol formal_type = formals->nth(k2)->get_type();
      if (!subtype(actual_type, formal_type)) {
	error = true;
	semant_error(this) << "In call method " << name << ", type " << actual_type 
			   << " of parameter " << formals->nth(k2)->get_name() 
			   << " does not conform to declared type " << formal_type << "." << std::endl;
      }
      k1 = actual->next(k1);
      k2 = formals->next(k2);
      if (actual->more(k1) ^ formals->more(k2)) {
	error = true;
	semant_error(this) << "Method " << name << " called with wrong number of arguments." << std::endl;
      }
    }
  }

  if (error) {
    type = Object;
  } else {
    type = method->get_return_type();
    if (type == SELF_TYPE)
      type = expr_type;
  }

  return type;
}

Symbol cond_class::typecheck() {
  if (pred->typecheck() != Bool)
    semant_error(this) << "Predicate of 'if' does not have type Bool." << std::endl;

  Symbol then_type = then_exp->typecheck();
  Symbol else_type = else_exp->typecheck();

  if (then_type == SELF_TYPE && else_type == SELF_TYPE)
    type = SELF_TYPE;
  else
    type = lca(then_type, else_type)->get_name();
  return type;
}

Symbol loop_class::typecheck() {
  if (pred->typecheck() != Bool) {
    semant_error(this) << "Loop condition does not have type Bool." << std::endl;
  }

  body->typecheck();
  type = Object;
  return type;
}

Symbol typcase_class::typecheck() {
  Symbol expr_type = expr->typecheck();
  std::set<Symbol> branch_type_decls;
  
  for (int i = cases->first(); cases->more(i); i = cases->next(i)) {
    branch_class* branch = static_cast<branch_class*>(cases->nth(i));
    if (branch_type_decls.find(branch->get_type_decl()) != branch_type_decls.end())
      semant_error(branch) << "Duplicate branch " << branch->get_type_decl() << " in case statement." << std::endl;
    else
      branch_type_decls.insert(branch->get_type_decl());
    Symbol branch_type = branch->typecheck();
    if (i == cases->first())
      type = branch_type;
    else if (type != SELF_TYPE || branch_type != SELF_TYPE)
      type = (lca(type, branch_type))->get_name();
  }

  return type;
}

Symbol block_class::typecheck() {
  for (int i = body->first(); body->more(i); i = body->next(i))
    type = body->nth(i)->typecheck();
  return type;
}

Symbol let_class::typecheck() {
  //std::cerr << "Checking let expr." << std::endl;
  objects.enterscope();
  objects.addid(identifier, new Symbol(type_decl));
  
  if (identifier == self) {
    semant_error(this) << "'self' cannot be bound in a 'let' expression." << std::endl;
  }

  Symbol init_type = init->typecheck();

  if (type_decl != SELF_TYPE && classMap.find(type_decl) == classMap.end()) {
    semant_error(this) << "Class " << type_decl 
		       << " of let-bound identifier " << identifier 
		       << " is undefined." << std::endl;
  } else if (init_type != No_type && !subtype(init_type, type_decl)) {
    semant_error(this) << "Inferred type " << init_type << " of initialization of " 
		       << identifier << " does not conform to identifier's declared type " << type_decl << std::endl;
  }

  type = body->typecheck();
  objects.exitscope();
  //std::cerr << "Checked let expr. " << std::endl;
  return type;
}

Symbol plus_class::typecheck() {
  Symbol type1 = e1->typecheck();
  Symbol type2 = e2->typecheck();
  if (type1 != Int || type2 != Int) {
    semant_error(this) << "non-Int arguments: " << type1 
		       << " + " << type2 << std::endl;
  }
  type = Int;
  return type;
}

Symbol sub_class::typecheck() {
  Symbol type1 = e1->typecheck();
  Symbol type2 = e2->typecheck();
  if (type1 != Int || type2 != Int) {
    semant_error(this) << "non-Int arguments: " << type1 
		       << " - " << type2 << std::endl;
  }
  type = Int;
  return type;
}

Symbol mul_class::typecheck() {
  Symbol type1 = e1->typecheck();
  Symbol type2 = e2->typecheck();
  if (type1 != Int || type2 != Int) {
    semant_error(this) << "non-Int arguments: " << type1 
		       << " * " << type2 << std::endl;
  }
  type = Int;
  return type;
}

Symbol divide_class::typecheck() {
  Symbol type1 = e1->typecheck();
  Symbol type2 = e2->typecheck();
  if (type1 != Int || type2 != Int) {
    semant_error(this) << "non-Int arguments: " << type1 
		       << " / " << type2 << std::endl;
  }
  type = Int;
  return type;
}

Symbol neg_class::typecheck() {
  type = e1->typecheck();
  if (type != Int) {
    semant_error(this) << "Argument of '~' has type " << type
		       << " instead of Int." << std::endl;
  }
  type = Int;
  return type;
}

Symbol lt_class::typecheck() {
  Symbol type1 = e1->typecheck();
  Symbol type2 = e2->typecheck();
  if (type1 != Int || type2 != Int) {
    semant_error(this) << "non-Int arguments: " << type1 
		       << " < " << type2 << std::endl;
  }
  type = Bool;
  return type;
}

Symbol eq_class::typecheck() {
  Symbol type1 = e1->typecheck();
  Symbol type2 = e2->typecheck();
  if ((isBasic(type1) || isBasic(type2)) && type1 != type2) {
    semant_error(this) << "Illegal comparison with a basic type." << std::endl;
  }
  type = Bool;
  return type;
}

Symbol leq_class::typecheck() {
  Symbol type1 = e1->typecheck();
  Symbol type2 = e2->typecheck();
  if (type1 != Int || type2 != Int) {
    semant_error(this) << "non-Int arguments: " << type1 
		       << " <= " << type2 << std::endl;
  }
  type = Bool;
  return type;
}

Symbol comp_class::typecheck() {
  type = e1->typecheck();
  if (type != Bool) {
    semant_error(this) << "Argument of 'not' has type " 
		       << type << " instead of Bool." << std::endl;
  }
  type = Bool;
  return type;
}

Symbol int_const_class::typecheck() {
  type = Int;
  return type;
}

Symbol bool_const_class::typecheck() {
  type = Bool;
  return type;
}

Symbol string_const_class::typecheck() {
  type = Str;
  return type;
}

Symbol new__class::typecheck() {
  if (type_name != SELF_TYPE && classMap.find(type_name) == classMap.end()) {
    semant_error(this) << "'new' used with undefined class " 
		       << type_name << "." << std::endl;
    type_name = Object;
  }
  type = type_name;
  return type;
}

Symbol isvoid_class::typecheck() {
  e1->typecheck();
  type = Bool;
  return type;
}

Symbol no_expr_class::typecheck() {
  type = No_type;
  return type;
}

Symbol object_class::typecheck() {
  if (name == self) {
    type = SELF_TYPE;
  } else if (objects.lookup(name)) {
    type = *objects.lookup(name);
  } else {
    semant_error(this) << "Undeclared identifier " << name << "." << std::endl;
    type = Object;
  }
  return type;
}

const char* method_class::constructor() {
  return "method";
}

const char* attr_class::constructor() {
  return "attr";
}

/*   This is the entry point to the semantic checker.

     Your checker should do the following two things:

     1) Check that the program is semantically correct
     2) Decorate the abstract syntax tree with type information
        by setting the `type' field in each Expression node.
        (see `tree.h')

     You are free to first do 1), make sure you catch all semantic
     errors. Part 2) can be done in a second stage, when you want
     to build mycoolc.
 */
void program_class::semant() {
    initialize_constants();

    // ClassTable constructor may do some semantic analysis
    GraphChecker checker(classes);

    if (errors()) {
      std::cerr << "Compilation halted due to static semantic errors." << std::endl;
      exit(1);
    }

    typecheck(classes);
    
    if (errors()) {
      std::cerr << "Compilation halted due to static semantic errors." << std::endl;
      exit(1);
    }
}


