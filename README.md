# COOL 
This repository contains partial implementation of COOL compiler.
COOL is statically-typed object-oriented expression-based language.
Following parts are implemented (main files showed below):

- [X] Lexer - maps character stream to token stram using set of regular expression; implemented using `flex` tool.
  - `assignments/Lexer/cool.flex` - main file, which contains set of regexps for lexing
- [X] Parser - builds AST (Abstract Syntax Tree) from given stream of tokens. Parser is implemented using `bison`, which uses shift-reduce parsing.
  - `assignments/Parser/cool.y` - contains the specification of the Context Free Grammar of COOL language
- [X] Type checker - performs compile-time checking of program for correctness. Checks for undeclared variables, types, including class attributes, formal parameters and let expressions. Builds inheritance graph and performs cycle checks.
  - `assignments/Typechecker/cool-tree.h` - declarations of AST parts
  - `assignments/Typechecker/semant.h` - declarations and definitions for semantic analyzer
  - `assignments/Typechecker/semant.cc` - main file with method definitions
- [ ] Assembly code generator
