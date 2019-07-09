/*
 *  The scanner definition for COOL.
 */

/*
 *  Stuff enclosed in %{ %} in the first section is copied verbatim to the
 *  output, so headers and global definitions are placed here to be visible
 * to the code in the file.  Don't remove anything that was here initially
 */
%{
#include <cool-parse.h>
#include <stringtab.h>
#include <utilities.h>

/* The compiler assumes these identifiers. */
#define yylval cool_yylval
#define yylex  cool_yylex

/* Max size of string constants */
#define MAX_STR_CONST 1025
#define YY_NO_UNPUT   /* keep g++ happy */

extern FILE *fin; /* we read from this file */

/* define YY_INPUT so we read from the FILE fin:
 * This change makes it possible to use this scanner in
 * the Cool compiler.
 */
#undef YY_INPUT
#define YY_INPUT(buf,result,max_size) \
	if ( (result = fread( (char*)buf, sizeof(char), max_size, fin)) < 0) \
		YY_FATAL_ERROR( "read() in flex scanner failed");

char string_buf[MAX_STR_CONST]; /* to assemble string constants */
char *string_buf_ptr;

extern int curr_lineno;
extern int verbose_flag;

extern YYSTYPE cool_yylval;

/*
 *  Add Your own definitions here
 */

bool is_valid(int len) {
  if (len >= MAX_STR_CONST) {
     cool_yylval.error_msg = "String constant too long";
     return false;	
  }

  return true;
}

int open = 0, str_len = 0;
int str_count = 0, int_count = 0, id_count = 0;

%}

/*
 * Define names for regular expressions here.
 */

DARROW          =>
ASSIGN		<-
LE		<=

a		[aA]
b		[bB]
c		[cC]
d		[dD]
e		[eE]
f		[fF]
g		[gG]
h		[hH]
i		[iI]
j		[jJ]
k		[kK]
l		[lL]
m		[mM]
n		[nN]
o		[oO]
p		[pP]
q		[qQ]
r		[rR]
s		[sS]
t		[tT]
u		[uU]
v		[vV]
w		[wW]
x		[xX]
y		[yY]
z		[zZ]

%x NESTED_COMMENT LINE_COMMENT IN_STRING STRING_ERROR

%%

\n			curr_lineno++;
[ \f\r\t\v]+		;

 /*
  *  Nested comments
  */
<NESTED_COMMENT>"(*"	open++;
<NESTED_COMMENT>"*)"	{ open--;
			  if (open == 0) {
			     BEGIN(INITIAL);
			  }
			}

"(*"		{ open = 1; BEGIN(NESTED_COMMENT); }
"*)"		{ cool_yylval.error_msg = "Unmatched *)"; return ERROR; }

<NESTED_COMMENT><<EOF>>	{ cool_yylval.error_msg = "EOF in comment"; 
			  BEGIN(INITIAL); 
			  return ERROR; 
			}

<NESTED_COMMENT>\\.+	;
<NESTED_COMMENT>\*	;
<NESTED_COMMENT>\(	;
<NESTED_COMMENT>[^\n\*\(]+	;
<NESTED_COMMENT>\n	curr_lineno++;

"--"			BEGIN(LINE_COMMENT);

<LINE_COMMENT>\n	{ curr_lineno++; BEGIN(INITIAL); }
<LINE_COMMENT>.		;

 /*
  *  String constants
  */

<IN_STRING>\"		{ cool_yylval.symbol = new StringEntry(strdup(string_buf), str_len, str_count++); BEGIN(INITIAL); return STR_CONST; }

<IN_STRING>\\n 		{ if (!is_valid(str_len)) {
			     BEGIN(STRING_ERROR);
			     return ERROR;
			  }
			  string_buf[str_len++] = '\n';
			}

<IN_STRING>\\f		{ if (!is_valid(str_len)) {
			     BEGIN(STRING_ERROR);
			     return ERROR;
			  }
			  string_buf[str_len++] = '\f';
			}

<IN_STRING>\\t		{ str_len++;
 			  if (!is_valid(str_len)) {
			     BEGIN(STRING_ERROR);
			     return ERROR;
			  }
			  string_buf[str_len++] = '\t';
			}

<IN_STRING>\\b		{ if (!is_valid(str_len)) {
			     BEGIN(STRING_ERROR);
			     return ERROR;
			  }
			  string_buf[str_len++] = '\b';
			} 

<STRING_ERROR>\"	BEGIN(INITIAL);
<STRING_ERROR>\n	{ curr_lineno++; BEGIN(INITIAL); }
<STRING_ERROR>.		;

\"			{ str_len = 0; BEGIN(IN_STRING); }


<IN_STRING>\\.		{ if (!is_valid(str_len)) {
			     BEGIN(STRING_ERROR);
			     return ERROR;
			  }
			  string_buf[str_len++] = yytext[1];
			}

<IN_STRING>\n		{ cool_yylval.error_msg = "Unterminated string constant"; curr_lineno++; BEGIN(INITIAL); return ERROR; }

<IN_STRING><<EOF>>	{ cool_yylval.error_msg = "EOF in string constant"; BEGIN(STRING_ERROR); return ERROR; }

<IN_STRING>\0		{ cool_yylval.error_msg = "String contains null character"; BEGIN(STRING_ERROR); return ERROR; }


<IN_STRING>.		{ if (!is_valid(str_len)) {
			     BEGIN(STRING_ERROR);
			     return ERROR;
			  }
			  string_buf[str_len++] = yytext[0];
			}

 /*
  *  The multiple-character operators.
  */
{DARROW}		{ return (DARROW); }
{ASSIGN}		{ return (ASSIGN); }
{LE}			{ return (LE); }

 /*
  * Keywords are case-insensitive except for the values true and false,
  * which must begin with a lower-case letter.
  */
{c}{l}{a}{s}{s}		{ return (CLASS); }
{e}{l}{s}{e}		{ return (ELSE); }
{f}{i}			{ return (FI); }
{i}{f}			{ return (IF); }
{i}{n}			{ return (IN); }
{i}{n}{h}{e}{r}{i}{t}{s}  { return (INHERITS); }
{l}{e}{t}		{ return (LET); }
{l}{o}{o}{p}		{ return (LOOP); }
{p}{o}{o}{l}		{ return (POOL); }
{t}{h}{e}{n}		{ return (THEN); }
{w}{h}{i}{l}{e}		{ return (WHILE); }
{c}{a}{s}{e}		{ return (CASE); }
{e}{s}{a}{c}		{ return (ESAC); }
{o}{f}			{ return (OF); }
{n}{e}{w}		{ return (NEW);	}
{i}{s}{v}{o}{i}{d}	{ return (ISVOID); }
{n}{o}{t}		{ return (NOT); }

t{r}{u}{e}		{ cool_yylval.boolean = true; return (BOOL_CONST); }
f{a}{l}{s}{e}		{ cool_yylval.boolean = false; return (BOOL_CONST); }

[0-9]+			{ cool_yylval.symbol = new IntEntry(strdup(yytext), strlen(yytext), int_count++); return (INT_CONST); }
[A-Z][a-zA-Z0-9_]*	{ cool_yylval.symbol = new IdEntry(strdup(yytext), strlen(yytext), id_count++); return (TYPEID); }
[a-z][a-zA-Z0-9_]*	{ cool_yylval.symbol = new IdEntry(strdup(yytext), strlen(yytext), id_count++); return (OBJECTID); } 


 /*
  *  String constants (C syntax)
  *  Escape sequence \c is accepted for all characters c. Except for 
  *  \n \t \b \f, the result is c.
  *
  */
"{"			{ return '{'; }
"}"			{ return '}'; }
";"			{ return ';'; }
":"			{ return ':'; }
","			{ return ','; }
"."			{ return '.'; }
"("			{ return '('; }
")"			{ return ')'; }
"@"			{ return '@'; }
"+"			{ return '+'; }
"-"			{ return '-'; }
"*"			{ return '*'; }
"/"			{ return '/'; }
"~"			{ return '~'; }
"<"			{ return '<'; }
"="			{ return '='; }

.			{ cool_yylval.error_msg = strdup(yytext); return ERROR; }

%%
