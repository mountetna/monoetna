
# Vulcan

Archimedes is a "calculation service". It runs code in the browser in a
scripting language called Vulcan.

Currently "Vulcan" runs in an interpreter written in Ruby. This means Vulcan
borrows heavily from Ruby for a lot of its capabilities (mostly the type
system). Consequently Vulcan should not be thought of as a "programming
language", meant to produce complete programs, but rather as a "query
language" meant to return data.

Instead of using our own language, we might run python or R code. However, by
all accounts sandboxing Python is hard, so any code we run would risk people
breaking out of the sandbox. This is less likely to happen with the scripting
language since code has to make it through the lexer/parser, and it's unlikely
people could thus break out of the sandbox.

Here is what is defined in the language already:

## Types

The basic types:
  String - delimited by single quotes, which may be escaped by double-quoting (e.g. 'That''s all folks!')
  Number - a float
  true, false - boolean values
  nil - an empty value
  Vector - somewhere between an array and an associative array (Hash), this is
  a list of values with optional string labels. Labels may be strings
  ('label') or an identifier (label)
  Matrix - an array of arrays, with data in rows and columns, along with
  string row names and column names.

## Variables

Data may be bound to variables. There are two kinds of
variable, exported variables (@blah) and local variables
($blah).

## Operators

### Assignment

`=` -  Valid l-values for assignments are @export or $local variables.

`@name = 'Hercules'`

### Binary operators

Most mathematical operations `^, /, *, +, -, %` will map across vectors (i.e., @x == 2 will return a vector like [ @i == 2 for each @i in @x ]); vectors of unequal length will cause an error.

```
@product = 2 * 3
# @product => 6

@product = [ 1, 2, 3 ] * 3
# @product => [ 3, 6, 9 ]

@product = [ 1, 2, 3 ] * [ 4, 5, 6 ]
# @product => [ 4, 10, 18 ]
```

Comparison operators `==, !=, <, >, <=, >=` likewise span across vectors.

```
@comparison = 2 > 3
# @comparison => false

@comparison = [ 1, 2, 3 ] > 2
# @comparison => [ false, false, true ]

@comparison = [ 1, 2, 3 ] > [ 4, 5, 6]
# @comparison => [ false, false, false ]
```

### Unary operations

Negation `-, !` of numbers and booleans also apply to vectors.

```
@y = 4
@x = -@y
# @x => -4

@y = [ 1, 2, 3]
@x = -@y
# @x => [ -1, -2, -3 ]

@y = false
@x = !@y
# @x => true

@y = [ false, false, true ]
@x = !@y
# @x => [ true, true, false ]
```

### Access operators

The dollar operator `$` is a way to dereference Matrix columns or Vector elements
```
@x = [ ant: 1, bear: 2, cat: 3 ]
@y = @x$bear
# @y => 2

@x = [ ant: [ 1, 2, 3 ], bear: [ 4, 5, 6 ], cat: [ 7, 8, 9 ] ]
```

Brackets `[]` may dereference vector or matrix contents.

Vectors may be dereferenced using an integer, a vector of
integers, or a (label) string.

```
@x = [ ant: 1, bear: 2, cat: 3]
@y = @x[1]
# @y => 2

@y = @x[ [ 1, 2 ] ]
# @y => [ bear: 2, cat: 3 ]

@y = @x[ 'bear' ]
# @y => 2
```

You may reference a matrix column using brackets (similar to `$`):

```
@matrix = ::[
  ant: [ ice: 1, jump: 2, kid: 3 ],
  bear: [ 4, 5, 6 ],
  cat: [ 7, 8, 9 ]
]
@x = @matrix[ 'cat' ]
# @x => [ 7, 8, 9 ]
```

You may also slice the matrix by rows and columns, often useful
in conjunction with the `which` function:

```
@x = @matrix[ [1, 2], [ 1, 2 ] ]

# @x => ::[ bear: [ 5, 6 ], cat: [ 8, 9 ] ]
```

### Vector operators

The concatenation operator lets you join two vectors end-to-end:

```
@x = [ 1, 2, 3 ] .. [ 4, 5, 6 ]
# @x => [ 1, 2, 3, 4, 5, 6 ]
```

The transpose operator lets you treat a Vector as a column
for matrix math:

```
@x = :: [ [ 1, 2 ], [1, 2 ] ]
@y = @x + [ 2, 1 ]
@z = @x + ~[ 2, 1 ]

# @y => [ 3 2 ]
#       [ 4 3 ]
# @z => [ 3 3 ]
#       [ 3 3 ]
```

The map operator `%%` lets you iterate across the elements of a vector

```
@v = [ 1, 2, 3 ] %% ($value, $i, $label) => $value * $i
# @v => [ 0, 2, 6 ]
```

The reduce operator `>>` lets you reduce the elements of a vector to a value

```
@v = [ 1, 2, 3 ] >> ($s = 0, $value, $i, $label) => $s + $v
# @v => 6

@v = [ 'ant', 'bear', 'cat' ] >> ($s = '', $v) => $s + $v
# @v => 'antbearcat'
```

### Matrix operators

The bind operator lets you transform a vector of column vectors
into a matrix:
```
@x = :: [ a: [ i: 1, j: 2 ], b: [ 3, 4 ] ]
# @x =>     a b
#       i [ 1 3 ]
#       j [ 2 4 ]
```

You may also bind to an existing matrix:
```
@y = @x :: [ c: [ 5, 6 ] ]
# @x =>     a b c
#       i [ 1 3 5 ]
#       j [ 2 4 6 ]
```

Or you may bind in place:

```
@y ::= [ d: [ 7, 8 ] ]
# @y =>     a b c d
#       i [ 1 3 5 7 ]
#       j [ 2 4 6 8 ]
```

The transpose operator lets you transpose a matrix:

```
@z = ~ @x
# @z =>     i j
#       a [ 1 2 ]
#       b [ 3 4 ]
```

The map operator allows you to iterate across the rows of a matrix.  It takes a
function as its second argument:
```
@x = ::[ a: [ i: 1, j: 2 ], b: [ 3, 4 ] ]
@y = @x %% ($row, $row_number) => $row * 2 + $row_number
# y =>     a b
       i [ 2 4 ]
       j [ 4 9 ]
```

## Functions

A function declaration may be bound to a variable:

```
$multiply = ($arg1, $arg2) => {
  return $arg1 * $arg2
}
# or
$multiply = ($arg1, $arg2) => $arg1 * $arg2

@product = $multiply(2, 3)
# @product => 6
```

The function has local scope, and read-only access to the parent scope. Thus,
you may do this:

```
$remultiply = ($arg1, $arg2) => {
  $x = $multiply($arg1, $arg2)
  return $x * $x
}
```

but not:

```
@y = 5
$assign_y = ($arg1, $arg2) => {
  @y = $arg1 * $arg2
}
```

You may make default assignments to the arguments:
```
@rigged_multiply = ($x = 2, $y = 4) => $x * $y

@a = @rigged_multiply(3)
# @a => 12

@b = @rigged_multiply()
# @b => 8

@c = @rigged_multiply(1, 1)
# @c => 1
```

## Macros

A macro declaration is bound to a variable.  Macros are invoked
with string arguments, which are pasted into the macro
definition and then evaluated.

```
$macro = { $arg1 * %1 }

$arg1 = 2
$arg2 = 3
@product = $macro( '$arg2' )
```


## Conditionals

The ternary operator allows conditional expressions:

```
$s = 1
$y = $s == 1 ? 'red' : 'green'
# $y => 'red'
```

You can also use if / else statements:

```
if ($s == 1)
  $y = 2
else
  $y = 'red'

# $y => 2

if ($s == 1 ) {
  @x = 2
  @y = 3
} else {
  @x = 4
  @y = 6
}

# @x => 2
# @y => 3
```
