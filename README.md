# pyvdrm
[![Build
Status](https://travis-ci.org/jeff-k/pyvdrm.svg?branch=master)](https://travis-ci.org/jeff-k/pyvdrm)

Tools for interpreting drug resistance mutations in viral amino acid sequences

A DRM calling algorithm's `Rule`s should be specified in a context free grammar. It is the intention of this
package to supply the semantics for the most commonly used operators on sets of
`Mutation`s, but these can be overloaded.

For example, the syntax of the Stanford HIVdb ASI2 algorithm is specified in Backus-Naur
form:

```
...
listitems -> comma residue
selectlist -> residue listitems*
...
scoreitem booleancondition mapper min? number | max l_par scorelist r_par
```

The semantics for these operators is implemented in the `drm` module.

With these syntax and semantics, we can define an _algorithm_ as a set of rules
that score a list of mutations.

## Supported Operations

There are distinct operations one may want to apply to a rule:

  - **Validation** Malformed rules should be detected before the mutations are 
    evaluated
  - **Printing/Transformation** Displaying rules with or without arguments is
    useful for debugging and reporting but additionally for translating a set of
    rules into other machine readable formats independent of semantics
  - **Evaluation** Given an set of positional mutations, a rule returns a score

## API Synopsis

  - **Algorithm**: the domain specific implementation of a rule grammar
  - **Rule**: an instance of the algorithm, initialized with
  - **Mutation**: the `vcf` module supports 
  - **Score**: the result of evaluating a rule with a set of mutations

An _environment_ is a set of Mutations on which a Rule is applied.

## Implementing a DRM algorithm

1. Define the rule syntax with pyparsing
2. Implement the semantics for the rule's operators
3. Define test cases

## Application Notes

The three application domains are:

  - The clinician that evaluates mutations to determine drug resistance
  - The bioinformatician that determines drug resistant mutations
  - The data scientist that derives DRM algorithms

The intention of this package is to facilitate use of non-ambiguous human
readable drug resistance rules but the rules are also machine readable. For
example, a machine learning algorithm should be able to derive these rules and
close the gap between research and clinical use.
