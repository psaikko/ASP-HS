# ASP-HS : A solver for optimization in answer set programming

ASP-HS uses a hybrid SAT-IP implicit hitting set approach to solve ASP optimization problems.
The solver works on variable-free programs so a grounder (e.g. gringo) is required.

ASP-HS uses Wasp (included) as its underlying ASP core extractor and IBM CPLEX as its hitting set optimizer. 
CPLEX must be installed separately before compiling ASP-HS. A free academic license for CPLEX can be obtained through the IBM Academic Initiative. ASP-HS has been tested with CPLEX 12.7.

Run the `configure.py` script before compiling ASP-HS to set the necessary IP solver filepaths.

To compile
```
make clean && make
```

A small test suite can be run with
```
make test
```

ASP-HS can read a ground ASP instance from file as
```
./bin/asp-hs file
```
or through the standard input stream when given parameter `-`
```
./bin/asp-hs - <file
```

For information on solver parameters, run
```
./bin/asp-hs --help
```
