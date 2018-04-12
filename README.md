# Matrix Compute Engine
A simple matrix compute engine, optimised by Intel SSE instructions. 

## Install and Compile
To install the project, simply clone it first
```
git clone https://github.com/ShuminKong/MatrixComputeEngine.git
```
And then navigate to the directory of this project to compile it.
```
cd MatrixComputeEngine
make all
```
## Usage
Run the program with parameters passed as follow.
```
./matrix <length of matrix> <maximum number of threads>
```
And then, enter the commands as follow. 
```
SET <key> = identity
SET <key> = random <seed>
SET <key> = uniform <value>
SET <key> = sequence <start> <step>
SET <key> = cloned <matrix>
SET <key> = sorted <matrix>
SET <key> = rotated <matrix>
SET <key> = reversed <matrix>
SET <key> = transposed <matrix>
SET <key> = matrix.add <matrix a> <matrix b>
SET <key> = matrix.mul <matrix a> <matrix b>
SET <key> = matrix.pow <matrix> <exponent>
SET <key> = matrix.conv <matrix> <kernel>
SET <key> = scalar.add <matrix> <scalar>
SET <key> = scalar.mul <matrix> <scalar>
SHOW <key>
SHOW <key> row <number>
SHOW <key> column <number>
SHOW <key> element <row> <column>
COMPUTE sum <key>
COMPUTE trace <key>
COMPUTE minimum <key>
COMPUTE maximum <key>
COMPUTE determinant <key>
COMPUTE frequency <key> <value>
```
## Acknowledgement
I would like to ackowledge University of Sydney for providing the ```main.c``` and ```matrix.h``` file, to which the license declared in LICENSE file does not apply.
