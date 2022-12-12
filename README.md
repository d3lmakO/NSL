# NSL
Exercises for the Numerical Simulation Laboratory course (year: 2021/2022).


## Compilation

To compile a single exercise run:
```bash
make ex*
```
with * a number between 1 and 10.

To compile all exercises run:
```bash
make all
```
If you want to clean object files and executable files run:
```bash
make clean
```
or:
```bash
make clean_ex*
```

## Execution

To run a single exercise follow this (example with ex9):
```bash
cd ex9/
./9_1 
```
## WARNING

- In ex8 data file `8_2_psi2.txt` not present beacuse it exceeded 100MB (You need to run the exercise).
- In ex10 you need to include `"mpi.h"`.

## Visualization

Please to visualize the final data analysis for all exercises run
```bash
jupyter-lab
```
before the ex* directories.
