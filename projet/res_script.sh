#!/bin/bash



# for FILE in res/h/128/*; do
# 	diff $FILE res/h/128/resh128s.txt >> res/h/128/diff/diff.txt
# done

echo 128

./src/sand 0 c 128 s > res/c/128/s.txt
./src/sand 0 c 128 S > res/c/128/S.txt
./src/sand 0 c 128 u > res/c/128/u.txt
./src/sand 0 c 128 U > res/c/128/U.txt

./src/sand 0 h 128 s > res/h/128/s.txt
./src/sand 0 h 128 S > res/h/128/S.txt
./src/sand 0 h 128 u > res/h/128/u.txt
./src/sand 0 h 128 U > res/h/128/U.txt

echo 512

./src/sand 0 c 512 s > res/c/512/s.txt
./src/sand 0 c 512 S > res/c/512/S.txt
./src/sand 0 c 512 u > res/c/512/u.txt
./src/sand 0 c 512 U > res/c/512/U.txt

./src/sand 0 h 512 s > res/h/512/s.txt
./src/sand 0 h 512 S > res/h/512/S.txt
./src/sand 0 h 512 u > res/h/512/u.txt
./src/sand 0 h 512 U > res/h/512/U.txt

echo parallel

OMP_NUM_THREADS=4 ./src/sand 0 c 128 F > res/c/128/F4.txt
OMP_NUM_THREADS=4 ./src/sand 0 c 128 P > res/c/128/P4.txt
OMP_NUM_THREADS=4 ./src/sand 0 c 128 t > res/c/128/t4.txt

OMP_NUM_THREADS=4 ./src/sand 0 h 128 F > res/h/128/F4.txt
OMP_NUM_THREADS=4 ./src/sand 0 h 128 P > res/h/128/P4.txt
OMP_NUM_THREADS=4 ./src/sand 0 h 128 t > res/h/128/t4.txt

OMP_NUM_THREADS=4 ./src/sand 0 c 512 F > res/c/512/F4.txt
OMP_NUM_THREADS=4 ./src/sand 0 c 512 P > res/c/512/P4.txt
OMP_NUM_THREADS=4 ./src/sand 0 c 512 t > res/c/512/t4.txt

OMP_NUM_THREADS=4 ./src/sand 0 h 512 F > res/h/512/F4.txt
OMP_NUM_THREADS=4 ./src/sand 0 h 512 P > res/h/512/P4.txt
OMP_NUM_THREADS=4 ./src/sand 0 h 512 t > res/h/512/t4.txt

echo parallel8

OMP_NUM_THREADS=8 ./src/sand 0 c 128 F > res/c/128/F8.txt
OMP_NUM_THREADS=8 ./src/sand 0 c 128 P > res/c/128/P8.txt
OMP_NUM_THREADS=8 ./src/sand 0 c 128 t > res/c/128/t8.txt

OMP_NUM_THREADS=8 ./src/sand 0 h 128 F > res/h/128/F8.txt
OMP_NUM_THREADS=8 ./src/sand 0 h 128 P > res/h/128/P8.txt
OMP_NUM_THREADS=8 ./src/sand 0 h 128 t > res/h/128/t8.txt

OMP_NUM_THREADS=8 ./src/sand 0 c 512 F > res/c/512/F8.txt
OMP_NUM_THREADS=8 ./src/sand 0 c 512 P > res/c/512/P8.txt
OMP_NUM_THREADS=8 ./src/sand 0 c 512 t > res/c/512/t8.txt

OMP_NUM_THREADS=8 ./src/sand 0 h 512 F > res/h/512/F8.txt
OMP_NUM_THREADS=8 ./src/sand 0 h 512 P > res/h/512/P8.txt
OMP_NUM_THREADS=8 ./src/sand 0 h 512 t > res/h/512/t8.txt

echo parallel16

OMP_NUM_THREADS=16 ./src/sand 0 c 128 F > res/c/128/F16.txt
OMP_NUM_THREADS=16 ./src/sand 0 c 128 P > res/c/128/P16.txt
OMP_NUM_THREADS=16 ./src/sand 0 c 128 t > res/c/128/t16.txt

OMP_NUM_THREADS=16 ./src/sand 0 h 128 F > res/h/128/F16.txt
OMP_NUM_THREADS=16 ./src/sand 0 h 128 P > res/h/128/P16.txt
OMP_NUM_THREADS=16 ./src/sand 0 h 128 t > res/h/128/t16.txt

OMP_NUM_THREADS=16 ./src/sand 0 c 512 F > res/c/512/F16.txt
OMP_NUM_THREADS=16 ./src/sand 0 c 512 P > res/c/512/P16.txt
OMP_NUM_THREADS=16 ./src/sand 0 c 512 t > res/c/512/t16.txt

OMP_NUM_THREADS=16 ./src/sand 0 h 512 F > res/h/512/F16.txt
OMP_NUM_THREADS=16 ./src/sand 0 h 512 P > res/h/512/P16.txt
OMP_NUM_THREADS=16 ./src/sand 0 h 512 t > res/h/512/t16.txt

echo parallel24

OMP_NUM_THREADS=24 ./src/sand 0 c 128 F > res/c/128/F24.txt
OMP_NUM_THREADS=24 ./src/sand 0 c 128 P > res/c/128/P24.txt
OMP_NUM_THREADS=24 ./src/sand 0 c 128 t > res/c/128/t24.txt

OMP_NUM_THREADS=24 ./src/sand 0 h 128 F > res/h/128/F24.txt
OMP_NUM_THREADS=24 ./src/sand 0 h 128 P > res/h/128/P24.txt
OMP_NUM_THREADS=24 ./src/sand 0 h 128 t > res/h/128/t24.txt

OMP_NUM_THREADS=24 ./src/sand 0 c 512 F > res/c/512/F24.txt
OMP_NUM_THREADS=24 ./src/sand 0 c 512 P > res/c/512/P24.txt
OMP_NUM_THREADS=24 ./src/sand 0 c 512 t > res/c/512/t24.txt

OMP_NUM_THREADS=24 ./src/sand 0 h 512 F > res/h/512/F24.txt
OMP_NUM_THREADS=24 ./src/sand 0 h 512 P > res/h/512/P24.txt
OMP_NUM_THREADS=24 ./src/sand 0 h 512 t > res/h/512/t24.txt
