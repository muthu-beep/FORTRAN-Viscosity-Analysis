This Folder contains all the Files required to make 2D Simulations

## Commands
```console
adi@adi~$ pgfortran -fast binewcg.f
adi@adi~$ ./a.out
```
To draw streamlines in the channel:
```console
adi@adi~$ pgfortran -c main.f
adi@adi~$ pgfortran -c velocity.f
adi@adi~$ pgfortran main.o velocity.o
adi@adi~$ pgfortran main.o velocity.o -o exe
adi@adi~$ ./exe 
```
