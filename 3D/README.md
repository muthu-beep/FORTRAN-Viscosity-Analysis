This Folder contains all the Files required to make the 3D simulations

## Commands
```console
adi@adi~$ pgfortran -c ini.f mmesh.f
adi@adi~$ pgfortran ini.o mmesh.o
adi@adi~$ ./a.out
adi@adi~$ ./lnk
adi@adi~$ ./nameofexecutable > results (it overwrites the file results)
adi@adi~$ ./nameofexecutable >> results (it adds to the file results) 
```
