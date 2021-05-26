# bi290-project
our project can **solve the problem of finding LCA in the tree-liked GO plot** (which absolutely isn't a tree struct) using the jumping method (or multiplication method)

we use the slim GO file *goslim_generic.obo* as a demonstration, you can download any GO file to analyze in the GO database :http://geneontology.org/

we should first preprocess the data using the instruction:

`cat goslim_generic.obo | grep -Po "(?<!_)id: GO:[0-9]{7}|is_a: GO:[0-9]{7}" | awk -f std.awk >slim_extract`

you can change the processed file name, and you can change the file name in the source file main.c, make sure the macro definition N is suitable for the data size

before using makefile to compile and link the source files ,you should first creat some directories to better create the dynamic shared lib:

`mkdir bin`

`mkdir shared`

`mkdir lib`

add the environment  variable of the include path using `export C_INCLUDE_PATH=./include:${C_INCLUDE_PATH}`

then you can call `make` to build the lib file libgolca.so (make sure now you are in the src directory)

using `make move` to move the .o file to shared/ , lib file to lib/ and binary file to bin/ 

before executing , make sure the dynamic lib path is added by using `export LD_LIBRARY_PATH=./lib${LD_LIBRARY_PATH}`

see more details in manual.

