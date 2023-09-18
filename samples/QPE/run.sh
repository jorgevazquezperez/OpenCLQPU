rm -r build/
EX_NAME=""quantum.out""
cmake . -B build -DMY_VARIABLE:STRING=${EX_NAME}
cd build
make

blue=$(tput setaf 3)
normal=$(tput sgr0)
printf "${blue}Execution of ${EX_NAME}${normal}\n\n"
./${EX_NAME} 2> errores.txt