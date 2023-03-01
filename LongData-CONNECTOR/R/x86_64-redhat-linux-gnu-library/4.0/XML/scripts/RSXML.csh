if(`test -n "-lxml2 -lz -llzma -lm -ldl"`) then

if(${?LD_LIBRARY_PATH}) then
    setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:-lxml2 -lz -llzma -lm -ldl
else
   setenv LD_LIBRARY_PATH -lxml2 -lz -llzma -lm -ldl
endif

endif
