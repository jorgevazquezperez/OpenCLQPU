# OpenQCL
Proyecto que tiene como finalidad acabar inmplementando una QPU como parte del estándar de OpenCL.

## Estado actual
Ahora mismo hemos conseguido implementar el malloc en el _kernel_ empleando el modelo KMA (Spliet et al., [2014](https://doi.org/10.1145/2588768.2576781)), el archivo de prueba se encuentra en `samples/kernel_malloc`.

He tenido problemas con el printf dentro del _kernel_ los cuales ya he solucionado (se debía a la necesidad por parte de OpenCL de llamar a `clFinish()` cuando acaba de operar).

Estos dos _samples_ he intentado ejecutarlos por separado y comprobar que funcionan correctamente porque los voy a necesitar 100% operativos en el _kernel_ de `src/quantum`. El `malloc` porque será el sustituto natural del `malloc` de gcc y el `printf` porque lo utilizaré a lo largo del proyecto como forma de debuggear los _kernel_ de OpenCL junto con la variable de entorno `POCL_DEBUG`.

Añadimos la función `rand` al kernel usando el algoritmo del paper (Langdon, William, [2009](http://www0.cs.ucl.ac.uk/staff/ucacbbl/ftp/papers/langdon_2009_CIGPU.pdf)).

Hecho un _sample_ de la matrices en el kernel tanto de tamaño fijo 4x4 como de tamaño arbitrario. Estas últimas al final no las hemos empleado en el kernel principal, pero las tenemos implementadas por si acaso fuese necesario en un futuro.

Ya he ejecutado el código contenido en `src/quantum` y he conseguido crear un circuito con una puerta Hadamard y leer el resultado. Ya tenemos la prueba de concepto usando Qulacs.

## Próximos pasos
Hay que realizar un _driver_ de una supuesta QPU. Para ello por ahora el mejor ejemplo que tenemos de por dónde deberíamos de ir es este: [OpenCL Implementation using pocl for zedboard](https://github.com/abisheksethu/opencl-implementation)
