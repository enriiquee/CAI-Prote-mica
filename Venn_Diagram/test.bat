@echo off
color 30
echo ============================================================================
echo =                                                                          =
echo =                        VENN DIAGRAM PROGRAM                              =
echo =                                                                          =
echo ============================================================================
  echo.                     
echo.
echo Cargando... 
%~d0\ENRIQUE\R\R-3.4.0\bin\i386\Rscript.exe %~d0\ENRIQUE\Programas\Venn_Diagram\test_accesion.R
echo Finalizado. Presiona cualquier tecla para salir
pause>nul
exit

