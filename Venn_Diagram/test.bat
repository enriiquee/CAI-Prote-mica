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
"C:\Program Files\R\R-3.4.0\bin\i386\Rscript.exe" %cd%\test_accesion.R
echo Finalizado. Presiona cualquier tecla para salir
pause>nul
exit
C:\Program Files\R\R-3.4.0\bin\i386