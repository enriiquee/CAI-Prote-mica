@echo off
color 30
echo ============================================================================
echo =                                                                          =
echo =                        FILTER PROTEIN PILOT                              =
echo =                                                                          =
echo ============================================================================
  echo.                     
echo.
echo Cargando... 
"C:\Program Files\R\R-3.4.0\bin\i386\Rscript.exe" %cd%\Filter_PP.R
echo Finalizado. Presiona cualquier tecla para salir
exit
