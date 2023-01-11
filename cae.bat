@echo off

%~dp0%bin\python\python.exe .\src\cae.py
if %ERRORLEVEL% neq 0 goto catch
pause >nul
exit /b 0

:catch
echo Something went wrong.
pause >nul
exit /b 1