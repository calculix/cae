@echo off

python %~dp0%src\cae.py
if %ERRORLEVEL% neq 0 goto catch
pause >nul
exit /b 0

:catch
echo Command 'python' is needed to run CAE.
echo Please, install Python v3.
pause >nul
exit /b 1