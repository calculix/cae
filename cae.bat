@echo off
IF "%~1" == "-inp" (
    start "" "%~dp0\src\cae.exe" -inp %2
) ELSE (
    start "" "%~dp0\src\cae.exe" %*
)
