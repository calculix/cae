@echo off

set DIR=%~dp0%
set CAE=%DIR%src\cae.exe

if exist %CAE% (
    echo Running binary.
    start "" %CAE% %*
) else (
    echo Binary does not exist. Running source code.
    set CAE=%DIR%src\cae.py
    start "" pythonw %CAE% %*
)
