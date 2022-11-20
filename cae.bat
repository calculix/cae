@echo off
python -V | findstr "3."
if %ERRORLEVEL%==0 python %~dp0%src\cae.py
else echo "Python 3 needed."