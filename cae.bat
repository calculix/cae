SET qwe=%1
IF "%qwe:~0,2%" == "-a" (
  del "*.res"
  del "*.prt"
  del "*.mdl"
  del "*.stt"
  del "*.sim"
  del "*.inp"
)
IF "%qwe:~0,3%" == "-aa" (
  del "centroids.txt"
  del "abaqusMacros.py"
)

src/cae.exe