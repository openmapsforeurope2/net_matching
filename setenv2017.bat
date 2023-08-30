@echo OFF

if not "%VISUALCVARS%" == "TRUE" (

	set currentDir=%CD%

	IF EXIST "C:\Program Files (x86)\Microsoft Visual Studio\2017\Enterprise\VC\Auxiliary\Build\vcvarsall.bat" (
		call "C:\Program Files (x86)\Microsoft Visual Studio\2017\Enterprise\VC\Auxiliary\Build\vcvarsall.bat" x86
	) else (
		call "C:\Program Files (x86)\Microsoft Visual Studio\2017\Professional\VC\Auxiliary\Build\vcvarsall.bat" x86
	)

	d:
	%currentDir:~0,2%
	cd %currentDir%
	
	set VISUALCVARS=TRUE
)
