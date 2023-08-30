@echo OFF

call setenv2017.bat
call %SOCLE_EXTERN_ROOT%/win32-msvc2017/library_config.bat

IF "%BOOST_ROOT%"=="" set BOOST_ROOT=%SOCLE_EXTERN_ROOT%\win32-msvc2017\boost_1_61_0

SET MY_DIR=%CD%

rmdir /Q /S solution-msvc2017
mkdir solution-msvc2017

cd solution-msvc2017  

cmake -G"Visual Studio 15 2017" -DCMAKE_INSTALL_PREFIX="%MY_DIR%" ..\

cd ..