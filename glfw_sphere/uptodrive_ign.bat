@echo off
setlocal
set SRCD=D:\work\lecture\src_ign\acg\glfw_sphere_ign
set DSTD=D:\drive\code\mac\glfw_sphere
set SRC1=%SRCD%\*.cpp
set SRC2=%SRCD%\*.h
set SRC3=%SRCD%\imgui
set SRC4=%SRCD%\inc
set SRC5=%SRCD%\lib
set SRC6=%SRCD%\bin
set SRC7=%SRCD%\Makefile
set DST1=%DSTD%\
set DST3=%DSTD%\imgui
set DST4=%DSTD%\inc
set DST5=%DSTD%\lib
set DST6=%DSTD%\bin
echo copy files to %DEST%
copy "%SRC1%" "%DST1%" /V /Y
copy "%SRC2%" "%DST1%" /V /Y
xcopy "%SRC3%" "%DST3%" /V /E /Y /I 
xcopy "%SRC4%" "%DST4%" /V /E /Y /I 
xcopy "%SRC5%" "%DST5%" /V /E /Y /I 
xcopy "%SRC6%" "%DST6%" /V /E /Y /I 
copy "%SRC7%" "%DST1%" /V /Y
endlocal
rem pause