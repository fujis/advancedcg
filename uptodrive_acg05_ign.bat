@echo off
setlocal
set SRCD=D:\work\lecture\src_ign\acg\acg05_mls_deformation
set DSTD=D:\drive\code\mac\Advanced05

set SRC1=%SRCD%\Advanced05\*.cpp
set SRC2=%SRCD%\Advanced05\*.h
set SRC3=%SRCD%\Advanced05\imgui
set SRC4=%SRCD%\include
set SRC6=%SRCD%\bin
set SRC7=%SRCD%\Advanced05\Makefile

set DST1=%DSTD%\Advanced05
set DST3=%DSTD%\Advanced05\imgui
set DST4=%DSTD%\include
set DST6=%DSTD%\bin

echo copy files to %DEST%
copy "%SRC1%" "%DST1%" /V /Y
copy "%SRC2%" "%DST1%" /V /Y
xcopy "%SRC3%" "%DST3%" /V /E /Y /I 
xcopy "%SRC4%" "%DST4%" /V /E /Y /I 
xcopy "%SRC6%" "%DST6%" /V /E /Y /I 
copy "%SRC7%" "%DST1%" /V /Y
endlocal
rem pause