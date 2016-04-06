@ECHO OFF &SETLOCAL
SET "Species=nsho"
set "Myvar=c:\nco\ncrcat.exe "
for /l %%a in (1979,1,2013) do call set "Myvar=%%Myvar%% WSI_DISPLAY_%%Species%%.%%a.nc"
call set "Myvar=%%Myvar%% WSI_DISPLAY_%%Species%%.1979_2013.nc"
ECHO %Myvar%