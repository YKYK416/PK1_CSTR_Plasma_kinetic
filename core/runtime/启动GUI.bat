@echo off
chcp 65001 >nul
rem ============================================================
rem ZDPlasKin Console GUI launcher (double-click)
rem %~dp0 = this bat file directory — portable with the folder
rem Python fallback: python_path.txt -> where python -> py -3
rem ============================================================
cd /d "%~dp0"
set "PY="

rem 1) interpreter recorded by last successful GUI run
if exist "%~dp0python_path.txt" (
  for /f "usebackq delims=" %%i in ("%~dp0python_path.txt") do (
    if not defined PY if exist "%%i" set "PY=%%i"
  )
)

rem 2) python on PATH
if not defined PY (
  for /f "delims=" %%i in ('where python 2^>nul') do (
    if not defined PY set "PY=%%i"
  )
)

rem 3) Windows py launcher
if not defined PY (
  where py >nul 2>nul && set "PY=py -3"
)

if not defined PY (
  echo [ERROR] Python interpreter not found.
  echo Tried: python_path.txt / where python / py -3 — all failed.
  echo Install Python 3 or run the GUI once on this machine.
  pause
  exit /b 1
)

echo Using interpreter: %PY%
if /i "%PY%"=="py -3" (
  py -3 "%~dp0zdp_gui.py"
) else (
  "%PY%" "%~dp0zdp_gui.py"
)
if errorlevel 1 pause
