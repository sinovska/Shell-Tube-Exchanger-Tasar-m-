@echo off
chcp 65001 >nul
cd /d "%~dp0"
if not exist app.py (
  echo HATA: Once zip'i "Tumunu ayikla" ile bir klasore cikar.
  pause & exit /b
)
where py >nul 2>nul && (set PY=py) || (set PY=python)
%PY% masaustu.py
if errorlevel 1 pause
