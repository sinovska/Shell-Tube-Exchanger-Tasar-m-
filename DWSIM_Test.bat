@echo off
chcp 65001 >nul
cd /d "%~dp0"
where py >nul 2>nul && (set PY=py) || (set PY=python)
echo Gerekli paketler kontrol ediliyor...
%PY% -m pip install --quiet pythonnet pywin32
echo.
echo DWSIM testi basliyor, 1-2 dakika surebilir...
echo.
%PY% dwsim_dogrulama.py
echo.
echo Bitti. Sonuclar "dwsim_sonuclar.csv" dosyasinda.
pause
