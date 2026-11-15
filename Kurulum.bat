@echo off
chcp 65001 >nul
cd /d "%~dp0"
if not exist app.py (
  echo HATA: Dosyalar zip icinden calistiriliyor. Once zip'e sag tiklayip "Tumunu ayikla" de,
  echo sonra cikan klasordeki Kurulum.bat'i calistir.
  pause & exit /b
)
where py >nul 2>nul && (set PY=py) || (set PY=python)
%PY% --version >nul 2>nul
if errorlevel 1 (
  echo HATA: Python bulunamadi.
  echo python.org/downloads adresinden Python 3.11 veya 3.12 kur.
  echo Kurulumun ilk ekraninda "Add python.exe to PATH" kutusunu MUTLAKA isaretle.
  pause & exit /b
)
%PY% --version
%PY% -m pip install --upgrade pip
%PY% -m pip install -r requirements.txt
%PY% -m pip install pywebview
echo.
echo Kurulum bitti. Simdi Baslat.bat'a cift tikla.
pause
