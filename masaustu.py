"""Eşanjör uygulamasını tarayıcı sekmesi yerine ayrı bir pencerede açar.

Sıra: 1) pywebview penceresi  2) Edge/Chrome "uygulama modu"  3) varsayılan tarayıcı.
Hata olursa hata.txt dosyasına yazılır.
"""
import os
import shutil
import socket
import subprocess
import sys
import tempfile
import time
import traceback
import webbrowser

KLASOR = os.path.dirname(os.path.abspath(__file__))
BASLIK = "E-102 Eşanjör Tasarımı"


def bos_port():
    with socket.socket() as s:
        s.bind(("127.0.0.1", 0))
        return s.getsockname()[1]


def sunucu_hazir(port, sunucu, sure=90):
    bitis = time.time() + sure
    while time.time() < bitis:
        if sunucu.poll() is not None:          # streamlit hemen kapandıysa
            return False
        try:
            with socket.create_connection(("127.0.0.1", port), timeout=1):
                return True
        except OSError:
            time.sleep(0.3)
    return False


def tarayici_bul():
    adaylar = [
        shutil.which("msedge"), shutil.which("chrome"),
        r"C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe",
        r"C:\Program Files\Microsoft\Edge\Application\msedge.exe",
        r"C:\Program Files\Google\Chrome\Application\chrome.exe",
        os.path.expandvars(r"%LOCALAPPDATA%\Google\Chrome\Application\chrome.exe"),
    ]
    return next((a for a in adaylar if a and os.path.exists(a)), None)


def pencere_ac(url):
    try:
        import webview
        webview.settings["ALLOW_DOWNLOADS"] = True       # indirme butonları çalışsın
        webview.create_window(BASLIK, url, width=1400, height=900, min_size=(1000, 700),
                              text_select=True)          # yazılar seçilip kopyalanabilsin
        webview.start()
        return
    except Exception as e:
        print("pywebview açılamadı, tarayıcı uygulama moduna geçiliyor:", e)

    tarayici = tarayici_bul()
    if tarayici:
        profil = tempfile.mkdtemp(prefix="esanjor_")   # ayrı profil: pencere kapanınca süreç biter
        p = subprocess.Popen([tarayici, f"--app={url}", f"--user-data-dir={profil}",
                              "--window-size=1400,900", "--no-first-run"])
        p.wait()
        return

    webbrowser.open(url)
    input("Uygulama tarayıcıda açıldı. Kapatmak için bu pencerede Enter'a bas...")


def main():
    port = bos_port()
    print("Hesap sunucusu başlatılıyor, lütfen bekle...")
    sunucu = subprocess.Popen(
        [sys.executable, "-m", "streamlit", "run", os.path.join(KLASOR, "app.py"),
         "--server.port", str(port), "--server.headless", "true"],
        cwd=KLASOR)
    try:
        if not sunucu_hazir(port, sunucu):
            raise RuntimeError("Streamlit başlamadı. Kurulum.bat'ı çalıştırdığından emin ol "
                               "(yukarıdaki satırlarda hata yazıyor olabilir).")
        pencere_ac(f"http://127.0.0.1:{port}")
    finally:
        sunucu.terminate()


if __name__ == "__main__":
    try:
        main()
    except Exception:
        hata = traceback.format_exc()
        print("\nHATA:\n" + hata)
        with open(os.path.join(KLASOR, "hata.txt"), "w", encoding="utf-8") as f:
            f.write(hata)
        input("Bu pencereyi kapatmak için Enter'a bas...")
