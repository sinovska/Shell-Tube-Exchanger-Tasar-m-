# E-102 Kısmi Yoğuşturucu Tasarım Uygulaması

Ondokuz Mayıs Üniversitesi Kimya Mühendisliği Bölümü, **Kimya Mühendisliğinde Tasarım I** dersi kapsamında hazırlanmıştır.

İzopropil alkolden (IPA) aseton üreten bir tesiste, reaktörden çıkan gazı soğutma suyuyla soğutan ve yoğuşturan
**E-102 kabuk-tüp eşanjörünü** tasarlar. Gazdaki aseton, IPA ve su yoğuşur; hidrojen yoğuşmaz. Bu yüzden eşanjör bir
**kısmi yoğuşturucudur** uygulama bu etkileri hesaba katarak
**TEMA standardına uygun, teknik olarak en güvenilir** tasarımı otomatik olarak bulur.

## Ne yapar?

- **Termodinamik:** Gaz için Peng–Robinson, sıvı için NRTL. Isı yükü, çiy noktası, yoğuşan ve kaçan miktar.
- **Bölgesel ısıl hesap:** Eşanjör 15 bölgeye ayrılır; her bölge kendi koşullarıyla hesaplanır.
- **Isı transferi:**
  - Tüp tarafı: Gnielinski + Colebrook.
  - Kabuk tarafı: Bell–Delaware.
  - Yoğuşma: Nusselt filmi; hidrojen direnci için Colburn–Hougen.
- **Basınç düşüşleri:** Tüp tarafı ve iki fazlı kabuk tarafı, en kötü durum dahil.
- **Mekanik:**
  - TEMA tipi (ısıl genleşme kontrolü), et ve ayna kalınlıkları (TEMA Ek A), nozullar ve ρv² sınırları.
  - Tüpler aynaya tek tek yerleştirilerek sayılır.
- **Maliyet:** Yatırım, soğutma suyu, pompa elektriği ve yıllık toplam maliyet.
- **Otomatik tasarım:**
  - **Denenenler:** tüp çapı ve boyu, geçiş sayısı, kabuk çapı, baffle aralığı ve kesimi, tüp dizilimi (30°/90°/45°), tüp adımı, su sıcaklıkları.
  - **Sıralama:** (1) TEMA zorunlu kuralları → (2) kritik ölçütler (belirsizliklere dayanıklılık) → (3) toplam ölçüt sayısı → (4) maliyet.
- **Çıktılar:** TEMA spesifikasyon föyü (PDF), tüm hesaplar (Excel), grafikler, adım adım hesap föyü.
- **DWSIM doğrulaması (isteğe bağlı):** Termodinamik sonuçlar DWSIM'de bağımsız olarak hesaplanıp karşılaştırılır.

## Kurulum

Python 3.10 veya üstü gerekir.

**Windows (en kolay):**
1. Depoyu indirin (Code → Download ZIP) ve bir klasöre açın.
2. `Kurulum.bat` dosyasına çift tıklayın; gerekli paketler kurulur. Yalnızca ilk seferde yapılır.
3. `Baslat.bat` dosyasına çift tıklayın; uygulama masaüstü penceresinde açılır.

**Komut satırı (tüm sistemler):**
```bash
pip install -r requirements.txt
streamlit run app.py
```

## Kullanım

1. **Girdileri kontrol edin (sol menü):** proses akışı (yıllık besleme, bileşim, sıcaklıklar, basınçlar), soğutma suyu
   sıcaklığı ve izin verilen basınç düşüşleri. Varsayılan değerler projedeki tasarım koşullarıdır.
2. **"Benim için tasarla"ya basın.** Uygulama bütün seçenekleri dener ve en uygun tasarımı girdilere yazar (5–6 dakika).
3. **Sonuçları inceleyin:**
   - **Otomatik tasarım:** seçilen tasarım, "Ne seçildi, neden?" tablosu ve en iyi 5 aday.
   - **Özet, Mekanik, Maliyet ve diğer sekmeler:** ayrıntılı hesaplar ve TEMA kontrolleri.
   - **Hesap adımları:** her denklem, sayılar yerine konmuş haliyle.
4. **TEMA Rapor/Hesaplamalar sekmesinden** PDF ve Excel çıktılarını indirin.
5. İsterseniz bir girdiyi değiştirip **Hesapla** ile sonucun nasıl değiştiğini görün.

**Görünümü düzenleme:** Sol menüdeki "Görünümü düzenle" ile kod yazmadan yazıları, sekmeleri, görselleri ve renkleri
değiştirebilirsiniz. Ayarlar `gorunum.json` dosyasına, görseller `gorseller/` klasörüne kaydedilir.

## DWSIM doğrulaması (isteğe bağlı)

DWSIM **zorunlu değildir**; uygulama DWSIM olmadan tam olarak çalışır.

DWSIM ile karşılaştırma yapmak isterseniz:
- Bilgisayar **Windows** olmalı ve **DWSIM 8 veya üstü** kurulu olmalıdır (ücretsiz: https://dwsim.org).
- `Kurulum.bat` gerekli Python paketlerini (`pythonnet`, `pywin32`) de kurar.
- DWSIM kuruluysa uygulama bunu kendisi bulur. Sol menüde "Tasarımdan sonra DWSIM ile doğrula" kutusu çıkar;
  işaretliyken "Benim için tasarla" sonunda karşılaştırma kendiliğinden yapılır.
- **Karşılaştırılanlar:** ısı yükü, çiy noktası, çıkış gaz ve sıvı debileri, su tarafı enerji dengesi, akışkan özellikleri.
- **Bilgi amaçlı:** DWSIM'in eşanjör modülü yoğuşmayan gazlı kısmi yoğuşmayı modellemediği için geometri (U, basınç düşüşü)
  karşılaştırması yalnızca **DWSIM karşılaştırma** sekmesinde, bilgi amaçlı yapılır.
- Elle çalıştırmak için `DWSIM_Test.bat` kullanılabilir; sonuçlar `dwsim_sonuclar.csv` dosyasına yazılır.

## Dosyalar

| Dosya | İçerik |
|---|---|
| `app.py` | Arayüz: girdiler, sekmeler, otomatik tasarım, çıktılar |
| `esanjor_core.py` | Hesap motoru (arayüzden bağımsız; `import esanjor_core` ile ayrıca kullanılabilir) |
| `tema_foy.py` | TEMA spesifikasyon föyü (PDF ve Excel) |
| `grafikler.py` | Grafikler |
| `gorunum.py`, `gorunum.json` | Kod yazmadan görünüm düzenleme ve kaydedilmiş ayarlar |
| `gorseller/` | Uygulamaya eklenen görseller |
| `dwsim_dogrulama.py` | DWSIM ile otomatik karşılaştırma betiği |
| `masaustu.py` | Masaüstü penceresi |
| `Kurulum.bat`, `Baslat.bat`, `DWSIM_Test.bat` | Windows için kurulum, başlatma ve DWSIM testi |
| `ornek_yogusma_egrisi.csv` | Yoğuşma eğrisini tablo olarak girmek için şablon |
| `requirements.txt` | Gerekli Python paketleri |

> **Not:** İngilizce TEMA formu çıktısı için TEMA'nın boş formu (`tema_form.pdf`) gerekir. Telif hakkı TEMA'ya ait
> olduğu için depoda yoktur; elinizde varsa uygulama klasörüne koymanız yeterlidir. Türkçe föy ve Excel çıktısı
> bu dosya olmadan da çalışır.

## Varsayımlar ve sınırlar

- **Varsayımlar:**
  - Kararlı rejim; çevreye ısı kaybı yok.
  - Reaktörde yan ürün oluşmaz; gaz eşanjöre tamamen buhar olarak girer.
  - Her noktada faz dengesi vardır; hidrojen yoğuşmaz.
  - Kabuk basıncı doğrusal düşer; yoğuşmadan doğan küçük basınç kazancı ihmal edilir (güvenli tarafta).
- **Varsayılan tasarım değerleri:**
  - İzin verilen basınç düşüşü: kabuk 20 kPa, tüp 35 kPa. Güvenlik payı: %15–35.
  - Kirlenme direnci: 0,0002 m²K/W. Malzeme: karbon çeliği, 3 mm korozyon payı. Standart: TEMA 9. baskı, Sınıf B.
- **Sınırlar:**
  - Akış kaynaklı titreşim hesaplanmaz.
  - En belirsiz büyüklük, hidrojenli yoğuşmada kabuk tarafı ısı transfer katsayısıdır; kritik ölçütler bu belirsizliği sınar.
  - Maliyet bağıntıları ve fiyatlar ön tasarım içindir.
  - Sonuçlar eğitim amaçlı ön tasarım niteliğindedir; imalattan önce ayrıntılı mekanik hesap gerekir.

## Hazırlayan

Sinem N. ÖZYÜREK (21060414) 
Ondokuz Mayıs Üniversitesi, Kimya Mühendisliği Bölümü, 2026
