"""
E-102 tasarımını DWSIM'de otomatik olarak kurar, çözer ve uygulamanın sonuçlarıyla karşılaştırır.

Ne yapar?
  1. Termodinamik test  : Cooler ile 120 → 40 °C soğutma (ısı yükü, çıkış gaz/sıvı debisi) ve çiy noktası.
  2. Eşanjör testi      : Heat Exchanger, "Shell and Tube Rating" modunda. Girdi dosyasında "tasarimlar"
                          listesi varsa en uygun ilk tasarımların hepsi ayrı ayrı denenir.
  3. U·A testi (ek)     : Heat Exchanger, uygulamanın U'su ve gerekli alanı verilerek. DWSIM tek bir LMTD
                          kullandığı için bölgesel hesabın etkisini gösterir (bilgi amaçlı).
  Sonuçları ekrana tablo olarak yazar, dwsim_sonuclar.csv dosyasına kaydeder ve
  DWSIM'de açılabilen E-102_dwsim.dwxmz dosyasını oluşturur.

Gereksinimler (yalnızca Windows):
  - DWSIM kurulu olmalı (8 veya üstü).
  - 64 bit Python 3.10–3.12 ve:  pip install pythonnet pywin32

Çalıştırma:
  python dwsim_dogrulama.py                      # varsayılan tasarım değerleriyle
  python dwsim_dogrulama.py dwsim_girdi.json     # uygulamadan indirilen güncel değerlerle
  DWSIM bulunamazsa:  set DWSIM_BIN=C:\\Program Files\\DWSIM   (kurulum klasörü)
"""
import csv
import json
import os
import sys

# ----------------------------------------------------------------------------
# 1) TASARIM DEĞERLERİ (uygulamanın son tasarımı). Uygulamanın DWSIM sekmesinden
#    "dwsim_girdi.json" indirilirse bu değerlerin yerine o dosya kullanılır.
# ----------------------------------------------------------------------------
GIRDI = {
    "bilesenler": {"Acetone": 10488.8, "Isopropanol": 1205.9, "Hydrogen": 364.1, "Water": 1691.3},  # kg/h
    "proses": {"T_giris_C": 120.0, "T_cikis_C": 40.0, "P_giris_bar": 2.8, "P_cikis_bar": 2.6},
    "su": {"T_giris_C": 20.0, "P_bar": 4.0, "debi_kgh": 262279.0},
    "geometri": {"Ds_mm": 650.0, "B_mm": 650.0, "baffle_kesim_yuzde": 25.0, "Do_mm": 25.4, "Di_mm": 22.1,
                 "L_m": 6.096, "N": 318, "gecis": 2, "Pt_mm": 31.75, "dizilim": "ucgen",
                 "kirlenme_kabuk": 0.0002, "kirlenme_tup": 0.0002, "k_duvar": 60.0, "puruzluluk_mm": 0.045},
    # NRTL: (bileşen i, bileşen j, A_ij, A_ji [cal/mol], alfa). A_ij = (g_ij - g_jj), DWSIM'in birimiyle.
    "nrtl": [["Water", "Acetone", 1110.9, 500.7, 0.3],
             ["Water", "Isopropanol", 1617.6, 11.1, 0.3],
             ["Acetone", "Isopropanol", 419.7, -20.7, 0.3]],
    "uygulama": {"Q_kW": 3045.0, "T_ciy_C": 85.1, "gaz_cikis_kgh": 2672.0, "sivi_cikis_kgh": 11078.0,
                 "su_cikis_C": 30.0, "U_kirli": 713.1, "A_m2": 154.69, "A_gerekli_m2": 115.65,
                 "dP_tup_kPa": 14.59, "dP_kabuk_kPa": 15.07},
}

# Uygulamanın kullandığı ad → DWSIM veri tabanındaki olası adlar
AD_ESLEME = {"Isopropanol": ["Isopropanol", "2-Propanol", "Isopropyl alcohol"],
             "Acetone": ["Acetone"], "Hydrogen": ["Hydrogen"], "Water": ["Water"]}


# ----------------------------------------------------------------------------
# 2) DWSIM'e BAĞLANMA
# ----------------------------------------------------------------------------
def dwsim_klasoru():
    """DWSIM kurulum klasörünü bulur; bulamazsa None döndürür."""
    adaylar = [os.environ.get("DWSIM_BIN", ""),
               os.path.expandvars(r"%LOCALAPPDATA%\DWSIM"), os.path.expandvars(r"%LOCALAPPDATA%\DWSIM8"),
               os.path.expandvars(r"%LOCALAPPDATA%\DWSIM9"),
               r"C:\Program Files\DWSIM", r"C:\Program Files\DWSIM8", r"C:\Program Files (x86)\DWSIM"]
    for k in adaylar:
        if k and os.path.exists(os.path.join(k, "DWSIM.Automation.dll")):
            return k
    return None


def net_yukle(DW):
    """DWSIM'in .NET kütüphanelerini yükler (yalnızca betik çalıştırılınca; uygulama içe aktarırken değil)."""
    global Array, Double, String, Dictionary, Directory, Automation3, Settings, ObjectType, NRTL_IPData
    global Cooler, HeatExchangerCalcMode
    import pythonnet
    try:
        pythonnet.load("netfx")          # DWSIM .NET Framework uygulamasıdır
    except Exception:
        pass
    import clr
    try:
        import pythoncom                 # bazı sürümlerde COM başlatılması gerekir
        pythoncom.CoInitialize()
    except ImportError:
        pass
    for alt in ("", "unitops", "ppacks", "extenders"):
        sys.path.append(os.path.join(DW, alt))
    for dll in ("CapeOpen", "DWSIM.Automation", "DWSIM.Interfaces", "DWSIM.GlobalSettings", "DWSIM.SharedClasses",
                "DWSIM.Thermodynamics", "DWSIM.UnitOperations", "DWSIM.Inspector", "System.Buffers"):
        yol = os.path.join(DW, dll + ".dll")
        if os.path.exists(yol):
            clr.AddReference(yol)
    from System import Array, Double, String
    from System.Collections.Generic import Dictionary
    from System.IO import Directory
    from DWSIM.Automation import Automation3
    from DWSIM.GlobalSettings import Settings
    from DWSIM.Interfaces.Enums.GraphicObjects import ObjectType
    from DWSIM.Thermodynamics.PropertyPackages.Auxiliary import NRTL_IPData
    from DWSIM.UnitOperations.UnitOperations import Cooler, HeatExchangerCalcMode
    Directory.SetCurrentDirectory(DW)


def somut(o):
    """pythonnet arayüz nesnesini gerçek sınıfa çevirir."""
    return getattr(o, "__implementation__", o)


# ----------------------------------------------------------------------------
# 3) AKIŞ ŞEMASINI KURMA
# ----------------------------------------------------------------------------
def akis_semasi_kur(g):
    interf = Automation3()
    sim = interf.CreateFlowsheet()

    # Bileşenler (DWSIM'deki adlarıyla)
    adlar = {}
    for ad in g["bilesenler"]:
        for aday in AD_ESLEME.get(ad, [ad]):
            if sim.AvailableCompounds.ContainsKey(aday):
                sim.AddCompound(aday)
                adlar[ad] = aday
                break
        else:
            sys.exit(f"'{ad}' DWSIM veri tabanında bulunamadı.")
    print("Bileşenler:", ", ".join(adlar.values()))

    # NRTL paketi ve uygulamanın parametreleri
    pp = somut(sim.CreateAndAddPropertyPackage("NRTL"))
    ipd = pp.m_uni.InteractionParameters
    for i, j, aij, aji, alfa in g["nrtl"]:
        i, j = adlar.get(i, i), adlar.get(j, j)
        if ipd.ContainsKey(j) and ipd[j].ContainsKey(i):       # ters sıradaki kayıt varsa sil
            ipd[j].Remove(i)
        if not ipd.ContainsKey(i):
            ipd.Add(i, Dictionary[String, NRTL_IPData]())          # boş iç sözlük
        p = NRTL_IPData()
        p.ID1, p.ID2, p.A12, p.A21, p.alpha12 = i, j, float(aij), float(aji), float(alfa)
        if ipd[i].ContainsKey(j):
            ipd[i].Remove(j)
        ipd[i].Add(j, p)
    print("NRTL parametreleri girildi.")
    return interf, sim, adlar


def akim(sim, ad, x, y, T_C, P_bar, debiler_kgh, adlar):
    ms = somut(sim.AddObject(ObjectType.MaterialStream, x, y, ad).GetAsObject())
    sira = [c.Name for c in ms.Phases[0].Compounds.Values]
    toplam = sum(debiler_kgh.values())
    w = [0.0] * len(sira)
    for k, v in debiler_kgh.items():
        w[sira.index(adlar.get(k, k))] = v / toplam
    ms.SetOverallMassComposition(Array[Double](w))
    ms.SetTemperature(T_C + 273.15)
    ms.SetPressure(P_bar * 1e5)
    ms.SetMassFlow(toplam / 3600.0)
    return ms


def bos_akim(sim, ad, x, y):
    return somut(sim.AddObject(ObjectType.MaterialStream, x, y, ad).GetAsObject())


def baglan(sim, a, b, ia=-1, ib=-1):
    sim.ConnectObjects(a.GraphicObject, b.GraphicObject, ia, ib)


def faz_debisi(ms, faz):
    try:
        return ms.GetPhase(faz).Properties.massflow * 3600.0
    except Exception:
        return float("nan")


def bilesen_debisi(ms, faz, bilesen):
    try:
        return ms.GetPhase(faz).Compounds[bilesen].MassFlow * 3600.0
    except Exception:
        return float("nan")


def geometri_gir(hx, geo):
    hx.CalculationMode = HeatExchangerCalcMode.ShellandTube_Rating
    st = hx.STProperties
    st.Shell_NumberOfShellsInSeries = 1
    st.Shell_NumberOfPasses = 1
    st.Shell_Di = geo["Ds_mm"]
    st.Shell_BaffleSpacing = geo["B_mm"]
    st.Shell_BaffleCut = geo["baffle_kesim_yuzde"]
    st.Shell_BaffleType = 0                    # tek segmental
    st.Shell_BaffleOrientation = 1             # dikey kesim
    st.Shell_Fouling = geo["kirlenme_kabuk"]
    st.Shell_Fluid = 1                         # kabukta sıcak akışkan
    st.Tube_Fluid = 0                          # tüpte soğuk akışkan
    st.Tube_De = geo["Do_mm"]
    st.Tube_Di = geo["Di_mm"]
    st.Tube_Length = geo["L_m"]
    st.Tube_NumberPerShell = int(geo["N"])
    st.Tube_PassesPerShell = int(geo["gecis"])
    st.Tube_Pitch = geo["Pt_mm"]
    st.Tube_Layout = 0 if geo["dizilim"].startswith("u") else 2
    st.Tube_Fouling = geo["kirlenme_tup"]
    st.Tube_ThermalConductivity = geo["k_duvar"]
    st.Tube_Roughness = geo["puruzluluk_mm"]


# ----------------------------------------------------------------------------
# 4) TESTLER
# ----------------------------------------------------------------------------
def uyum_yuzde(a, d):
    """Benzerlik: 100 − bağıl fark (%), 0'ın altına inmez."""
    try:
        if a is None or d is None or a != a or d != d or a == 0:
            return None
        return max(0.0, 100.0 - abs(d - a) / abs(a) * 100.0)
    except Exception:
        return None


def main():
    for akis in (sys.stdout, sys.stderr):          # Windows'ta çıktı yönlendirilince cp1254 "Δ" yazamaz
        try:
            akis.reconfigure(encoding="utf-8", errors="replace")
        except Exception:
            pass
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    cikti = next((a.split("=", 1)[1] for a in sys.argv[1:] if a.startswith("--cikti=")), None)
    DW = dwsim_klasoru()
    if DW is None:
        sys.exit("DWSIM bulunamadı. Kurulum klasörünü DWSIM_BIN ortam değişkeniyle belirtin, ör.:\n"
                 "  set DWSIM_BIN=C:\\Program Files\\DWSIM")
    print(f"DWSIM klasörü: {DW}")
    net_yukle(DW)
    g = GIRDI
    dosya = args[0] if args else os.path.join(os.path.dirname(os.path.abspath(__file__)), "dwsim_girdi.json")
    if os.path.exists(dosya):
        with open(dosya, encoding="utf-8") as f:
            g = json.load(f)
        print(f"Girdiler '{os.path.basename(dosya)}' dosyasından okundu.")
    pr, su, geo, app = g["proses"], g["su"], g["geometri"], g["uygulama"]
    interf, sim, adlar = akis_semasi_kur(g)
    su_debi = {"Water": su["debi_kgh"]}

    # --- Test 1: soğutucu (termodinamik)
    p_in = akim(sim, "PROSES-1", 50, 50, pr["T_giris_C"], pr["P_giris_bar"], g["bilesenler"], adlar)
    p_out = bos_akim(sim, "PROSES-1-CIKIS", 250, 50)
    e1 = somut(sim.AddObject(ObjectType.EnergyStream, 150, 120, "Q-SOGUTUCU").GetAsObject())
    c1 = somut(sim.AddObject(ObjectType.Cooler, 150, 50, "SOGUTUCU").GetAsObject())
    baglan(sim, p_in, c1); baglan(sim, c1, p_out); baglan(sim, c1, e1)
    c1.CalcMode = Cooler.CalculationMode.OutletTemperature
    c1.OutletTemperature = pr["T_cikis_C"] + 273.15
    c1.DeltaP = (pr["P_giris_bar"] - pr["P_cikis_bar"]) * 1e5

    # Çiy noktası: aynı gazı basınç düşüşü olmadan buhar oranı 1'e kadar soğut
    d_in = akim(sim, "PROSES-2", 50, 200, pr["T_giris_C"], pr["P_giris_bar"], g["bilesenler"], adlar)
    d_out = bos_akim(sim, "CIY-NOKTASI", 250, 200)
    e2 = somut(sim.AddObject(ObjectType.EnergyStream, 150, 270, "Q-CIY").GetAsObject())
    c2 = somut(sim.AddObject(ObjectType.Cooler, 150, 200, "CIY-HESABI").GetAsObject())
    baglan(sim, d_in, c2); baglan(sim, c2, d_out); baglan(sim, c2, e2)
    c2.CalcMode = Cooler.CalculationMode.OutletVaporFraction
    c2.OutletVaporFraction = 1.0
    c2.DeltaP = 0.0

    # --- Test 2: kabuk-tüp eşanjör, rating modu (bir veya birden fazla tasarım)
    tasarimlar = g.get("tasarimlar") or [dict(ad="Seçilen", su=su, geometri=geo, uygulama=app)]
    rating = []
    for n, t in enumerate(tasarimlar):
        y = 350 + 250 * n
        hi = akim(sim, f"PROSES-T{n + 1}", 50, y, pr["T_giris_C"], pr["P_giris_bar"], g["bilesenler"], adlar)
        si = akim(sim, f"SU-T{n + 1}", 50, y + 100, t["su"]["T_giris_C"], su["P_bar"], {"Water": t["su"]["debi_kgh"]}, adlar)
        ho, so = bos_akim(sim, f"PROSES-T{n + 1}-CIKIS", 250, y), bos_akim(sim, f"SU-T{n + 1}-CIKIS", 250, y + 100)
        hxn = somut(sim.AddObject(ObjectType.HeatExchanger, 150, y + 50, f"E-102-T{n + 1}").GetAsObject())
        baglan(sim, hi, hxn, -1, 0); baglan(sim, si, hxn, -1, 1)
        baglan(sim, hxn, ho, 0, -1); baglan(sim, hxn, so, 1, -1)
        geometri_gir(hxn, t["geometri"])
        rating.append((t, hxn, ho, so))

    # --- Test 3: aynı eşanjör, U ve A verilerek
    y3 = 350 + 250 * len(tasarimlar)
    h_in2 = akim(sim, "PROSES-4", 50, y3, pr["T_giris_C"], pr["P_giris_bar"], g["bilesenler"], adlar)
    s_in2 = akim(sim, "SU-GIRIS-2", 50, y3 + 100, su["T_giris_C"], su["P_bar"], su_debi, adlar)
    h_out2, s_out2 = bos_akim(sim, "PROSES-4-CIKIS", 250, y3), bos_akim(sim, "SU-CIKIS-2", 250, y3 + 100)
    hx2 = somut(sim.AddObject(ObjectType.HeatExchanger, 150, y3 + 50, "E-102-UA").GetAsObject())
    baglan(sim, h_in2, hx2, -1, 0); baglan(sim, s_in2, hx2, -1, 1)
    baglan(sim, hx2, h_out2, 0, -1); baglan(sim, hx2, s_out2, 1, -1)
    hx2.CalculationMode = HeatExchangerCalcMode.CalcBothTemp_UA
    hx2.OverallCoefficient = app["U_kirli"]
    hx2.Area = app.get("A_gerekli_m2", app["A_m2"])
    hx2.HotSidePressureDrop = (pr["P_giris_bar"] - pr["P_cikis_bar"]) * 1e5
    hx2.ColdSidePressureDrop = app["dP_tup_kPa"] * 1000

    # --- Çöz
    Settings.SolverMode = 0
    print("DWSIM çözülüyor...")
    hatalar = interf.CalculateFlowsheet2(sim)
    if hatalar is not None and len(list(hatalar)) > 0:
        print("DWSIM uyarı/hataları:")
        for h in hatalar:
            print("  -", h.Message if hasattr(h, "Message") else h)

    kayit = os.path.join(os.path.dirname(os.path.abspath(__file__)), "E-102_dwsim.dwxmz")
    try:
        interf.SaveFlowsheet(sim, kayit, True)
        print(f"DWSIM dosyası kaydedildi: {kayit}")
    except Exception as e:
        print("DWSIM dosyası kaydedilemedi:", e)

    # --- Sonuçlar
    h2 = adlar.get("Hydrogen", "Hydrogen")
    satirlar = [
        # (test, büyüklük, birim, uygulama, DWSIM, kabul sınırı, tür: 'yuzde' veya 'mutlak')
        ("Termodinamik", "Isı yükü", "kW", app["Q_kW"], abs(c1.DeltaQ or 0.0), 3, "yuzde"),
        ("Termodinamik", "Çiy noktası", "°C", app["T_ciy_C"], d_out.GetTemperature() - 273.15, 2, "mutlak"),
        ("Termodinamik", "Çıkış gaz debisi", "kg/h", app["gaz_cikis_kgh"], faz_debisi(p_out, "Vapor"), 5, "yuzde"),
        ("Termodinamik", "Çıkış sıvı debisi", "kg/h", app["sivi_cikis_kgh"], faz_debisi(p_out, "OverallLiquid"), 5, "yuzde"),
        ("Termodinamik", "Sıvıdaki hidrojen", "kg/h", 0.0, bilesen_debisi(p_out, "OverallLiquid", h2), None, ""),
        ("Termodinamik", "Gazdaki aseton", "kg/h", app.get("aseton_gazda_kgh", float("nan")),
         bilesen_debisi(p_out, "Vapor", adlar.get("Acetone", "Acetone")), None, ""),
        ("U·A verilerek", "Isı yükü", "kW", app["Q_kW"], hx2.Q or float("nan"), None, ""),
        ("U·A verilerek", "Proses çıkış sıcaklığı", "°C", g["proses"]["T_cikis_C"], h_out2.GetTemperature() - 273.15, None, ""),
        ("U·A verilerek", "Su çıkış sıcaklığı", "°C", app["su_cikis_C"], s_out2.GetTemperature() - 273.15, None, ""),
    ]
    for t, hxn, ho, so in rating:
        ta, ua = str(t.get("ad", "Eşanjör"))[:17], t["uygulama"]
        satirlar += [
            (ta, "Isı yükü", "kW", ua["Q_kW"], hxn.Q or float("nan"), None, ""),
            (ta, "Proses çıkış sıcaklığı", "°C", pr["T_cikis_C"], ho.GetTemperature() - 273.15, None, ""),
            (ta, "Su çıkış sıcaklığı", "°C", ua["su_cikis_C"], so.GetTemperature() - 273.15, None, ""),
            (ta, "U (kirli)", "W/m²K", ua.get("U_kirli", float("nan")), hxn.OverallCoefficient or float("nan"), None, ""),
            (ta, "Tüp tarafı ΔP", "kPa", ua["dP_tup_kPa"], hxn.ColdSidePressureDrop / 1000, None, ""),
            (ta, "Kabuk ΔP", "kPa", ua["dP_kabuk_kPa"], hxn.HotSidePressureDrop / 1000, None, ""),
            (ta, "Kabuk ΔP (sınır)", "kPa", g.get("kabuk_dP_siniri_kPa", 20.0), hxn.HotSidePressureDrop / 1000,
             0, "sinir"),
        ]
    tablo, cizgiler = [], []
    for test, ad, birim, a, d, sinir, tur in satirlar:
        fark = d - a
        if tur == "sinir":
            fark_s, ok = f"{d / a * 100:.0f}% sınır", d <= a
        elif tur == "yuzde" and a:
            fark_s, ok = f"%{100 * fark / a:+.1f}", abs(100 * fark / a) <= sinir
        elif tur == "mutlak":
            fark_s, ok = f"{fark:+.1f}", abs(fark) <= sinir
        else:
            fark_s, ok = f"{fark:+.1f}", None
        if tur == "sinir":
            sonuc = "SINIR İÇİNDE" if ok else "SINIRI AŞIYOR"
        else:
            sonuc = "bilgi" if ok is None else ("UYGUN" if ok else "FARK VAR")
        cizgiler.append(f"{test:<18}{ad:<24}{birim:<7}{a:>10.1f}{d:>11.1f}{fark_s:>10}  {sonuc}")
        uy = None if tur == "sinir" or ad in ("Sıvıdaki hidrojen",) else uyum_yuzde(a, d)
        tablo.append(dict(test=test, buyukluk=ad, birim=birim, uygulama=round(a, 2), dwsim=round(d, 2),
                          fark=fark_s, sonuc=sonuc, uyum=None if uy is None else round(uy, 1)))

    if cikti:                                      # önce sonuç dosyası: ekrana yazma hatası sonucu kaybettirmesin
        with open(cikti, "w", encoding="utf-8") as f:
            json.dump(tablo, f, ensure_ascii=False, indent=1)
    try:
        print()
        print(f"{'Test':<18}{'Büyüklük':<24}{'Birim':<7}{'Uygulama':>10}{'DWSIM':>11}{'Fark':>10}  Sonuç")
        print("-" * 92)
        print("\n".join(cizgiler))
    except Exception as e:
        print("Tablo ekrana yazılamadı:", e)
    csv_yol = os.path.join(os.path.dirname(os.path.abspath(__file__)), "dwsim_sonuclar.csv")
    with open(csv_yol, "w", newline="", encoding="utf-8-sig") as f:
        yaz = csv.DictWriter(f, fieldnames=list(tablo[0]), delimiter=";")
        yaz.writeheader()
        yaz.writerows(tablo)
    print(f"\nSonuçlar kaydedildi: {csv_yol}")
    print("Not: DWSIM kabuk tarafını Tinker yöntemiyle, yoğuşma filmini ve hidrojenin difüzyon direncini ayrıca\n"
          "     modellemeden hesaplar; U ve basınç düşüşlerinde uygulamayla birebir uyum beklenmez.\n"
          "     U·A testinde DWSIM uç sıcaklıklardan tek bir LMTD kullanır; yoğuşma eğrisi doğrusal olmadığından\n"
          "     ısı yükünü fazla bulması beklenir (bölgesel hesabın gerekçesi).")


if __name__ == "__main__":
    main()
