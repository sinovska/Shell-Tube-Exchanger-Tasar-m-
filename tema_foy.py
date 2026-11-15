"""
TEMA ısı değiştirici spesifikasyon föyü (TEMA Specification Sheet, 1–64 satır düzeni).

Tek veri kaynağı: collect() uygulamanın hesap sonuçlarından föy değerlerini toplar.
Çıktılar:
    pdf_tr()   Türkçe föy (A4, 2 sayfa)
    to_excel() 1. sayfa TEMA föyü, 2. sayfa tüm hesaplanan veriler, sonra ayrıntı tabloları
    to_html()  uygulamada önizleme
"""
import io
import math
import os

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))

STEEL_RHO = 7850.0                                    # kg/m³
PLATES = [5, 6, 8, 10, 12, 14, 16, 18, 20, 25, 30]    # standart levha kalınlıkları, mm
# TEMA RCB-3.13, karbon çeliği kabuk asgari kalınlığı (sınıf B/C), iç çap üst sınırı mm -> mm
TEMA_MIN_SHELL = [(305, 6.4), (737, 7.9), (991, 9.5), (1524, 11.1), (2032, 12.7), (99999, 12.7)]
BWG = {2.77: 12, 2.41: 13, 2.11: 14, 1.83: 15, 1.65: 16, 1.47: 17, 1.24: 18, 0.89: 20}
LAYOUT_DEG = {"triangular": "30°", "square": "90°", "rotated square": "45°"}
MATERIALS = {"Karbon çeliği": ("SA-179", "SA-516 Gr.70"), "Paslanmaz 304": ("SA-213 TP304", "SA-240 304"),
             "Paslanmaz 316": ("SA-213 TP316", "SA-240 316"), "Paslanmaz 321": ("SA-213 TP321", "SA-240 321")}
TEMA_CODE = {"Sabit ayna (BEM)": "BEM", "U-tüp (BEU)": "BEU", "Yüzer kafa, split-ring (AES)": "AES",
             "Yüzer kafa, pull-through (AET)": "AET"}

# Föyde kullanıcının değiştirebileceği seçimler (varsayılanlar karbon çeliği, kimya prosesi için tipik)
DEFAULT_CHOICES = dict(
    customer="", address="", location="", job_no="", ref_no="", proposal_no="",
    shell_fluid="Reaktör çıkış buharı", tube_fluid="Soğutma suyu",
    code="ASME VIII Div.1", tema_class="B", code_stamp=True, mawp=True, inspection=True,
    radiography="Nokta (E = 0,85)", special_service="Yok (yanıcı: H₂, aseton)", pwht=False,
    flange_class="150", facing="RF", vent_shell='1"', drain_shell='1"', vent_tube='3/4"', drain_tube='3/4"',
    CA_tube=3.0, baffle_cut=25, cut_dir="Dikey", tube_joint="Genişletme + oluk", supports="Semer (saddle)", n_supports=2,
    mat_tube="SA-179", mat_shell="SA-516 Gr.70", mat_tubesheet="SA-516 Gr.70", mat_channel="SA-516 Gr.70",
    mat_baffle="SA-285 Gr.C", mat_tierod="SA-36", mat_flange="SA-105", mat_neck="SA-106 Gr.B",
    gasket="Spiral, grafit", bolting="SA-193 B7 / SA-194 2H",
    insulation="Yalnız personel koruması (başkalarınca)", paint="Proje boya şartnamesine göre", mdmt="",
)

# =====================================================================
# Yardımcı hesaplar
# =====================================================================
def latent_heat(curve):
    """Yoğuşan kütleye göre ağırlıklı ortalama gizli ısı [kJ/kg]."""
    num = den = 0.0
    for i in range(len(curve) - 1):
        a, b = curve.iloc[i], curve.iloc[i + 1]
        dmL = (b.m_L - a.m_L) / 3600
        if dmL <= 1e-9:
            continue
        dT = a["T"] - b["T"]
        mV, mL = 0.5 * (a.m_V + b.m_V) / 3600, 0.5 * (a.m_L + b.m_L) / 3600
        cpV, cpL = 0.5 * (a.cp_V + b.cp_V), 0.5 * (a.cp_L + b.cp_L)
        num += ((b.Q_cum - a.Q_cum) - (mV * cpV + mL * cpL) * dT) / dmL * dmL
        den += dmL
    return num / den if den > 0 else np.nan


def weights(geo, t_shell_mm):
    """Yaklaşık ağırlıklar [kg]: tüpler + kabuk; ayna, kanal, nozul ve baffle için %30 pay."""
    Do, Di, L, Ds = geo["Do"] / 1000, (geo["Do"] - 2 * geo["t"]) / 1000, geo["L"] / 1000, geo["Ds"] / 1000
    tubes = math.pi / 4 * (Do**2 - Di**2) * L * geo["N"] * STEEL_RHO
    shell = math.pi * (Ds + t_shell_mm / 1000) * t_shell_mm / 1000 * L * STEEL_RHO
    empty = 1.30 * (tubes + shell)
    water = (math.pi / 4 * Ds**2 * L - math.pi / 4 * Do**2 * L * geo["N"] + math.pi / 4 * Di**2 * L * geo["N"]) * 1000
    return dict(empty=empty, full=empty + water, bundle=1.15 * tubes)


def shell_thickness(t_req, Ds):
    """Hesaplanan (KP dahil) kalınlık ile TEMA asgari kalınlığının büyüğü, standart levhaya yuvarlanır."""
    t_min = next(t for d, t in TEMA_MIN_SHELL if Ds <= d)
    return next((p for p in PLATES if p >= max(t_req, t_min)), PLATES[-1]), t_min


def _bundle_rhov2(inl, out, sg):
    """TEMA RCB-4.624: darbe plakası yokken demet giriş/çıkış alanı = baffle bölmesindeki tüpler arası çapraz akış
    alanı (Kern A_s). Çıkışta iki faz homojen yoğunluk. Döndürür: (giriş, çıkış) ρv², kg/(m·s²)."""
    G = (inl.m_V + inl.m_L) / 3600 / sg["As"]
    x = out.m_V / (out.m_V + out.m_L)
    rho_out = 1 / (x / out.rho_V + (1 - x) / out.rho_L) if out.m_L > 0 else out.rho_V
    return G**2 / inl.rho_V, G**2 / rho_out


def zonal_U(r):
    """Bölgesel sonuçlardan eşdeğer MTD ve U değerleri (föy ve özet aynı tanımı kullanır)."""
    z = r["zones"]
    UA = (z.Ud * z.A_m2).sum()
    mtd = r["Q_kW"] * 1000 / UA                      # F dahil
    return dict(mtd=mtd, U_req=r["Q_kW"] * 1000 / (r["A"] * mtd), U_dirty=UA / r["A_req"],
                U_clean=(z.Uc * z.A_m2).sum() / r["A_req"])


# =====================================================================
# 1) Föy verilerini topla
# =====================================================================
def collect(R, inp, ch=None):
    """R: uygulamanın hesap sonucu. inp: arayüz girdileri (basınçlar, sınırlar, proje bilgileri...).
    ch: DEFAULT_CHOICES biçiminde kullanıcı seçimleri. Döndürür: föyün tüm değerleri (sayılar ham)."""
    ch = {**DEFAULT_CHOICES, **(ch or {})}
    r, curve, m = R["r"], R["curve"], R["mech"]
    g, sg, ts, w, z = r["geo"], r["sg"], r["ts"], r["water"], r["zones"]
    ends = R.get("ends") or {}
    inl, out = curve.iloc[0], curve.iloc[-1]
    U = zonal_U(r)
    t_shell, t_min = shell_thickness(m["t_shell_req"], g["Ds"])
    noz = [int(x) for x in m["nozzles"]["NPS (inç)"]]
    tema = TEMA_CODE.get(inp.get("head_type"), "BEM")
    T_d_shell = math.ceil((r["T_hot"][0] + 30) / 5) * 5
    T_d_tube = math.ceil((max(w["T_in"], w["T_out"]) + 30) / 5) * 5
    d = dict(
        ch=ch, proj=inp.get("proj", {}), tema=tema, shells=r["shells"],
        Ds=g["Ds"], L=g["L"], A=r["A"], A_req=r["A_req"], OD=r["over_design"],
        # performans
        m_sh=inl.m_V + inl.m_L, mV_in=inl.m_V, mV_out=out.m_V, mL_in=inl.m_L, mL_out=out.m_L,
        nc_in=ends.get("nc_in"), nc_out=ends.get("nc_out"), m_w=w["m"],
        T_sh=(r["T_hot"][0], r["T_hot"][1]), T_t=(w["T_in"], w["T_out"]),
        rhoV=(inl.rho_V, out.rho_V), muV=(inl.mu_V * 1000, out.mu_V * 1000), cpV=(inl.cp_V, out.cp_V), kV=(inl.k_V, out.k_V),
        rhoL_out=out.rho_L if out.m_L > 0 else None, muL_out=out.mu_L * 1000 if out.m_L > 0 else None,
        cpL_out=out.cp_L if out.m_L > 0 else None, kL_out=out.k_L if out.m_L > 0 else None,
        rho_w=w["rho"], mu_w=w["mu"] * 1000, cp_w=w["cp"], k_w=w["k"],
        lam=latent_heat(curve), T_dew=R["T_dew"],
        P_sh_in=inp["P_in_bar"] * 100, P_t_in=inp["P_w_bar"] * 100,
        v_sh=(inl.m_V / 3600 / sg["As"] / inl.rho_V, out.m_V / 3600 / sg["As"] / out.rho_V), v_t=ts["v"],
        dp_s=(inp["dp_shell_max"], r["dp_shell_kPa"]), dp_t=(inp["dp_tube_max"], r["dp_tube_kPa"]),
        Rf=(r["Rf_shell"] * 1000, w["Rf"] * 1000), h=(r["h_shell"], ts["h"]),
        Q=r["Q_kW"], **U,
        # yapı
        P_d=(inp["P_d_shell"] * 100, inp["P_d_tube"] * 100), P_test=(1.3 * inp["P_d_shell"] * 100, 1.3 * inp["P_d_tube"] * 100),
        T_d=(T_d_shell, T_d_tube), passes=(r["shells"], g["passes"]), CA=(inp["CA_shell"], ch["CA_tube"]),
        noz=noz, rhov2_in=float(m["nozzles"].iloc[0]["ρv² (kg/m·s²)"]), rhov2_bundle=_bundle_rhov2(inl, out, sg),
        # konfigürasyon
        N=g["N"], Do=g["Do"], t=g["t"], bwg=BWG.get(round(g["t"], 2)), Pt=g["Pt"], layout=LAYOUT_DEG.get(g["layout"], g["layout"]),
        Db=sg["Db"] * 1000, n_baffle=sg["n_baffle"], B=g["B"],
        exp_joint=(tema == "BEM" and "körüğü" in m["tema"]),
        impingement=m["impingement"].startswith("GEREKLİ"),
        t_shell=t_shell, t_shell_min=t_min, t_shell_req=m["t_shell_req"], W=weights(g, t_shell),
        T_metal=(m["T_shell_metal"], m["T_tube_metal"]), span=(m["span_mm"], m["span_max_mm"]),
        t_ts=(m.get("tubesheet") or {}).get("T_sel_mm"), t_ts_req=(m.get("tubesheet") or {}).get("T_req_mm"),
        vib=(m.get("vib") or {}).get("V_Vc_max"), fn_min=(m.get("vib") or {}).get("fn_min"),
        # notlar
        n_zones=len(z), h_method=r["h_method"], thermo=inp.get("thermo_note", ""), comps=R.get("comps") or (),
        clearance=R.get("clearance", 15.0),
        shell_method=r.get("shell_method", "kern"), n_ss=(r.get("bd") or {}).get("n_ss"),
    )
    return d


# =====================================================================
# Biçimlendirme
# =====================================================================
def num(x, dec=1, lang="tr"):
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return "–"
    s = f"{x:,.{dec}f}"
    return s.replace(",", " ").replace(".", ",") if lang == "tr" else s


def yn(b, lang="tr"):
    return ("E" if b else "H") if lang == "tr" else ("Y" if b else "N")


# =====================================================================
# 2) Föy düzeni: bloklar -> satırlar -> hücreler (PDF, Excel, HTML aynı modeli çizer)
#    hücre = (metin, tür)  tür: n satır no · l etiket · b kalın etiket · c değer · bc kalın başlık · h bölüm · img
#    blok  = dict(widths=[...], rows=[[hücre, ...]], spans=[(satır, sütun1, sütun2)], light=[satır])
# =====================================================================
W8 = [16, 128, 66, 58, 66, 66, 74, 65]               # pt, toplam 539
W11 = [16, 92, 36, 44, 44, 52, 48, 44, 44, 52, 67]


class Block:
    def __init__(self, widths=W8):
        self.widths, self.rows, self.spans, self.light = widths, [], [], []

    def add(self, cells, spans=(), light=False):
        i = len(self.rows)
        self.rows.append(cells + [("", "c")] * (len(self.widths) - len(cells)))
        self.spans += [(i, a, b) for a, b in spans]
        if light:
            self.light.append(i)

    def section(self, no, title):
        self.add([(str(no), "n"), (title, "h")], [(1, 7)])

    def perf(self, no, label, unit, si, so, ti, to, merge=False):
        self.add([(str(no), "n"), (label, "l"), ("", "l"), (unit, "c"), (si, "c"), (so, "c"), (ti, "c"), (to, "c")],
                 [(1, 2)] + ([(4, 5), (6, 7)] if merge else []))

    def triple(self, no, a, av, b, bv, c="", cv=""):
        self.add([(str(no), "n"), (a, "l"), (av, "c"), (b, "l"), ("", "l"), (bv, "c"), (c, "l"), (cv, "c")], [(3, 4)])

    def wide(self, no, label, value):
        self.add([(str(no), "n"), (label, "b"), (value, "l")], [(2, 7)])

    def text(self, no, value):
        self.add([(str(no), "n"), (value, "l")], [(1, 7)])


def _header(blk, d):
    ch, pj = d["ch"], d["proj"]
    for i, (a, b, c, e) in enumerate([
            ("Müşteri:", ch["customer"] or pj.get("project", ""), "İş No.:", ch["job_no"]),
            ("Adres:", ch["address"], "Ref. / PO No.:", ch["ref_no"]),
            ("Tesis konumu:", ch["location"] or pj.get("plant", ""), "Teklif No.:", ch["proposal_no"]),
            ("Ünitenin hizmeti:", pj.get("service", ""), "Tarih / Rev.:", f"{pj.get('date', '')}  /  {pj.get('rev', '')}"),
            ("Hazırlayan:", pj.get("author", ""), "Kalem No.:", pj.get("item", ""))]):
        blk.add([(str(1 + i), "n"), (a, "b"), (b, "l"), ("", "l"), ("", "l"), (c, "b"), (e, "l")], [(2, 4), (6, 7)])


def layout(d):
    """Föyün iki sayfası: [[blok, ...], [blok, ...]]"""
    ch, n = d["ch"], num
    nz = d["noz"]
    fl = f'Class {ch["flange_class"]} {ch["facing"]}'
    # ---------------- SAYFA 1
    p = Block()
    _header(p, d)
    p.add([("6", "n"), ("Boyut (Ds / L):", "b"), (f"{n(d['Ds'], 0)} / {n(d['L'], 0)} mm", "c"), ("TEMA tipi:", "b"),
           (d["tema"], "c"), ("Yatay", "c"), ("Bağlantı:", "b"), (f"{d['shells']} seri · 1 paralel", "c")])
    p.add([("7", "n"), ("Yüzey/ünite (etkin):", "b"), (f"{n(d['A'] * d['shells'])} m²", "c"), ("Kabuk/ünite:", "b"),
           (str(d["shells"]), "c"), ("Yüzey/kabuk:", "b"), (f"{n(d['A'])} m²", "c")], [(6, 7)])
    p.section(8, "BİR ÜNİTENİN PERFORMANSI")
    p.add([("9", "n"), ("Akışkan yerleşimi", "b"), ("", "l"), ("", "l"), ("Kabuk tarafı", "bc"), ("", "l"), ("Tüp tarafı", "bc")],
          [(1, 3), (4, 5), (6, 7)], light=True)
    p.add([("10", "n"), ("", "l"), ("", "l"), ("", "l"), ("Giriş", "bc"), ("Çıkış", "bc"), ("Giriş", "bc"), ("Çıkış", "bc")],
          [(1, 3)], light=True)
    mw = n(d["m_w"], 0)
    p.perf(11, "Akışkan adı", "", ch["shell_fluid"], "", ch["tube_fluid"], "", merge=True)
    p.perf(12, "Toplam akışkan debisi", "kg/h", n(d["m_sh"], 0), "", mw, "", merge=True)
    p.perf(13, "   Buhar", "kg/h", n(d["mV_in"], 0), n(d["mV_out"], 0), "–", "–")
    p.perf(14, "   Sıvı", "kg/h", n(d["mL_in"], 0), n(d["mL_out"], 0), "–", "–")
    p.perf(15, "   Su buharı (steam)", "kg/h", "–", "–", "–", "–")
    p.perf(16, "   Su", "kg/h", "–", "–", mw, mw)
    p.perf(17, "   Yoğuşmayan gaz", "kg/h", n(d["nc_in"], 0), n(d["nc_out"], 0), "–", "–")
    p.perf(18, "Sıcaklık", "°C", n(d["T_sh"][0]), n(d["T_sh"][1]), n(d["T_t"][0]), n(d["T_t"][1]))
    vl = lambda v, l, dv=2, dl=0: f"{n(v, dv)} / {n(l, dl) if l is not None else '–'}"
    p.perf(19, "Yoğunluk (buhar / sıvı)", "kg/m³", vl(d["rhoV"][0], None), vl(d["rhoV"][1], d["rhoL_out"]),
           f"– / {n(d['rho_w'], 0)}", f"– / {n(d['rho_w'], 0)}")
    p.perf(20, "Viskozite (buhar / sıvı)", "mPa·s", vl(d["muV"][0], None, 4), vl(d["muV"][1], d["muL_out"], 4, 3),
           f"– / {n(d['mu_w'], 3)}", f"– / {n(d['mu_w'], 3)}")
    p.perf(21, "Özgül ısı (buhar / sıvı)", "kJ/(kg·°C)", vl(d["cpV"][0], None, 3), vl(d["cpV"][1], d["cpL_out"], 3, 3),
           f"– / {n(d['cp_w'], 3)}", f"– / {n(d['cp_w'], 3)}")
    p.perf(22, "Isıl iletkenlik (buhar / sıvı)", "W/(m·°C)", vl(d["kV"][0], None, 4), vl(d["kV"][1], d["kL_out"], 4, 3),
           f"– / {n(d['k_w'], 3)}", f"– / {n(d['k_w'], 3)}")
    p.perf(23, "Gizli ısı (yoğuşan kütleye göre ort.)", "kJ/kg @ °C",
           f"{n(d['lam'], 0)} @ {n(d['T_dew'])}–{n(d['T_sh'][1])}", "", "–", "", merge=True)
    p.perf(24, "Giriş basıncı (mutlak)", "kPa(a)", n(d["P_sh_in"], 0), "", n(d["P_t_in"], 0), "", merge=True)
    p.perf(25, "Hız", "m/s", n(d["v_sh"][0]), n(d["v_sh"][1]), n(d["v_t"], 2), n(d["v_t"], 2))
    p.perf(26, "Basınç düşüşü, izin verilen / hesaplanan", "kPa", f"{n(d['dp_s'][0])} / {n(d['dp_s'][1])}", "",
           f"{n(d['dp_t'][0])} / {n(d['dp_t'][1])}", "", merge=True)
    p.perf(27, "Kirlenme direnci (min.)", "(m²·°C)/kW", n(d["Rf"][0], 2), "", n(d["Rf"][1], 2), "", merge=True)
    p.perf(28, "Ortalama film katsayısı", "W/(m²·°C)", n(d["h"][0], 0), "", n(d["h"][1], 0), "", merge=True)
    p.add([("29", "n"), ("Aktarılan ısı:", "b"), (f"{n(d['Q'], 0)} kW", "c"), ("MTD (bölgesel, düzeltilmiş):", "b"),
           ("", "l"), ("", "l"), (f"{n(d['mtd'])} °C", "c")], [(3, 5), (6, 7)])
    p.add([("30", "n"), ("U, gerekli:", "b"), (f"{n(d['U_req'], 0)} W/(m²·°C)", "c"), ("Kirli:", "b"), (n(d["U_dirty"], 0), "c"),
           ("Temiz:", "b"), (n(d["U_clean"], 0), "c"), (f"Aşırı tas. %{n(d['OD'])}", "c")])
    p.section(31, "BİR KABUĞUN YAPISI")
    p.add([("32", "n"), ("", "l"), ("", "l"), ("", "l"), ("Kabuk tarafı", "bc"), ("", "l"), ("Tüp tarafı", "bc")],
          [(1, 3), (4, 5), (6, 7)], light=True)
    p.perf(33, "Tasarım basıncı iç / dış / test", "kPa(g)", f"{n(d['P_d'][0], 0)} / – / {n(d['P_test'][0], 0)}", "",
           f"{n(d['P_d'][1], 0)} / – / {n(d['P_test'][1], 0)}", "", merge=True)
    p.perf(34, "Tasarım sıcaklığı / MDMT", "°C", f"{n(d['T_d'][0], 0)} / {ch['mdmt'] or '–'}", "",
           f"{n(d['T_d'][1], 0)} / {ch['mdmt'] or '–'}", "", merge=True)
    p.perf(35, "Geçiş sayısı", "", "1", "", str(d["passes"][1]), "", merge=True)
    p.perf(36, "Korozyon payı", "mm", n(d["CA"][0]), "", n(d["CA"][1]), "", merge=True)
    p.perf(37, "Gelecek işleme payı (ayna, flanş, kapak)", "mm", "–", "", "–", "", merge=True)
    p.perf(38, "Asgari radyografi", "", ch["radiography"], "", ch["radiography"], "", merge=True)
    p.perf(39, "Özel hizmet (öldürücü, ekşi vb.)", "", ch["special_service"], "", "Yok", "", merge=True)
    p.perf(40, "Kaynak sonrası ısıl işlem (PWHT)", "", "Evet" if ch["pwht"] else "Hayır", "", "Evet" if ch["pwht"] else "Hayır", "", merge=True)
    conn = [("Giriş", f'1 × {nz[0]}" {fl}', f'1 × {nz[2]}" {fl}'), ("Çıkış", f'1 × {nz[1]}" {fl}', f'1 × {nz[2]}" {fl}'),
            ("Ara", "–", "–"), ("Havalık", f'1 × {ch["vent_shell"]} {fl}', f'1 × {ch["vent_tube"]} {fl}'),
            ("Drenaj", f'1 × {ch["drain_shell"]} {fl}', f'1 × {ch["drain_tube"]} {fl}'), ("Diğer", "–", "–")]
    for k, (lab, s_, t_) in enumerate(conn):
        p.add([(str(41 + k), "n"), (["Bağlantılar:", "Adet / boyut / sınıf,", "tip ve yüzey", "", "", ""][k], "l"), (lab, "l"),
               ("", "l"), (s_, "c"), ("", "c"), (t_, "c")], [(2, 3), (4, 5), (6, 7)])
    p.text(47, "")
    p.section(48, "KABUK BAŞINA EŞANJÖR KONFİGÜRASYONU")
    p.triple(49, "Tüp sayısı:", str(d["N"]), "Tüp dış çapı (mm):", n(d["Do"], 2), "Tüp eti (mm) / BWG:",
             f"{n(d['t'], 2)} / {d['bwg'] or '–'}")
    p.triple(50, "Tüp tipi (düz/kanatlı/diğer):", "Düz", "Tüp adımı (mm):", n(d["Pt"], 2), "Tüp eti min./ort.:", "Min.")
    p.triple(51, "Tüp boyu (m):", n(d["L"] / 1000, 3), "Tüp düzeni (30/60/90/45°):", d["layout"])
    p.triple(52, "Tüp–ayna bağlantısı:", ch["tube_joint"], "Kabuk çapı (mm), iç/dış:", f"{n(d['Ds'], 0)} (iç)", "Gövde flanşı tipi:",
             "Yok" if d["tema"] == "BEM" else "Kaynak boyunlu")
    p.triple(53, "Çapraz baffle tipi:", "Tek segmental", "Çapraz baffle sayısı:", str(d["n_baffle"]), "Baffle aralığı (mm):", n(d["B"], 0))
    p.triple(54, "Baffle kesimi (% çap):", str(ch["baffle_cut"]), "Dikey / yatay kesim:", ch["cut_dir"], "Giriş aralığı (mm):", n(d["B"], 0))
    p.triple(55, "Boyuna baffle sızdırmazlığı:", "Yok", "Tüp destek tipi:", "–", "U-büküm desteği:", "Yok" if d["tema"] != "BEU" else "Var")
    p.triple(56, "Genleşme körüğü gerekli (E/H):", yn(d["exp_joint"]), "Körük tipi:", "Teyit edilecek" if d["exp_joint"] else "Yok",
             "Bypass sızdırmazlığı:", "Teyit edilecek")
    p.triple(57, "Çarpma koruması gerekli (E/H):", yn(d["impingement"]), "Koruma tipi:", "Plaka" if d["impingement"] else "Yok",
             "Yedek conta takımı:", "1")
    p.triple(58, "Mesnet tipi:", ch["supports"], "Mesnet sayısı:", str(ch["n_supports"]))
    p.triple(59, "ρV² giriş nozulu (kg/(m·s²)):", n(d["rhov2_in"], 0), "ρV² demet girişi:", n(d["rhov2_bundle"][0], 0),
             "ρV² demet çıkışı:", n(d["rhov2_bundle"][1], 0))
    p.triple(60, "Toplam ağırlık/kabuk, boş (kg):", "≈ " + n(d["W"]["empty"], 0), "Su dolu (kg):", "≈ " + n(d["W"]["full"], 0),
             "Demet (kg):", "≈ " + n(d["W"]["bundle"], 0))
    hm = {"colburn_hougen": "Colburn–Hougen", "sbg": "SBG", "fixed": "sabit h"}.get(d["h_method"], d["h_method"])
    for k, t in enumerate([
            f"Notlar: Isıl tasarım {d['n_zones']} bölgeli; kabuk h: {'Bell–Delaware' if d['shell_method'] == 'bell_delaware' else 'Kern'} (gaz) + Nusselt (film) + {hm}; tüp h: Gnielinski. "
            f"Gerekli alan {n(d['A_req'])} m².",
            f"Kabuk ΔP: {'Bell–Delaware' if d['shell_method'] == 'bell_delaware' else 'Kern'} + iki faz çarpanı (C = 8), bölgesel toplam + nozullar. Termodinamik: {d['thermo'] or '–'}.",
            f"Kabuk eti: hesap {n(d['t_shell_req'], 1)} mm (KP dahil), TEMA asgari {n(d['t_shell_min'], 1)} mm → {d['t_shell']} mm. "
            + (f"Ayna {d['t_ts']} mm (A.131). " if d.get("t_ts") else "")
            + (f"Titreşim V/Vc maks. {n(d['vib'], 2)}. " if d.get("vib") is not None else "")
            + "Ağırlıklar yaklaşıktır.",
            "ρV² demet değerleri TEMA RCB-4.624 (darbe plakasız) alanıyla; sınır 5953 (RCB-4.62). MDMT hesaplanmamıştır."]):
        p.text(61 + k, t)
    page1 = [p]

    # ---------------- SAYFA 2
    a = Block()
    _header(a, d)
    a.text(6, "")
    a.section(7, "EŞANJÖR ŞEMASI")
    a.add([("8–27", "n"), ("", "img")], [(1, 7)])
    a.section(28, "TASARIM VE KOD GEREKSİNİMLERİ")
    a.triple(29, "Tasarım kodu:", ch["code"], "Kod damgası gerekli:", yn(ch["code_stamp"]), "MAWP hesabı (E/H):", yn(ch["mawp"]))
    a.triple(30, "TEMA sınıfı (R, C, B):", ch["tema_class"], "Müşteri muayenesi (E/H):", yn(ch["inspection"]), "Çevrimsel hiz. (E/H):", "H")
    a.triple(31, "Başlangıç UT (E/H):", "H", "Başlangıç girdap akımı (E/H):", "H")
    a.section(33, "TASARIM, İŞLETME VE ARIZA KOŞULLARI")
    c = Block(W11)
    c.add([("34", "n"), ("", "l"), ("", "l"), ("Kabuk akışkan sıc.", "bc"), ("", "l"), ("Basınç", "bc"), ("Kab. metal", "bc"),
           ("Tüp akışkan sıc.", "bc"), ("", "l"), ("Basınç", "bc"), ("Tüp metal", "bc")], [(3, 4), (7, 8)], light=True)
    c.add([("35", "n"), ("Durum", "bc"), ("Çevrim", "bc"), ("Giriş (°C)", "bc"), ("Çıkış (°C)", "bc"), ("(kPa)", "bc"),
           ("sıc. (°C)", "bc"), ("Giriş (°C)", "bc"), ("Çıkış (°C)", "bc"), ("(kPa)", "bc"), ("sıc. (°C)", "bc")], light=True)
    c.add([("36", "n"), ("Tasarım", "l"), ("–", "c"), (n(d["T_d"][0], 0), "c"), ("–", "c"), (f"{n(d['P_d'][0], 0)} (g)", "c"), ("–", "c"),
           (n(d["T_d"][1], 0), "c"), ("–", "c"), (f"{n(d['P_d'][1], 0)} (g)", "c"), ("–", "c")])
    c.add([("37", "n"), ("Normal işletme", "l"), ("–", "c"), (n(d["T_sh"][0], 0), "c"), (n(d["T_sh"][1], 0), "c"),
           (f"{n(d['P_sh_in'], 0)} (a)", "c"), (n(d["T_metal"][0], 0), "c"), (n(d["T_t"][0], 0), "c"), (n(d["T_t"][1], 0), "c"),
           (f"{n(d['P_t_in'], 0)} (a)", "c"), (n(d["T_metal"][1], 0), "c")])
    for k in range(38, 43):
        c.add([(str(k), "n")])
    b = Block()
    b.section(44, "MALZEMELER")
    b.add([("45", "n"), ("Tüpler:", "b"), (ch["mat_tube"], "l"), ("Kabuk tarafı", "bc"), ("", "l"), ("", "l"), ("Tüp tarafı", "bc")],
          [(3, 5), (6, 7)])
    fixed = d["tema"] in ("BEM",)
    fh = d["tema"] in ("AES", "AET")
    mats = [("Sabit ayna:", ch["mat_tubesheet"], "Kabuk:", ch["mat_shell"], "Kanal / bonnet:", ch["mat_channel"]),
            ("Yüzer ayna:", ch["mat_tubesheet"] if fh else "Yok", "Kabuk kapağı:", ch["mat_shell"] if fh else "Yok",
             "Kanal kapağı:", "Yok (bonnet)" if d["tema"].startswith("B") else ch["mat_channel"]),
            ("Yüzer kafa kapağı:", ch["mat_shell"] if fh else "Yok", "Gövde flanşları:", "Yok" if fixed else ch["mat_flange"],
             "Gövde flanşları:", ch["mat_flange"]),
            ("Çapraz baffle'lar:", ch["mat_baffle"], "Contalar:", "Yok" if fixed else ch["gasket"], "Contalar:", ch["gasket"]),
            ("Boyuna baffle:", "Yok", "Nozul flanşları:", f"{ch['mat_flange']} WN {ch['facing']}", "Nozul flanşları:",
             f"{ch['mat_flange']} WN {ch['facing']}"),
            ("Bağlantı çubukları, aralayıcılar:", ch["mat_tierod"], "Nozul boyunları:", ch["mat_neck"], "Nozul boyunları:", ch["mat_neck"]),
            ("Mesnetler:", ch["mat_tierod"], "Genleşme körüğü:", "Teyit edilecek" if d["exp_joint"] else "Yok")]
    for k, row in enumerate(mats):
        b.triple(46 + k, *row)
    for k, (lab, v) in enumerate([("İzolasyon tipi / kalınlığı:", ch["insulation"]), ("Dış cıvatalar:", ch["bolting"]),
                                  ("İç cıvatalar:", ch["bolting"] if fh else "Yok"), ("Boya gereksinimleri:", ch["paint"])]):
        b.wide(53 + k, lab, v)
    b.text(57, "")
    b.add([("58", "n"), ("Diğer gereksinimler ve şartnameler:", "b"), ("", "l"),
           (f"TEMA 9. baskı (2007), Sınıf {ch['tema_class']}; {ch['code']}; ASME B16.5 flanşlar.", "l")], [(1, 2), (3, 7)])
    b.text(59, "Isıl tasarım öğrenci tasarım aracıyla yapılmıştır; yayından önce DWSIM ile doğrulanmalıdır.")
    return [page1, [a, c, b]]


# =====================================================================
# 3) Şema (yan görünüş + ayna yerleşimi)
# =====================================================================
def sketch_png(d, lang="tr", dpi=200):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle, FancyBboxPatch, Circle, Polygon
    from matplotlib.collections import PatchCollection
    import esanjor_core as ec

    T = (lambda tr, en_: tr) if lang == "tr" else (lambda tr, en_: en_)
    nf = (lambda x, dec=0: num(x, dec, lang))
    INK, MUTED, BLUE, ORANGE = "#1f1f1d", "#6b6b66", "#2a78d6", "#eb6834"
    with plt.rc_context({"font.family": "DejaVu Sans", "font.size": 6.5}):
        fig = plt.figure(figsize=(494 / 72, 212 / 72), dpi=dpi)
        ax = fig.add_axes([0.0, 0.0, 0.70, 1.0]); ax.set_xlim(0, 100); ax.set_ylim(0, 60); ax.axis("off")
        x0, x1, yb, yt = 25, 82, 18, 42
        ax.add_patch(Rectangle((x0, yb), x1 - x0, yt - yb, fill=False, lw=0.9, ec=INK))
        for x in (x0, x1):
            ax.add_patch(Rectangle((x - 0.8, yb - 2), 1.6, yt - yb + 4, fc="#d9d9d6", ec=INK, lw=0.6))
        ax.add_patch(FancyBboxPatch((x0 - 11, yb), 10, yt - yb, boxstyle="round,pad=0,rounding_size=4", fill=False, lw=0.9, ec=INK))
        rear_round = d["tema"][2] in "MU"
        ax.add_patch(FancyBboxPatch((x1 + 1, yb), 10, yt - yb, boxstyle="round,pad=0,rounding_size=4" if rear_round else "square,pad=0",
                                    fill=False, lw=0.9, ec=INK))
        npass = d["passes"][1]
        front, rear = (npass // 2, npass // 2 - 1) if npass > 1 else (0, 0)
        for k in range(front):
            y = yb + (yt - yb) * (k + 1) / (front + 1)
            ax.plot([x0 - 11, x0 - 1], [y, y], color=INK, lw=0.6)
        for k in range(max(rear, 0)):
            y = yb + (yt - yb) * (k + 1) / (rear + 1)
            ax.plot([x1 + 1, x1 + 11], [y, y], color=INK, lw=0.6)
        ax.text(x0 - 6, 30 if front % 2 == 0 else 35.5, d["tema"][0], ha="center", va="center", fontsize=7, fontweight="bold")
        ax.text(62, 47.5, d["tema"][1], ha="center", fontsize=7, fontweight="bold")
        ax.text(x1 + 6, 30 if rear % 2 == 0 else 35.5, d["tema"][2], ha="center", va="center", fontsize=7, fontweight="bold")
        nb = d["n_baffle"]
        for i in range(1, nb + 1):
            x = x0 + i * (x1 - x0) / (nb + 1)
            ax.plot([x, x], [yb + 1, yt - 1], color=MUTED, lw=0.7, ls=(0, (3, 2)))
        for y in (24, 30, 36):
            ax.plot([x0, x1], [y, y], color="#b9b9b4", lw=0.5)

        def nozzle(x, y, up, label, col, ha="center", dx=0):
            h = 5 if up else -5
            ax.add_patch(Rectangle((x - 1.6, y if up else y + h), 3.2, abs(h), fc="white", ec=INK, lw=0.7))
            ax.annotate(label, (x, y + h), xytext=(dx, 4 if up else -4), textcoords="offset points",
                        ha=ha, va="bottom" if up else "top", fontsize=6.2, color=col, fontweight="bold")
        nz, ch = d["noz"], d["ch"]
        nozzle(35, yt, True, f'S1  {nz[0]}" ' + T("buhar girişi", "vapor inlet"), ORANGE, ha="left", dx=-5)
        nozzle(76, yb, False, f'S2  {nz[1]}" ' + T("kondensat\n+ gaz çıkışı", "condensate\n+ gas outlet"), ORANGE, ha="left", dx=-5)
        nozzle(77, yt, True, f'V  {ch["vent_shell"]} ' + T("havalık", "vent"), MUTED)
        nozzle(30, yb, False, f'D  {ch["drain_shell"]} ' + T("drenaj", "drain"), MUTED, ha="left", dx=-5)
        out_rear = npass % 2 == 1
        for y_, up_, lab_, xx in ((yt, True, f'T2  {nz[2]}"\n' + T("su çıkışı", "water out"), x1 + 9.6 if out_rear else x0 - 6),
                                  (yb, False, f'T1  {nz[2]}"\n' + T("su girişi", "water in"), x0 - 6)):
            h_ = 5 if up_ else -5
            ax.add_patch(Rectangle((xx - 1.6, y_ if up_ else y_ + h_), 3.2, abs(h_), fc="white", ec=INK, lw=0.7))
            ax.text(xx - 2.6, y_ + h_ / 2, lab_, ha="right", va="center", fontsize=6.2, color=BLUE, fontweight="bold")
        for x in (x0 + 0.35 * (x1 - x0), x0 + 0.68 * (x1 - x0)):
            ax.add_patch(Polygon([[x - 4, yb - 7.5], [x + 4, yb - 7.5], [x + 2.5, yb], [x - 2.5, yb]], closed=True, fc="#e6e6e3", ec=INK, lw=0.6))
        ax.annotate("", (x0, 7), xytext=(x1, 7), arrowprops=dict(arrowstyle="<->", lw=0.6, color=INK))
        ax.text((x0 + x1) / 2, 4.2, f"{nf(d['L'])} mm " + T("(tüp boyu)", "(tube length)"), ha="center", fontsize=6.2)
        ax.annotate("", (x1 + 14, yb), xytext=(x1 + 14, yt), arrowprops=dict(arrowstyle="<->", lw=0.6, color=INK))
        ax.text(x1 + 15, 30, T("İç çap ", "ID ") + nf(d["Ds"]), va="center", fontsize=6.2, rotation=90)
        ax.text(53, 57, f"TEMA {d['tema']}  ·  {nf(d['Ds'])} × {nf(d['L'])}  ·  " +
                T(f"1 kabuk geçişi / {npass} tüp geçişi  ·  yatay", f"1 shell pass / {npass} tube passes  ·  horizontal"),
                ha="center", fontsize=7, fontweight="bold", color=INK)
        ax.text(53, 52.5, T(f"{nb} tek segmental baffle, {nf(d['B'])} mm aralık, %{ch['baffle_cut']} kesim  ·  2 semer mesnet",
                            f"{nb} single-segmental baffles, {nf(d['B'])} mm spacing, {ch['baffle_cut']}% cut  ·  2 saddles"),
                ha="center", fontsize=6.2, color=MUTED)
        ax2 = fig.add_axes([0.71, 0.08, 0.28, 0.84])
        Ds = d["Ds"]
        layout_ = {"30°": "triangular", "90°": "square", "45°": "rotated square"}.get(d["layout"], "triangular")
        nmax, pts, lanes = ec.tube_layout(Ds - d["clearance"], d["Do"], d["Pt"], layout_, npass)
        order = np.argsort(np.hypot(pts[:, 0], pts[:, 1]))
        used = np.zeros(nmax, bool); used[order[:min(d["N"], nmax)]] = True
        ax2.add_collection(PatchCollection([Circle(p_, d["Do"] / 2) for p_ in pts[used]], fc=BLUE, ec="none"))
        ax2.add_patch(Circle((0, 0), Ds / 2, fill=False, lw=0.8, ec=INK))
        for y in lanes["h"]:
            ax2.plot([-Ds / 2, Ds / 2], [y, y], color=ORANGE, lw=0.6)
        for x in lanes["v"]:
            ax2.plot([x, x], [-Ds / 2, Ds / 2], color=ORANGE, lw=0.6)
        lim = Ds / 2 * 1.045
        ax2.set_xlim(-lim, lim); ax2.set_ylim(-lim, lim); ax2.set_aspect("equal"); ax2.axis("off")
        ax2.set_title(T(f"Ayna: {d['N']} tüp × Ø{nf(d['Do'], 2)} mm, {d['layout']}, adım {nf(d['Pt'], 2)}",
                        f"Tubesheet: {d['N']} tubes × Ø{nf(d['Do'], 2)} mm, {d['layout']}, pitch {nf(d['Pt'], 2)}"), fontsize=5.8, pad=2)
        buf = io.BytesIO()
        fig.savefig(buf, format="png", dpi=dpi)
        plt.close(fig)
    return buf.getvalue()


# =====================================================================
# 4) Türkçe PDF
# =====================================================================
def _fonts():
    import matplotlib
    from reportlab.pdfbase import pdfmetrics
    from reportlab.pdfbase.ttfonts import TTFont
    if "DV" not in pdfmetrics.getRegisteredFontNames():
        fd = os.path.join(os.path.dirname(matplotlib.__file__), "mpl-data", "fonts", "ttf")
        pdfmetrics.registerFont(TTFont("DV", os.path.join(fd, "DejaVuSans.ttf")))
        pdfmetrics.registerFont(TTFont("DVB", os.path.join(fd, "DejaVuSans-Bold.ttf")))


def pdf_tr(d):
    from reportlab.lib import colors
    from reportlab.lib.pagesizes import A4
    from reportlab.lib.styles import ParagraphStyle
    from reportlab.platypus import SimpleDocTemplate, Table, TableStyle, Paragraph, Spacer, Image, PageBreak
    from xml.sax.saxutils import escape
    _fonts()
    FS = 6.2
    st = {"l": ParagraphStyle("l", fontName="DV", fontSize=FS, leading=FS + 1)}
    st["c"] = ParagraphStyle("c", parent=st["l"], alignment=1)
    st["b"] = ParagraphStyle("b", parent=st["l"], fontName="DVB")
    st["bc"] = ParagraphStyle("bc", parent=st["b"], alignment=1)
    st["h"] = ParagraphStyle("h", fontName="DVB", fontSize=FS + 0.4, leading=FS + 2, alignment=1, textColor=colors.white)
    st["n"] = ParagraphStyle("n", parent=st["c"], textColor=colors.HexColor("#6b6b66"))
    DARK, LIGHT, GRID = colors.HexColor("#2b2b29"), colors.HexColor("#eeeeea"), colors.HexColor("#9e9e99")
    img = Image(io.BytesIO(sketch_png(d, "tr")), width=505, height=214)

    def table(blk):
        data, sty = [], []
        for i, row in enumerate(blk.rows):
            out = []
            for j, (t, k) in enumerate(row):
                if k == "img":
                    out.append(img)
                elif k == "h":
                    out.append(Paragraph(escape(t), st["h"]))
                    sty.append(("BACKGROUND", (j, i), (-1, i), DARK))
                else:
                    out.append(Paragraph(escape(str(t)), st[k]) if t else "")
            data.append(out)
        sty += [("SPAN", (a, i), (b_, i)) for i, a, b_ in blk.spans]
        sty += [("BACKGROUND", (1, i), (-1, i), LIGHT) for i in blk.light]
        rh = [None if any(k == "img" for _, k in row) else 11.2 for row in blk.rows]
        t = Table(data, colWidths=blk.widths, rowHeights=rh)
        t.setStyle(TableStyle([("GRID", (0, 0), (-1, -1), 0.4, GRID), ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
                               ("TOPPADDING", (0, 0), (-1, -1), 0), ("BOTTOMPADDING", (0, 0), (-1, -1), 0.8),
                               ("LEFTPADDING", (0, 0), (-1, -1), 2.5), ("RIGHTPADDING", (0, 0), (-1, -1), 2.5)] + sty))
        return t

    def title():
        return [Paragraph("TEMA ISI DEĞİŞTİRİCİ SPESİFİKASYON FÖYÜ",
                          ParagraphStyle("t", fontName="DVB", fontSize=10.5, alignment=1, leading=12)),
                Paragraph("TEMA Heat Exchanger Specification Sheet formatının Türkçe karşılığı · SI birimleri",
                          ParagraphStyle("st", fontName="DV", fontSize=6.5, alignment=1, textColor=colors.HexColor("#6b6b66"), leading=9)),
                Spacer(1, 2)]

    item = d["proj"].get("item", "") or "Eşanjör"

    def footer(canvas, doc):
        canvas.setFont("DV", 6.5)
        canvas.setFillColor(colors.HexColor("#6b6b66"))
        canvas.drawRightString(A4[0] - 28, 18, f"sayfa {doc.page} / 2")
        canvas.drawString(28, 18, f"{item} · TEMA spesifikasyon föyü · uygulama hesaplarından otomatik doldurulmuştur")

    p1, p2 = layout(d)
    buf = io.BytesIO()
    doc = SimpleDocTemplate(buf, pagesize=A4, leftMargin=28, rightMargin=28, topMargin=18, bottomMargin=26,
                            title=f"{item} TEMA Spesifikasyon Föyü")
    story = title() + [table(b) for b in p1] + [PageBreak()] + title() + [table(b) for b in p2]
    doc.build(story, onFirstPage=footer, onLaterPages=footer)
    return buf.getvalue()


# =====================================================================
# 5) Tüm hesaplanan veriler (Excel 2. sayfa) — kaynak sütunuyla
# =====================================================================
def data_table(d, R):
    r, c, m = R["r"], R["cost"], R["mech"]
    lf, ts = r["lf"], r["ts"]
    H, G, S = "Hesap", "Girdi", "Seçim"
    rows = [
        ("Isıl", "Isı yükü Q", d["Q"], "kW", H), ("Isıl", "Çiy noktası", d["T_dew"], "°C", H),
        ("Isıl", "Ortalama gizli ısı", d["lam"], "kJ/kg", H), ("Isıl", "Soğutma suyu debisi", d["m_w"], "kg/h", H),
        ("Isıl", "LMTD (tek)", lf["LMTD"], "°C", H), ("Isıl", "F düzeltme faktörü", lf["F"], "-", H),
        ("Isıl", "R / P", f"{lf['R']:.3f} / {lf['P']:.3f}", "-", H),
        ("Isıl", "MTD (bölgesel, F dahil)", d["mtd"], "°C", H), ("Isıl", "Bölge sayısı", d["n_zones"], "-", G),
        ("Isıl", "h kabuk (alan ağırlıklı)", d["h"][0], "W/m²K", H), ("Isıl", "h kabuk (en düşük bölge)", r["h_shell_min"], "W/m²K", H),
        ("Isıl", "h tüp", d["h"][1], "W/m²K", H), ("Isıl", "Tüp Re / Pr / Nu", f"{ts['Re']:.0f} / {ts['Pr']:.2f} / {ts['Nu']:.0f}", "-", H),
        ("Isıl", "U temiz", d["U_clean"], "W/m²K", H), ("Isıl", "U kirli", d["U_dirty"], "W/m²K", H),
        ("Isıl", "U gerekli (tasarım alanına göre)", d["U_req"], "W/m²K", H),
        ("Isıl", "Gerekli alan", d["A_req"], "m²", H), ("Isıl", "Tasarım alanı", d["A"], "m²", H),
        ("Isıl", "Aşırı tasarım", d["OD"], "%", H),
        ("Isıl", "ε-NTU alanı (ters akış, F=1)", r["A_eNTU_cf"], "m²", H), ("Isıl", "LMTD alanı (ters akış, F=1)", r["A_LMTD_cf"], "m²", H),
        ("Hidrolik", "Tüp hızı", d["v_t"], "m/s", H), ("Hidrolik", "Kabuk gaz hızı giriş / çıkış", f"{d['v_sh'][0]:.1f} / {d['v_sh'][1]:.1f}", "m/s", H),
        ("Hidrolik", "Tüp ΔP (hesap / sınır)", f"{d['dp_t'][1]:.1f} / {d['dp_t'][0]:.1f}", "kPa", H),
        ("Hidrolik", "Kabuk ΔP (hesap / sınır)", f"{d['dp_s'][1]:.1f} / {d['dp_s'][0]:.1f}", "kPa", H),
        ("Hidrolik", "Kabuk ΔP yalnız gaz", r["dp_shell_gas_kPa"], "kPa", H),
        ("Geometri", "Tüp sayısı", d["N"], "-", G), ("Geometri", "Tüp Do × et", f"{d['Do']} × {d['t']}", "mm", G),
        ("Geometri", "Tüp boyu", d["L"], "mm", G), ("Geometri", "Geçiş sayısı", d["passes"][1], "-", G),
        ("Geometri", "Adım / dizilim", f"{d['Pt']} / {d['layout']}", "mm", G), ("Geometri", "Kabuk iç çapı", d["Ds"], "mm", G),
        ("Geometri", "Demet çapı (OTL)", d["Db"], "mm", H), ("Geometri", "Baffle aralığı / sayısı", f"{d['B']:.0f} / {d['n_baffle']}", "mm / -", G),
        ("Mekanik", "TEMA tipi", d["tema"], "-", H), ("Mekanik", "Metal sıcaklığı kabuk / tüp", f"{d['T_metal'][0]:.0f} / {d['T_metal'][1]:.0f}", "°C", H),
        ("Mekanik", "Nozullar kabuk giriş / çıkış / su", " / ".join(f'{x}"' for x in d["noz"]), "NPS", H),
        ("Mekanik", "ρv² kabuk giriş nozulu", d["rhov2_in"], "kg/m·s²", H),
        ("Mekanik", "Desteksiz tüp açıklığı / TEMA sınırı", f"{d['span'][0]:.0f} / {d['span'][1]:.0f}", "mm", H),
        ("Mekanik", "Kabuk eti hesap (KP dahil) / TEMA asgari / seçilen", f"{d['t_shell_req']:.2f} / {d['t_shell_min']} / {d['t_shell']}", "mm", H),
        ("Mekanik", "Ayna kalınlığı gerekli / seçilen (TEMA A.131–A.132, R-7.11)",
         f"{d['t_ts_req']:.1f} / {d['t_ts']}" if d.get("t_ts") else "-", "mm", H),
        ("Mekanik", "Titreşim: en büyük V/Vc (TEMA V-10) / en düşük fn",
         f"{d['vib']:.2f} / {d['fn_min']:.1f}" if d.get("vib") is not None else "-", "- / Hz", H),
        ("Mekanik", "Ağırlık boş / su dolu / demet", f"{d['W']['empty']:.0f} / {d['W']['full']:.0f} / {d['W']['bundle']:.0f}", "kg", H),
        ("Mekanik", "Tasarım basıncı kabuk / tüp", f"{d['P_d'][0]:.0f} / {d['P_d'][1]:.0f}", "kPa(g)", G),
        ("Mekanik", "Test basıncı kabuk / tüp (1,3 × tasarım)", f"{d['P_test'][0]:.0f} / {d['P_test'][1]:.0f}", "kPa(g)", H),
        ("Mekanik", "Tasarım sıcaklığı kabuk / tüp (işletme + 30 °C)", f"{d['T_d'][0]:.0f} / {d['T_d'][1]:.0f}", "°C", H),
        ("Maliyet", "Satın alma (2010)", c["Ce_2010"], "USD", H), ("Maliyet", "Satın alma (güncel)", c["Ce"], "USD", H),
        ("Maliyet", "Kurulum çarpanı", c["install_mult"], "-", G), ("Maliyet", "Kurulu yatırım", c["C_installed"], "USD", H),
        ("Maliyet", "CRF", c["crf"], "-", H), ("Maliyet", "Yıllık sermaye payı", c["annual_capital"], "USD/yıl", H),
        ("Maliyet", "Pompa gücü", c["P_pump_kW"], "kW", H), ("Maliyet", "Pompa elektriği", c["E_cost"], "USD/yıl", H),
        ("Maliyet", "Soğutma suyu", c["W_cost"], "USD/yıl", H), ("Maliyet", "Yıllık toplam maliyet (TAC)", c["TAC"], "USD/yıl", H),
    ]
    df = pd.DataFrame(rows, columns=["Grup", "Büyüklük", "Değer", "Birim", "Kaynak"])
    df["Değer"] = [round(float(v), 4) if isinstance(v, (int, float, np.floating, np.integer)) and not isinstance(v, bool) else v
                   for v in df["Değer"]]
    return df


# =====================================================================
# 6) Excel
# =====================================================================
def to_excel(d, R, extra_sheets=None):
    from openpyxl import Workbook
    from openpyxl.styles import Alignment, Border, Font, PatternFill, Side
    from openpyxl.utils import get_column_letter
    from openpyxl.drawing.image import Image as XLImage

    wb = Workbook()
    ws = wb.active
    ws.title = "TEMA föyü"
    pages = layout(d)
    # ortak sütun ızgarası: iki genişlik setinin sınırlarının birleşimi
    cum = lambda w: list(np.cumsum([0] + w))
    bounds = sorted(set(cum(W8)) | set(cum(W11)))
    for i in range(len(bounds) - 1):
        ws.column_dimensions[get_column_letter(i + 1)].width = (bounds[i + 1] - bounds[i]) / 5.2
    col_of = {tuple(w): [bounds.index(b) + 1 for b in cum(w)] for w in (W8, W11)}
    thin = Side(style="thin", color="9E9E99")
    border = Border(left=thin, right=thin, top=thin, bottom=thin)
    dark, light = PatternFill("solid", fgColor="2B2B29"), PatternFill("solid", fgColor="EEEEEA")
    fonts = {"l": Font(name="Arial", size=7.5), "c": Font(name="Arial", size=7.5), "b": Font(name="Arial", size=7.5, bold=True),
             "bc": Font(name="Arial", size=7.5, bold=True), "h": Font(name="Arial", size=7.5, bold=True, color="FFFFFF"),
             "n": Font(name="Arial", size=7, color="6B6B66"), "img": Font(name="Arial", size=7.5)}
    ncol = len(bounds) - 1
    row = 1
    breaks = []
    for pi, page in enumerate(pages):
        ws.cell(row, 1, "TEMA ISI DEĞİŞTİRİCİ SPESİFİKASYON FÖYÜ").font = Font(name="Arial", size=12, bold=True)
        ws.merge_cells(start_row=row, start_column=1, end_row=row, end_column=ncol)
        ws.cell(row, 1).alignment = Alignment(horizontal="center")
        row += 2
        for blk in page:
            cols = col_of[tuple(blk.widths)]
            spans = {}
            for i, a, b in blk.spans:
                spans.setdefault(i, []).append((a, b))
            for i, cells in enumerate(blk.rows):
                is_img = any(k == "img" for _, k in cells)
                covered = set()
                for a, b in spans.get(i, []):
                    covered |= set(range(a + 1, b + 1))
                for j, (t, k) in enumerate(cells):
                    if j in covered:
                        continue
                    end_j = next((b for a, b in spans.get(i, []) if a == j), j)
                    c1, c2 = cols[j], cols[end_j + 1] - 1
                    cell = ws.cell(row, c1, None if k == "img" else (t if t != "" else None))
                    cell.font = fonts[k]
                    cell.alignment = Alignment(horizontal="center" if k in ("c", "bc", "h", "n") else "left",
                                               vertical="center", wrap_text=True)
                    if c2 > c1:
                        ws.merge_cells(start_row=row, start_column=c1, end_row=row, end_column=c2)
                    for cc in range(c1, c2 + 1):
                        ws.cell(row, cc).border = border
                        if k == "h":
                            ws.cell(row, cc).fill = dark
                        elif i in blk.light and j > 0:
                            ws.cell(row, cc).fill = light
                if is_img:
                    ws.row_dimensions[row].height = 220
                    im = XLImage(io.BytesIO(sketch_png(d, "tr", dpi=110)))
                    im.width, im.height = 640, 275
                    ws.add_image(im, f"{get_column_letter(cols[1])}{row}")
                else:
                    ws.row_dimensions[row].height = 12
                row += 1
        if pi == 0:
            breaks.append(row - 1)
        row += 1
    from openpyxl.worksheet.pagebreak import Break
    for b in breaks:
        ws.row_breaks.append(Break(id=b))
    ws.page_setup.paperSize = ws.PAPERSIZE_A4
    ws.page_setup.orientation = "portrait"
    ws.page_setup.scale = 85                       # her föy sayfası bir A4 sayfasına sığar
    ws.print_options.horizontalCentered = True
    ws.page_margins.left = ws.page_margins.right = 0.4
    ws.page_margins.top = ws.page_margins.bottom = 0.4

    sheets = {"Hesaplanan veriler": data_table(d, R)}
    sheets.update({k: v for k, v in (extra_sheets or {}).items() if v is not None and len(v)})
    for name, df in sheets.items():
        w2 = wb.create_sheet(name[:31])
        df = df.reset_index(drop=True)
        for j, col in enumerate(df.columns, 1):
            cc = w2.cell(1, j, str(col))
            cc.font = Font(name="Arial", size=9, bold=True)
            cc.fill = PatternFill("solid", fgColor="DDE7F3")
            w2.column_dimensions[get_column_letter(j)].width = max(12, min(48, len(str(col)) + 4))
        for i, rec in enumerate(df.itertuples(index=False), 2):
            for j, v in enumerate(rec, 1):
                if isinstance(v, (np.floating, float)):
                    v = None if not np.isfinite(v) else round(float(v), 4)
                elif isinstance(v, np.integer):
                    v = int(v)
                elif not isinstance(v, (int, str)) and v is not None:
                    v = str(v)
                w2.cell(i, j, v).font = Font(name="Arial", size=9)
        w2.freeze_panes = "A2"
        if name == "Hesaplanan veriler":
            for col, wd in zip("ABCDE", (11, 52, 22, 12, 10)):
                w2.column_dimensions[col].width = wd
    buf = io.BytesIO()
    wb.save(buf)
    return buf.getvalue()


# =====================================================================
# 7) HTML önizleme
# =====================================================================
def to_html(d):
    from html import escape
    css = ("<style>.tf{border-collapse:collapse;width:100%;table-layout:fixed;font-family:Arial,sans-serif;font-size:11.5px;margin:0}"
           ".tf td{border:1px solid #9e9e99;padding:2px 4px;line-height:1.25}.tf .n{color:#6b6b66;text-align:center;font-size:10px}"
           ".tf .c,.tf .bc{text-align:center}.tf .b,.tf .bc{font-weight:600}.tf .h{background:#2b2b29;color:#fff;text-align:center;"
           "font-weight:700;letter-spacing:.03em}.tf .lt{background:#eeeeea}</style>")
    import base64
    img = base64.b64encode(sketch_png(d, "tr", dpi=130)).decode()
    out = [css]
    for pi, page in enumerate(layout(d)):
        out.append(f"<div style='font:700 13px Arial;text-align:center;margin:{6 if pi == 0 else 18}px 0 4px'>"
                   f"TEMA ISI DEĞİŞTİRİCİ SPESİFİKASYON FÖYÜ — sayfa {pi + 1}/2</div>")
        for blk in page:
            tot = sum(blk.widths)
            h = ["<table class='tf'><colgroup>" + "".join(f"<col style='width:{100 * w / tot:.2f}%'>" for w in blk.widths) + "</colgroup>"]
            spans = {}
            for i, a, b in blk.spans:
                spans.setdefault(i, []).append((a, b))
            for i, cells in enumerate(blk.rows):
                covered = set()
                for a, b in spans.get(i, []):
                    covered |= set(range(a + 1, b + 1))
                tr = "<tr>"
                for j, (t, k) in enumerate(cells):
                    if j in covered:
                        continue
                    span = next((b - a + 1 for a, b in spans.get(i, []) if a == j), 1)
                    cls = k + (" lt" if i in blk.light and j > 0 else "")
                    body = (f"<img src='data:image/png;base64,{img}' style='width:100%'>" if k == "img" else escape(str(t)))
                    tr += f"<td class='{cls}' colspan='{span}'>{body}</td>"
                h.append(tr + "</tr>")
            out.append("".join(h) + "</table>")
    return "".join(out)
