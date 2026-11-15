"""
Kabuk–Tüp Eşanjör Tasarım Uygulaması (kısmi yoğuşturucu)
Çalıştırma:  streamlit run app.py
"""
import io
import json
import math
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import streamlit as st

import esanjor_core as ec
import grafikler as gr
import tema_foy as tf
import gorunum as gv

GA = gv.yukle()    # görünüm ayarları (gorunum.json) — uygulamadaki "Görünümü düzenle" panelinden değişir
st.set_page_config(page_title=GA["pencere_adi"], page_icon=gv.gorsel("simge") or "🔥", layout="wide")

BLUE, ORANGE, GREEN, GRAY, INK = "#2a78d6", "#eb6834", "#1baf7a", "#8a8a85", "#3b3b38"
DEFAULT_COMPONENTS = pd.DataFrame({"Bileşen (thermo adı)": ["acetone", "isopropanol", "hydrogen", "water"],
                                   "Debi (kg/h)": [11400.0, 1300.0, 395.0, 1995.0]})
TUBE_SIZES = {15.88: 1.65, 19.05: 2.11, 25.4: 1.65, 31.75: 2.11, 38.1: 2.11}     # Do mm -> tipik et mm


# =====================================================================
# Önbellekli ağır hesaplar (flash)
# =====================================================================
@st.cache_resource(show_spinner=False)
def get_mixture(components, flows, liquid_model, nrtl_items):
    return ec.Mixture(list(components), list(flows), liquid_model, dict(nrtl_items))


@st.cache_data(show_spinner=False)
def get_curve(components, flows, liquid_model, nrtl_items, T_in, T_out, P_in, P_out, n_zones):
    mix = get_mixture(components, flows, liquid_model, nrtl_items)
    return ec.condensation_curve(mix, T_in, T_out, P_in, P_out, n_zones)


@st.cache_data(show_spinner=False)
def get_ends(components, flows, liquid_model, nrtl_items, T_in, T_out, P_in, P_out):
    """Föy için giriş/çıkış uç verileri: buhar mol kütleleri ve yoğuşmayan gaz debileri."""
    mix = get_mixture(components, flows, liquid_model, nrtl_items)
    nc = mix.noncond
    zs = np.array(mix.zs)

    def info(T, P):
        s_ = mix.flash(T, P)
        if s_.gas is None:
            return np.nan, 0.0
        mV = s_.VF * mix.n_tot * s_.gas.MW() * 3.6
        return s_.gas.MW(), float(mV * np.array(s_.gas.ws())[nc].sum())
    MW_in, nc_in = info(T_in, P_in)
    MW_out, nc_out = info(T_out, P_out)
    MW_nc = float((mix.MW[nc] * zs[nc]).sum() / zs[nc].sum()) if nc.any() else np.nan
    return dict(MW_V_in=MW_in, MW_V_out=MW_out, nc_in=nc_in, nc_out=nc_out, MW_nc=MW_nc)


@st.cache_data(show_spinner=False)
def get_outlet_vapor(components, flows, liquid_model, nrtl_items, T_out, P_out):
    """Çıkıştaki buhar fazında bileşen debileri [kg/h] (gazla kaçan ürün; DWSIM karşılaştırmasında da kullanılır)."""
    return get_mixture(components, flows, liquid_model, nrtl_items).vapor_flows(T_out, P_out)


def parse_nrtl(text):
    """Her satır: bileşen1, bileşen2, b12 [K], b21 [K], alfa"""
    items = []
    for line in text.strip().splitlines():
        if line.strip() and not line.strip().startswith("#"):
            p = [x.strip() for x in line.split(",")]
            items.append(((p[0], p[1]), (float(p[2]), float(p[3]), float(p[4]))))
    return tuple(items)


def fmt(x, d=1):
    return f"{x:,.{d}f}".replace(",", " ")


def draw_tubesheet(Ds, D_otl, Do, Pt, layout, passes, N_used=None):
    """Ayna çizimi: kabuk, dış tüp sınırı, geçiş şeritleri ve tüpler."""
    from matplotlib.collections import PatchCollection
    from matplotlib.patches import Circle
    n, pts, lanes = ec.tube_layout(D_otl, Do, Pt, layout, passes)
    order = np.argsort(np.hypot(pts[:, 0], pts[:, 1]))            # merkezden dışa
    used = np.zeros(n, bool)
    used[order[:N_used if N_used else n]] = True
    fig, ax = plt.subplots(figsize=(5.2, 5.2))
    ax.add_collection(PatchCollection([Circle(p, Do / 2) for p in pts[used]], facecolor=BLUE, edgecolor="white", lw=0.3))
    if (~used).any():
        ax.add_collection(PatchCollection([Circle(p, Do / 2) for p in pts[~used]], facecolor="#d9d9d6", edgecolor="white", lw=0.3))
    ax.add_patch(Circle((0, 0), Ds / 2, fill=False, color="#333", lw=1.5))
    ax.add_patch(Circle((0, 0), D_otl / 2, fill=False, color=GRAY, lw=1, ls="--"))
    for y in lanes["h"]:
        half = math.sqrt(max((Ds / 2) ** 2 - y**2, 0))
        ax.plot([-half, half], [y, y], color=ORANGE, lw=1.2)
    for x in lanes["v"]:
        ax.plot([x, x], [-Ds / 2, Ds / 2], color=ORANGE, lw=1.2)
    ax.set(xlim=(-Ds / 2 * 1.05, Ds / 2 * 1.05), ylim=(-Ds / 2 * 1.05, Ds / 2 * 1.05), aspect="equal")
    ax.axis("off")
    ax.set_title(f"Ds {Ds:.0f} mm · OTL {D_otl:.0f} mm · {passes} geçiş", fontsize=10, loc="left")
    return fig, n


# =====================================================================
# KENAR ÇUBUĞU — GİRDİLER
# =====================================================================
ss = st.session_state
if "pending" in ss:                       # otomatik tasarımın seçtiği değerler (widget'lar oluşmadan önce)
    ss.update(ss.pop("pending"))
for k, v in dict(N=736, passes=4, Do=25.4, t=1.65, layout="triangular", Pt=31.75, L=6000, B=750, Ds=1030).items():
    ss.setdefault(f"geo_{k}", v)


def K(key, default):
    """Anahtarlı girdi: varsayılanı bir kez yazar, anahtarı döndürür (otomatik tasarım bu anahtarlara yazar)."""
    ss.setdefault(key, default)
    return key

_g = gv.gorsel("girdiler")
if _g:
    _w = GA["girdiler_gorseli_genislik"]
    st.sidebar.image(_g, width="stretch" if _w >= 100 else int(3.0 * _w))
st.sidebar.title(gv.temizle("🔧 " + GA["girdiler_basligi"], GA))
duzenle_acik = st.sidebar.toggle(gv.temizle("🎨 Görünümü düzenle", GA),
                                 help="Başlık, sekmeler, görseller ve renkleri koda dokunmadan değiştir.")
gv.yaz(st.sidebar.caption, "m01", "Değerleri girin ve **Hesapla**'ya basın. Varsayılanlar örnektir.")

with st.sidebar.expander("0 · Proje bilgileri (teknik föy için)"):
    proj = dict(project=st.text_input("Müşteri / proje", key=K("proj_project", "Aseton üretim tesisi (IPA dehidrojenasyonu)")),
                plant=st.text_input("Tesis / konum", key=K("proj_plant", "")),
                service=st.text_input("Hizmet", key=K("proj_service", "Reaktör çıkışı kısmi yoğuşturucu")),
                item=st.text_input("Kalem no", key=K("proj_item", "E-102")),
                author=st.text_input("Hazırlayan", key=K("proj_author", "")),
                date=st.text_input("Tarih", key=K("proj_date", pd.Timestamp.today().strftime("%d.%m.%Y"))),
                rev=st.text_input("Revizyon", key=K("proj_rev", "0")))

with st.sidebar.expander("1 · Proses akışı (kabuk tarafı)", expanded=True):
    curve_source = st.radio("Yoğuşma eğrisi", ["Bileşimden hesapla (thermo)", "Tablodan yükle (CSV)"],
                            help="Simülatörden (DWSIM/Aspen) eğriniz varsa tablo seçeneği daha güvenilirdir.")
    comp_mode = st.radio("Bileşen debileri", ["Beslemeden hesapla (IPA–su)", "Elle gir"], key=K("comp_mode", "Beslemeden hesapla (IPA–su)"),
                         horizontal=True, help="Beslemeden: reaktör kütle denkliğiyle E-102 giriş bileşimi bulunur.")
    if comp_mode.startswith("Besleme"):
        c1, c2 = st.columns(2)
        feed_tpy = c1.number_input("Besleme (t/yıl)", min_value=1.0, step=1000.0, key=K("feed_tpy", 110000.0))
        feed_basis = c2.selectbox("Besleme tanımı", ["mixture", "ipa"], key=K("feed_basis", "mixture"),
                                  format_func={"mixture": "IPA–su karışımı (toplam)", "ipa": "Yalnız IPA"}.get)
        feed_w = c1.number_input("Beslemede IPA (% kütle)", min_value=50.0, max_value=100.0, key=K("feed_w", 87.7),
                                 help="IPA–su azeotropu: kütlece %87,7 IPA, 80,4 °C (1 atm).")
        feed_X = c2.number_input("Tek geçiş dönüşümü (%)", min_value=1.0, max_value=100.0, key=K("feed_X", 90.0))
        feed_h = c1.number_input("Çalışma (saat/yıl) ", min_value=1000.0, max_value=8760.0, key=K("feed_h", 8000.0))
        eff = ec.reactor_effluent(feed_tpy, feed_h, feed_w / 100, feed_X / 100, feed_basis)
        comp_df = eff[["Bileşen (thermo adı)", "Debi (kg/h)"]]
        st.dataframe(eff.round(1), hide_index=True, width="stretch")
        gv.yaz(st.caption, "m137", "Reaktör çıkışı = E-102 girişi: {toplam} kg/h · aseton üretimi ≈ {uretim} t/yıl",
               degerler=dict(toplam=fmt(eff["Debi (kg/h)"].sum(), 0), uretim=fmt(eff["Debi (kg/h)"].iloc[0] * feed_h / 1000, 0)))
    else:
        comp_df = st.data_editor(DEFAULT_COMPONENTS, num_rows="dynamic", width="stretch", key="comp")
    uploaded = None
    if curve_source.startswith("Tablo"):
        uploaded = st.file_uploader("Eğri tablosu (CSV)", type="csv",
                                    help="Kolonlar: " + ", ".join(ec.TABLE_COLUMNS) +
                                         ". yNC, Mc, DAB_m2s yoksa Colburn–Hougen yerine SBG kullanılır.")
    c1, c2 = st.columns(2)
    T_hot_in = c1.number_input("Giriş T (°C)", key=K("T_hot_in", 120.0))
    T_hot_out = c2.number_input("Çıkış T (°C)", key=K("T_hot_out", 40.0))
    P_in_bar = c1.number_input("Giriş P (bar, mutlak)", min_value=0.1, key=K("P_in_bar", 2.8))
    P_out_bar = c2.number_input("Çıkış P (bar, mutlak)", min_value=0.1, key=K("P_out_bar", 2.6))
    Rf_shell = st.number_input("Kabuk kirlenme direnci (m²K/W)", format="%.5f", key=K("Rf_shell", 0.00023))
    liquid_model = st.selectbox("Sıvı fazı modeli", ["nrtl", "ideal"], key=K("liquid_model", "nrtl"),
                                format_func={"nrtl": "NRTL (su–aseton–IPA seti yerleşik)", "ideal": "İdeal (karşılaştırma)"}.get,
                                help="Yerleşik NRTL seti: NIST ölçümleri + literatür (OSTI 799235), α = 0,3. "
                                     "IPA–su azeotropunu %87,5 kütle / 80,0 °C olarak verir.")
    nrtl_text = ""
    if liquid_model == "nrtl":
        nrtl_text = st.text_area("Ek / farklı NRTL parametreleri (isteğe bağlı)", "# bileşen1, bileşen2, b12 [K], b21 [K], alfa\n",
                                 help="τ = b/T. Girilen ikili yerleşik değerin yerine geçer.")
    n_zones = st.slider("Bölge sayısı", 6, 40, 15)
    names = [c for c in comp_df.iloc[:, 0].dropna().astype(str)]
    product = "acetone" if "acetone" in names else (names[0] if names else None)

with st.sidebar.expander("2 · Soğutma suyu (tüp tarafı)"):
    c1, c2 = st.columns(2)
    w_T_in = c1.number_input("Giriş T (°C) ", key=K("w_T_in", 25.0))
    ss.setdefault("w_T_out", 35.0)
    w_T_out = c2.number_input("Çıkış T (°C) ", key="w_T_out")
    w_m = st.number_input("Debi (kg/h) — 0: ısı dengesinden", min_value=0.0, key=K("w_m", 0.0))
    c1, c2 = st.columns(2)
    w_rho = c1.number_input("ρ (kg/m³)", value=995.0)
    w_cp = c2.number_input("Cp (kJ/kg·K)", value=4.18)
    w_mu = c1.number_input("μ (cP)", value=0.80, format="%.3f")
    w_k = c2.number_input("k (W/m·K)", value=0.615, format="%.3f")
    w_Rf = st.number_input("Su kirlenme direnci (m²K/W)", format="%.5f", key=K("w_Rf", 0.0002))
    w_P = st.number_input("Su giriş basıncı (bar, mutlak)", key=K("w_P", 4.0))

with st.sidebar.expander("3 · Kabuk h yöntemi"):
    h_method = st.selectbox("Yöntem", ["colburn_hougen", "sbg", "fixed"], key=K("h_method", "colburn_hougen"),
                            format_func={"colburn_hougen": "Colburn–Hougen (önerilen)", "sbg": "SBG (basit)",
                                         "fixed": "Sabit h (elle)"}.get)
    h_fixed = st.number_input("Sabit h (W/m²K) — karşılaştırmada da kullanılır", value=450.0)
    h_mult = st.number_input("Hesaplanan h çarpanı (emniyet)", step=0.05, key=K("h_mult", 1.0))
    shell_method = st.selectbox("Kabuk tarafı tek faz yöntemi", ["bell_delaware", "kern"], key=K("shell_method", "bell_delaware"),
                                format_func={"bell_delaware": "Bell–Delaware (kaçak + bypass, önerilen)",
                                             "kern": "Kern (karşılaştırma)"}.get,
                                help="Gaz fazı h'si ve kabuk ΔP'nin tek fazlı referansı.")
    n_ss_auto = st.checkbox("Sızdırmazlık şeridi sayısı otomatik", key=K("n_ss_auto", True),
                            help="Demet–kabuk boşluğu > 20 mm ise her 6 çapraz akış sırasına bir çift (TEMA RCB-4.8).")
    n_ss = None if n_ss_auto else st.number_input("Sızdırmazlık şeridi çifti", min_value=0, max_value=20, key=K("n_ss", 0))

with st.sidebar.expander("4 · Geometri"):
    c1, c2 = st.columns(2)
    c1.selectbox("Tüp dış çapı Do (mm)", list(TUBE_SIZES), key="geo_Do")
    c2.number_input("Et kalınlığı (mm)", key="geo_t", step=0.1)
    c1.number_input("Tüp sayısı N", key="geo_N", step=4, min_value=10)
    c2.selectbox("Tüp geçiş sayısı", [1, 2, 4, 6, 8], key="geo_passes")
    c1.selectbox("Dizilim", ["triangular", "square", "rotated square"], key="geo_layout")
    c2.number_input("Tüp adımı Pt (mm)", key="geo_Pt", step=0.5)
    c1.number_input("Tüp boyu L (mm)", key="geo_L", step=100)
    c2.number_input("Baffle aralığı B (mm)", key="geo_B", step=10)
    st.number_input("Kabuk iç çapı Ds (mm)", key="geo_Ds", step=10)
    head_type = st.selectbox("Arka kafa tipi (demet boşluğu için)", list(ec.BUNDLE_CLEARANCE), key=K("head_type", "Sabit ayna (BEM)"),
                             on_change=lambda: ss.update(clearance=ec.BUNDLE_CLEARANCE[ss["head_type"]]))
    clearance = st.number_input("Demet–kabuk çap boşluğu (mm)", key=K("clearance", 15.0),
                                help="Towler & Sinnott Şekil 19.12'den okuyun.")
    c1, c2 = st.columns(2)
    shells = c1.number_input("Seri kabuk sayısı", min_value=1, max_value=4, key=K("shells", 1))
    k_wall = c2.number_input("Duvar k (W/m·K)", value=60.0)
    K_pass = st.number_input("Tüp dönüş kaybı (hız yükü/geçiş)", key=K("K_pass", 2.5))

with st.sidebar.expander("5 · Tasarım sınırları"):
    c1, c2 = st.columns(2)
    v_min = c1.number_input("Tüp hızı min (m/s)", key=K("v_min", 1.0))
    v_max = c2.number_input("Tüp hızı max (m/s)", key=K("v_max", 2.5))
    dp_tube_max = c1.number_input("Tüp ΔP max (kPa)", key=K("dp_tube_max", 35.0))
    dp_shell_max = c2.number_input("Kabuk ΔP max (kPa)", key=K("dp_shell_max", 20.0))
    od_min = c1.number_input("Aşırı tasarım min (%)", key=K("od_min", 10.0))
    od_max = c2.number_input("Aşırı tasarım max (%)", key=K("od_max", 30.0))

with st.sidebar.expander("6 · Mekanik"):
    c1, c2 = st.columns(2)
    P_d_shell = c1.number_input("Kabuk tasarım P (barg)", value=3.5)
    P_d_tube = c2.number_input("Tüp tasarım P (barg)", value=5.0)
    S_shell = c1.number_input("Kabuk izin ver. gerilme (MPa)", value=138.0, help="SA-516 Gr.70 ≈ 138 MPa")
    S_tube = c2.number_input("Tüp izin ver. gerilme (MPa)", value=92.0, help="SA-179 ≈ 92 MPa")
    E_joint = c1.number_input("Kaynak verimi E", value=0.85)
    CA_shell = c2.number_input("Korozyon payı (mm)", value=3.0)

with st.sidebar.expander("7 · Maliyet"):
    corr = st.selectbox("Maliyet bağıntısı", ["Otomatik (TEMA tipine göre)"] + list(ec.COST_CORRELATIONS))
    material = st.selectbox("Malzeme", list(ec.MATERIAL_FACTORS))
    cepci_now = st.number_input("Güncel CEPCI", value=800.0,
                                help="2023 yıllık ortalaması 797,9; 2024 ~800 civarı. Dersinizin kullandığı değeri girin.")
    gv.yaz(st.caption, "m02", "Baz: Towler & Sinnott, Ocak 2010, CEPCI = 532,9")
    c1, c2 = st.columns(2)
    i_rate = c1.number_input("Faiz oranı (%)", value=10.0) / 100
    years = c2.number_input("Ömür (yıl)", value=10, min_value=1)
    hours = c1.number_input("Çalışma (saat/yıl)", value=8000)
    pump_eff = c2.number_input("Pompa verimi", value=0.70)
    elec_price = c1.number_input("Elektrik ($/kWh)", value=0.10, format="%.3f")
    water_price = c2.number_input("Soğutma suyu ($/m³)", value=0.02, format="%.3f")
    fx_rate = c1.number_input("Kur (TL/USD), 0: gösterme", value=0.0)
    with st.popover("Kurulum faktörleri"):
        inst = {k: st.number_input(k, value=v, key=f"inst_{k}") for k, v in ec.INSTALL_FACTORS.items()}

def standard_settings():
    """Otomatik tasarımın kullandığı standart ayarlar (TEMA / Towler & Sinnott). Proses girdilerine dokunmaz."""
    dp_proc = max((ss["P_in_bar"] - ss["P_out_bar"]) * 100, 5.0)
    return dict(Rf_shell=0.0002, w_Rf=0.0002, w_m=0.0, h_method="colburn_hougen", h_mult=1.0, K_pass=2.5, shells=1,
                v_min=1.0, v_max=2.5, dp_tube_max=35.0, dp_shell_max=round(dp_proc, 1), od_min=15.0, od_max=30.0,
                head_type="Sabit ayna (BEM)", clearance=15.0)


def start_auto():
    ss.update(standard_settings())
    for k in ("opt_res", "scan", "sizing_res", "auto_report"):
        ss.pop(k, None)
    ss["auto_stage"] = "optimize"


st.sidebar.divider()
st.sidebar.button("🧭 Benim için tasarla", type="primary", width="stretch", on_click=start_auto,
                  help="Yalnız proses verisini (bölüm 1–2) girin. Geri kalan her şeyi uygulama TEMA/Towler aralıklarından "
                       "seçer, tüm standart geometrileri dener ve en uygun sağlam tasarımı girdilere yazar (~1 dk).")
run = st.sidebar.button("▶ Hesapla (kendi girdilerimle)", width="stretch")


# =====================================================================
# HESAP
# =====================================================================
def compute():
    comps = tuple(comp_df.iloc[:, 0].dropna().astype(str))
    flows = tuple(float(x) for x in comp_df.iloc[:, 1].dropna())
    nrtl_items = parse_nrtl(nrtl_text) if liquid_model == "nrtl" else ()
    tube = dict(m=w_m or None, T_in=w_T_in, T_out=w_T_out, rho=w_rho, cp=w_cp, mu=w_mu / 1000, k=w_k, Rf=w_Rf)
    geo = {k: ss[f"geo_{k}"] for k in ("N", "passes", "Do", "t", "layout", "Pt", "L", "B", "Ds")}
    geo["N"] = int(geo["N"])
    # demet çapı: taramayla aynı tanım (kabuk iç çapı − boşluk); tüpler sığmıyorsa en küçük sığdıran sınır
    _Db = float(geo["Ds"]) - float(clearance)
    if ec.tube_layout(_Db, geo["Do"], geo["Pt"], geo["layout"], geo["passes"])[0] >= geo["N"]:
        geo["Db"] = _Db
    else:
        geo["Db"] = ec.otl_for_tubes(geo["N"], geo["Do"], geo["Pt"], geo["layout"], geo["passes"])
    kw = dict(h_method=h_method, h_fixed=h_fixed, h_mult=h_mult, Rf_shell=Rf_shell, k_wall=k_wall,
              shell_method=shell_method, n_ss=n_ss,
              K_pass=K_pass, shells=int(shells))
    use_table = curve_source.startswith("Tablo")

    if use_table:
        if uploaded is None:
            raise ValueError("Tablo seçeneği için CSV dosyası yükleyin.")
        curve, T_dew = ec.curve_from_table(pd.read_csv(uploaded))
        if comps and curve["D_AB"].isna().all():          # yNC varsa D_AB bileşimden tamamlanır
            curve = ec.fill_diffusion(curve, get_mixture(comps, flows, "ideal", ()))
        curve2 = None
    else:
        args = (comps, flows, liquid_model, nrtl_items)
        curve, T_dew = get_curve(*args, T_hot_in, T_hot_out, P_in_bar, P_out_bar, n_zones)
        curve2, _ = get_curve(*args, T_hot_in, T_hot_out, P_in_bar, P_out_bar, 2 * n_zones)

    r = ec.rate_condensing(tube, geo, curve, **kw)
    ends = None if use_table else get_ends(comps, flows, liquid_model, nrtl_items, T_hot_in, T_hot_out, P_in_bar, P_out_bar)
    mech = dict(P_design_shell_barg=P_d_shell, P_design_tube_barg=P_d_tube, S_shell_MPa=S_shell,
                S_tube_MPa=S_tube, E_joint=E_joint, CA_shell_mm=CA_shell, alpha=12e-6)
    m = ec.mechanical(r, curve, T_dew, mech)
    corr_name = corr
    if corr.startswith("Otomatik"):
        corr_name = "Yüzer kafa" if m["tema"].startswith("AES") else "U-tüp (sabit ayna için de yaklaşık)"
    cost = dict(correlation=corr_name, f_m=ec.MATERIAL_FACTORS[material], cepci_now=cepci_now, cepci_base=532.9,
                install=inst, i_rate=i_rate, years=int(years), hours=hours, elec_price=elec_price,
                pump_eff=pump_eff, water_price_m3=water_price, fx_rate=fx_rate or None)
    c = ec.cost_estimate(r, cost)

    # karşılaştırma ve doğrulama
    single = ec.rate_single_lmtd(r, h_fixed)
    methods = []
    for mm in ("fixed", "sbg", "colburn_hougen"):
        rr = r if mm == h_method else ec.rate_condensing(tube, geo, curve, **dict(kw, h_method=mm))
        methods.append((mm, rr))
    r2 = ec.rate_condensing(tube, geo, curve2, **kw) if curve2 is not None else None
    sens = ec.sensitivity(tube, geo, curve, r, **kw)
    hv = ec.hydraulic_verification(tube, geo, curve, r, kw)
    vap = None if use_table else get_outlet_vapor(comps, flows, liquid_model, nrtl_items,
                                                  float(curve["T"].iloc[-1]), float(curve["P"].iloc[-1]))
    prod_V = vap.get(product) if vap and product else None
    return dict(r=r, curve=curve, T_dew=T_dew, mech=m, ends=ends, use_table=use_table, cost=c, cost_in=cost, single=single, methods=methods,
                r2=r2, sens=sens, hv=hv, prod_V=prod_V, tube=tube, geo=geo, kw=kw, comps=comps, flows=flows,
                clearance=clearance, limits=dict(v=(v_min, v_max), dp_t=dp_tube_max, dp_s=dp_shell_max, od=(od_min, od_max)))


def mech_inputs():
    return dict(P_design_shell_barg=P_d_shell, P_design_tube_barg=P_d_tube, S_shell_MPa=S_shell, S_tube_MPa=S_tube,
                E_joint=E_joint, CA_shell_mm=CA_shell, alpha=12e-6)


def cost_fn_from(Rx):
    """Aday geometri -> maliyet (TEMA tipine göre bağıntı otomatik)."""
    def f(rr):
        cc = dict(Rx["cost_in"])
        if corr.startswith("Otomatik"):
            mm = ec.mechanical(rr, Rx["curve"], Rx["T_dew"], mech_inputs())
            cc["correlation"] = "Yüzer kafa" if mm["tema"].startswith("AES") else "U-tüp (sabit ayna için de yaklaşık)"
        return ec.cost_estimate(rr, cc)
    return f


WATER_IN_START, WATER_IN_STEP = 20, 5      # soğutma suyu giriş sıcaklıkları 20 °C'den 5'er artırılarak denenir


def auto_design(Rx, bar):
    """Soğutma suyu girişi (20, 25, 30 … °C) × su çıkışı × tüm standart geometriler denenir.
    Kısıtları sağlayan adaylar üç kademede değerlendirilir: TEMA zorunlu hükümleri, kritik ölçütler (K2, K3, K9, K10)
    ve toplam puan (K1–K11); eşitlikte en düşük yıllık maliyet."""
    curve_, tube_, kw_ = Rx["curve"], Rx["tube"], Rx["kw"]
    lim_ = dict(Rx["limits"], od=(15.0, 35.0))
    T_hot_out = float(curve_["T"].iloc[-1])
    T_ins = list(range(WATER_IN_START, int(T_hot_out) - 5 + 1, WATER_IN_STEP))
    T_ins = [T for T in T_ins if T + 5 < T_hot_out]
    frames, T_lists = [], {}
    for i, Ti in enumerate(T_ins):
        tb = dict(tube_, T_in=float(Ti))
        T_list = list(range(Ti + 5, 45 + 1))
        T_lists[Ti] = T_list
        allc_i, _ = ec.optimize(tb, curve_, T_list, cost_fn_from(Rx), kw_, lim_, clearance=ec.BUNDLE_CLEARANCE["Sabit ayna (BEM)"],
                                progress=lambda f, T, Ti=Ti, i=i: bar.progress((i + f) / len(T_ins) * 0.6,
                                    text=f"2/3 · Su girişi {Ti} °C, çıkışı {T} °C: tüm geometriler deneniyor…"))
        if not allc_i.empty:
            frames.append(allc_i.assign(T_su_giris=float(Ti)))
    if not frames:
        raise ValueError("Denenen su giriş sıcaklıklarının hiçbirinde kısıtları sağlayan tasarım bulunamadı.")
    allc = pd.concat(frames, ignore_index=True)
    ev = ec.technical_evaluation(tube_, curve_, allc, kw_, lim_,
                                 progress=lambda f: bar.progress(0.6 + 0.4 * f, text="2/3 · Adaylar teknik ölçütlerle sınanıyor…"))
    row = ev.iloc[0]
    Ti_sel = float(row.T_su_giris)
    geo_ = dict(N=int(row.N), passes=int(row.gecis), Do=float(row.Do), t=TUBE_SIZES[float(row.Do)], layout="triangular",
                Pt=round(1.25 * float(row.Do), 2), L=int(row.L), B=int(row.B), Ds=int(row.Ds))
    tube_o = dict(tube_, T_in=Ti_sel, T_out=float(row.T_su_cikis), m=None)
    r_ = ec.rate_condensing(tube_o, geo_, curve_, **kw_)
    mech_ = ec.mechanical(r_, curve_, Rx["T_dew"], mech_inputs())
    head = "Yüzer kafa, split-ring (AES)" if mech_["tema"].startswith("AES") else "Sabit ayna (BEM)"
    clr = ec.BUNDLE_CLEARANCE[head]
    while ec.tube_layout(geo_["Ds"] - clr, geo_["Do"], geo_["Pt"], "triangular", geo_["passes"])[0] < geo_["N"]:
        geo_["Ds"] += 10                                         # tüpler kabuğa sığana kadar
    sel = ev[ev.T_su_giris == Ti_sel]
    best = sel.loc[sel.groupby("T_su_cikis")["yillik_toplam_USD"].idxmin()].sort_values("T_su_cikis").reset_index(drop=True)
    sweep = ec.tac_vs_tubes_sweep(tube_o, curve_, geo_["Do"], geo_["t"], geo_["passes"], (3048, 4877, 6096),
                                  cost_fn_from(Rx), kw_, lim_, clearance=clr, N_ref=geo_["N"] * geo_["L"] / 3048)
    pending = {f"geo_{k}": v for k, v in geo_.items()}
    pending.update(w_T_in=Ti_sel, w_T_out=float(row.T_su_cikis), head_type=head, clearance=clr, od_min=15.0, od_max=35.0)
    by_inlet = ev.groupby("T_su_giris").agg(aday=("puan", "size"), en_yuksek_kritik=("kritik", "max"),
                                            en_yuksek_puan=("puan", "max"),
                                            en_ucuz_TAC=("yillik_toplam_USD", "min")).reset_index()
    opt_res = dict(allc=ev, best=best, top=ev.head(5), T_list=T_lists[int(Ti_sel)], h_rob=0.8, r_opt=r_, sweep=sweep)
    return dict(pending=pending, opt_res=opt_res, top=ev.head(5), row=row, rank=1, robust_ok=bool(row.K2), geo=geo_, head=head, tema=mech_,
                n_feasible=len(ev), T_list=T_lists[int(Ti_sel)], T_ins=T_ins, T_in=Ti_sel, limits=lim_, by_inlet=by_inlet)



# =====================================================================
# DWSIM doğrulaması (yalnız DWSIM kurulu Windows bilgisayarda)
# =====================================================================
import platform as _platform
import subprocess as _subprocess
import sys as _sys
import tempfile as _tempfile

import dwsim_dogrulama as dwv

_DW_AD = {"acetone": "Acetone", "isopropanol": "Isopropanol", "2-propanol": "Isopropanol",
          "hydrogen": "Hydrogen", "water": "Water"}
_CAS_AD = {"7732-18-5": "Water", "67-64-1": "Acetone", "67-63-0": "Isopropanol"}


def dwsim_klasoru_bul():
    return dwv.dwsim_klasoru() if _platform.system() == "Windows" else None


def _dw_geometri(g_, Rf_s, Rf_t, k_w):
    return {"Ds_mm": float(g_["Ds"]), "B_mm": float(g_["B"]), "baffle_kesim_yuzde": 25.0,
            "Do_mm": float(g_["Do"]), "Di_mm": round(float(g_["Do"]) - 2 * float(g_["t"]), 3),
            "L_m": float(g_["L"]) / 1000, "N": int(g_["N"]), "gecis": int(g_["passes"]),
            "Pt_mm": float(g_["Pt"]), "dizilim": "ucgen" if str(g_.get("layout", "tri")).startswith("tri") else "kare",
            "kirlenme_kabuk": float(Rf_s), "kirlenme_tup": float(Rf_t), "k_duvar": float(k_w), "puruzluluk_mm": 0.045}


def aday_geo(row):
    Do_ = float(row.Do)
    return dict(N=int(row.N), passes=int(row.gecis), Do=Do_, t=TUBE_SIZES[Do_], layout="triangular",
                Pt=round(1.25 * Do_, 2), L=int(row.L), B=int(row.B), Ds=int(row.Ds))


def dwsim_girdisi(R_, adaylar=None):
    """DWSIM betiğinin girdisi: bileşenler, akımlar, seçili geometri, NRTL ve uygulama sonuçları.
    adaylar verilirse (DataFrame) her biri ayrı bir eşanjör testi olarak eklenir."""
    r_, cv_, w_ = R_["r"], R_["curve"], R_["r"]["water"]
    g = {
        "bilesenler": {_DW_AD.get(c.lower(), c.title()): round(float(fl), 2) for c, fl in zip(R_["comps"], R_["flows"])},
        "proses": {"T_giris_C": float(cv_["T"].iloc[0]), "T_cikis_C": float(cv_["T"].iloc[-1]),
                   "P_giris_bar": float(cv_["P"].iloc[0]), "P_cikis_bar": float(cv_["P"].iloc[-1])},
        "su": {"T_giris_C": float(w_["T_in"]), "P_bar": float(w_P), "debi_kgh": round(float(w_["m"]), 1)},
        "geometri": _dw_geometri(r_["geo"], r_["Rf_shell"], w_["Rf"], r_["k_wall"]),
        "nrtl": [[_CAS_AD[i_], _CAS_AD[j_], round(bij * 1.98721, 1), round(bji * 1.98721, 1), alf]
                 for (i_, j_), (bij, bji, alf) in ec.NRTL_DEFAULT.items()],
        "uygulama": {"Q_kW": round(float(r_["Q_kW"]), 1), "T_ciy_C": round(float(R_["T_dew"]), 2),
                     "gaz_cikis_kgh": round(float(cv_.m_V.iloc[-1]), 1), "sivi_cikis_kgh": round(float(cv_.m_L.iloc[-1]), 1),
                     "su_cikis_C": round(float(w_["T_out"]), 2), "U_kirli": round(float(tf.zonal_U(r_)["U_dirty"]), 1),
                     "A_m2": round(float(r_["A"]), 2), "A_gerekli_m2": round(float(r_["A_req"]), 2),
                     "dP_tup_kPa": round(float(r_["dp_tube_kPa"]), 2), "dP_kabuk_kPa": round(float(r_["dp_shell_kPa"]), 2),
                     "aseton_gazda_kgh": None if R_.get("prod_V") is None else round(float(R_["prod_V"]), 1)},
        "kabuk_dP_siniri_kPa": float(R_["limits"]["dp_s"]),
    }
    if adaylar is not None and len(adaylar):
        g["tasarimlar"] = [{"ad": f"{i + 1}) {int(x.N)}/{int(x.Ds)}",
                            "su": {"T_giris_C": float(x.T_su_giris), "debi_kgh": round(float(x.su_debisi_kgh), 1)},
                            "geometri": _dw_geometri(aday_geo(x), r_["Rf_shell"], w_["Rf"], r_["k_wall"]),
                            "uygulama": {"Q_kW": round(float(r_["Q_kW"]), 1), "su_cikis_C": float(x.T_su_cikis),
                                         "dP_tup_kPa": round(float(x.dP_tup), 2), "dP_kabuk_kPa": round(float(x.dP_kabuk), 2)}}
                           for i, x in enumerate(adaylar.itertuples())]
    return g


def dwsim_calistir(girdi, timeout=900):
    """Betiği ayrı bir işlem olarak çalıştırır (DWSIM hatası uygulamayı çökertmez). Döndürür: satırlar, kayıt."""
    betik = os.path.join(os.path.dirname(os.path.abspath(__file__)), "dwsim_dogrulama.py")
    gizli = getattr(_subprocess, "CREATE_NO_WINDOW", 0)
    if _subprocess.run([_sys.executable, "-c", "import pythonnet"], capture_output=True, creationflags=gizli).returncode:
        _subprocess.run([_sys.executable, "-m", "pip", "install", "--quiet", "pythonnet", "pywin32"],
                        capture_output=True, creationflags=gizli)
    with _tempfile.TemporaryDirectory() as td:
        gfile, cfile = os.path.join(td, "girdi.json"), os.path.join(td, "cikti.json")
        with open(gfile, "w", encoding="utf-8") as f_:
            json.dump(girdi, f_, ensure_ascii=False)
        try:
            ortam = dict(os.environ, PYTHONIOENCODING="utf-8", PYTHONUTF8="1")
            p_ = _subprocess.run([_sys.executable, betik, gfile, f"--cikti={cfile}"], capture_output=True, text=True,
                                 encoding="utf-8", errors="replace", timeout=timeout, creationflags=gizli,
                                 cwd=os.path.dirname(betik), env=ortam)
            sat = json.load(open(cfile, encoding="utf-8")) if os.path.exists(cfile) else None
            return dict(satirlar=sat, log=(p_.stdout or "") + (p_.stderr or ""))
        except _subprocess.TimeoutExpired:
            return dict(satirlar=None, log="DWSIM süre sınırı içinde yanıt vermedi.")


def dwsim_termo_uyum(sat):
    v_ = [x["uyum"] for x in sat if x["test"] == "Termodinamik" and x.get("uyum") is not None
          and x["buyukluk"] in ("Isı yükü", "Çiy noktası", "Çıkış gaz debisi", "Çıkış sıvı debisi")]
    return sum(v_) / len(v_) if v_ else None


def _yuv(v_, n):
    return None if v_ is None or v_ != v_ else round(float(v_), n)


def dwsim_adaylari_degerlendir(sat, adaylar, R_):
    """Her aday için DWSIM sonucu: ısı yükü karşılanıyor mu, basınç düşüşleri sınırda mı."""
    Q_req, lim_ = float(R_["r"]["Q_kW"]), R_["limits"]
    rows = []
    for i, x in enumerate(adaylar.itertuples()):
        ad = f"{i + 1}) {int(x.N)}/{int(x.Ds)}"[:17]

        def v(b):
            return next((s["dwsim"] for s in sat if s["test"] == ad and s["buyukluk"] == b), float("nan"))
        Q_, Tout, U_, dpt, dps = v("Isı yükü"), v("Proses çıkış sıcaklığı"), v("U (kirli)"), v("Tüp tarafı ΔP"), v("Kabuk ΔP")
        neden = []
        if Q_ != Q_:
            neden.append("DWSIM sonucu yok")
        elif not Q_ >= 0.98 * Q_req:
            neden.append(f"ısı yükü yetersiz ({Q_:.0f} kW)")
        if not dps <= lim_["dp_s"]:
            neden.append(f"kabuk ΔP {dps:.1f} > {lim_['dp_s']:.0f} kPa")
        if not dpt <= lim_["dp_t"]:
            neden.append(f"tüp ΔP {dpt:.1f} > {lim_['dp_t']:.0f} kPa")
        rows.append(dict(Sıra=i + 1, Tasarım=f"{int(x.N)} tüp, Ø{x.Do} × {x.L / 1000:.2f} m, {int(x.gecis)} geçiş, Ds {int(x.Ds)}",
                         Kritik=f"{int(x.kritik)}/{len(ec.CRITICAL_CRITERIA)}" if hasattr(x, "kritik") else "–",
                         Puan=f"{int(x.puan)}/{len(ec.TECH_CRITERIA)}" if hasattr(x, "puan") else "–",
                         **{"Yıllık maliyet ($)": round(float(x.yillik_toplam_USD))},
                         **{"DWSIM ısı yükü (kW)": _yuv(Q_, 0), "DWSIM proses çıkışı (°C)": _yuv(Tout, 1),
                            "Kabuk ΔP: uyg. / DWSIM": f"{x.dP_kabuk:.1f} / {dps:.1f}",
                            "Tüp ΔP: uyg. / DWSIM": f"{x.dP_tup:.1f} / {dpt:.1f}", "DWSIM U": _yuv(U_, 0)},
                         Sonuç="✅ GEÇTİ" if not neden else "❌ " + "; ".join(neden), gecti=not neden))
    return pd.DataFrame(rows)


def dwsim_adaylari(n=5):
    """DWSIM'de denenecek tasarımlar: otomatik tasarımın teknik sıralaması ya da tarama sonucu (en düşük TAC)."""
    O_ = ss.get("opt_res")
    if not O_ or O_["allc"].empty:
        return None
    ev_ = O_["allc"]
    if "kritik" not in ev_:
        ev_ = ev_.sort_values("yillik_toplam_USD")
    if "T_su_giris" not in ev_:
        ev_ = ev_.assign(T_su_giris=float(ss["res"]["tube"]["T_in"]))
    return ev_.drop_duplicates(["N", "Ds", "L", "gecis", "B", "Do"]).head(n).reset_index(drop=True)


stage = ss.get("auto_stage")
if stage == "optimize":
    ss["auto_stage"] = None
    ok = False
    try:
        with st.spinner("1/3 · Yoğuşma eğrisi hesaplanıyor…"):
            R0 = compute()
        bar = st.progress(0.0, text="2/3 · Tasarımlar deneniyor…")
        rep_ = auto_design(R0, bar)
        bar.empty()
        ss["pending"] = rep_.pop("pending")
        ss["opt_res"] = rep_.pop("opt_res")
        ss["auto_report"] = rep_
        ss["auto_stage"] = "final"
        ok = True
    except Exception as e:
        st.error(f"Otomatik tasarım yapılamadı: {e}")
    if ok:
        st.rerun()
elif stage == "final":
    with st.spinner("3/3 · Seçilen tasarım ayrıntılı hesaplanıyor…"):
        ss["res"] = compute()
    if "opt_res" in ss:
        ss["opt_res"]["fp"] = None                    # aşağıda opt_fingerprint ile işaretlenir
        ss["opt_fp_pending"] = True
    ss["auto_stage"] = None
    ss["auto_done"] = True

if run:
    try:
        with st.spinner("Hesaplanıyor (ilk çalıştırmada flash hesabı birkaç saniye sürer)…"):
            ss["res"] = compute()
            ss.pop("scan", None)
    except Exception as e:
        st.error(f"Hesap yapılamadı: {e}")

_g = gv.gorsel("baslik")
if _g:
    st.image(_g, width=GA["baslik_gorseli_genislik"])
st.title(GA["baslik"])
if GA["alt_yazi"].strip():
    st.markdown(GA["alt_yazi"])
if duzenle_acik:
    with st.container(border=True):
        gv.duzenleme_paneli(GA)
def tanitim_goster():
    """İlk sayfa: uygulamanın amacı, kullanımı, yöntemi, seçim kuralı, doğrulama, varsayımlar ve sınırlamalar."""
    gv.yaz(st.markdown, "m80", """
### E-102 Kısmi Yoğuşturucu Tasarım Uygulaması

Bu uygulama, izopropil alkolden (IPA) aseton üretilen tesiste reaktör çıkış gazını soğutan **E-102 kabuk-tüp
eşanjörünü** tasarlar. Proses gazı (aseton, IPA, su buharı ve hidrojen) kabuk tarafında kısmen yoğuşur, soğutma suyu
tüplerden akar. Hidrojen bu koşullarda yoğuşmadığı için eşanjör bir **kısmi yoğuşturucudur**. Hesaplar ve konstrüksiyon
kararları **TEMA Standartları (9. baskı, Sınıf B)** esas alınarak yapılır.
""")
    gv.yaz(st.markdown, "m81", """
#### Nasıl kullanılır?
1. Sol menüde proses ve soğutma suyu verilerini kontrol edin; varsayılan değerler proje verisidir.
2. **Benim için tasarla**'ya basın. Uygulama yüzlerce tasarımı dener ve teknik olarak en uygun olanı seçer.
3. Sonuçları sekmelerden inceleyin: **Hesap adımları** her denklemi sayılarıyla, **Bölgeler & doğrulama** bağımsız
   kontrolleri, **Tasarım taraması** alternatifleri, **DWSIM karşılaştırma** simülatörle kıyası gösterir.
4. TEMA föyünü (PDF), tüm hesapları (Excel) ve grafikleri **TEMA föyü** sekmesinden indirin.

Kendi tasarımınızı denemek için sol menüdeki girdileri değiştirip **Hesapla**'ya basabilirsiniz.
""")
    c1_, c2_ = st.columns(2)
    with c1_:
        gv.yaz(st.markdown, "m82", """
#### Hesap yöntemi
| Adım | Yöntem |
|---|---|
| Reaktör çıkışı | Kütle denkliği |
| Faz dengesi | Peng–Robinson (gaz), NRTL (sıvı) |
| Tüp tarafı | Gnielinski (h), Colebrook (ΔP), geçiş başına 2,5 hız yükü |
| Kabuk tarafı h | Bell–Delaware (gaz), Nusselt (yoğuşma filmi), Colburn–Hougen (hidrojen etkisi) |
| Gerekli alan | Bölgesel LMTD ve F, ε-NTU ile kontrol |
| Kabuk ΔP | Bell–Delaware, bölgesel iki faz çarpanı (C = 8), nozullar |
| Tüp sayısı | Tüplerin aynaya geometrik yerleşimi (geçiş şeritleri dahil) |
| Mekanik | TEMA 9. baskı (ısıl genleşme, ρv², destek açıklığı, titreşim, ayna), ASME VIII Div. 1 |
| Maliyet | Towler & Sinnott, CEPCI |
""")
    with c2_:
        gv.yaz(st.markdown, "m83", """
#### Tasarım nasıl seçilir?
Soğutma suyu girişi 20, 25 ve 30 °C için su çıkışı ve standart geometriler (tüp çapı, boyu, geçiş sayısı, kabuk çapı,
baffle aralığı) denenir. Uygun tasarımlar üç kademede sıralanır:

1. **TEMA zorunlu hükümleri:** desteksiz tüp boyu ve tüp demeti giriş/çıkış ρv² sınırı.
2. **Kritik ölçütler:** kabuk katsayısı %20 düşük çıksa bile yeterli pay (K2), su 5 °C sıcak gelse de yeterlilik (K3),
   gaz hiç yoğuşmasa bile kabuk basınç düşüşünün sınırda kalması (K9), akış kaynaklı titreşimin olmaması (K10).
3. **İyi uygulama ölçütleri:** güvenlik payı, basınç düşüşü payları, su hızı, tüp desteği, F faktörü, titreşimde %20 pay.

Eşitlik durumunda yıllık maliyeti en düşük tasarım seçilir.
""")
    gv.yaz(st.markdown, "m92", """
#### Hesaplar nasıl doğrulanır?
| Katman | Neyle karşılaştırılır |
|---|---|
| Denklemlerin uygulanışı | **Hesap adımları** sekmesinde her denklem sayılarıyla; el hesabıyla birebir aynı çıkmalı |
| Tüp sayısı | Phadke (1984) tam sayım yöntemi, Perry ve VDI bağıntıları |
| Bell–Delaware | Düzeltme çarpanları Taborek grafiklerinin bağımsız sayısallaştırmasıyla (ht kütüphanesi) |
| Alan | Bölgesel LMTD ↔ bölgesel ε-NTU (aynı çıkmalı), bölge sayısı iki katına çıkarılınca değişmemeli |
| Belirsizlik | Kern ↔ Bell–Delaware, iki faz sabiti C = 0–20, kabuk h ±%20 |
| Termodinamik | DWSIM (aynı NRTL parametreleriyle): ısı yükü, çiy noktası, çıkış gaz ve sıvı debileri |
""")
    c1_, c2_ = st.columns(2)
    with c1_:
        gv.yaz(st.markdown, "m84", """
#### Varsayımlar
- Kararlı rejim; çevreye ısı kaybı yoktur.
- Reaktörde yalnızca IPA → aseton + H₂ tepkimesi gerçekleşir; yan ürünler ihmal edilir.
- Eşanjöre giren akım tamamen buhardır; her noktada faz dengesi vardır ve hidrojen yoğuşmaz.
- Kabuk tarafı basıncı girişten çıkışa doğrusal olarak düşer (yoğuşma eğrisi bu basınçlarla hesaplanır).
- Kirlenme dirençleri ve boşluklar TEMA'dan alınır; malzeme karbon çeliğidir.
""")
    with c2_:
        gv.yaz(st.markdown, "m85", """
#### Sınırlamalar
- Kabuk tarafı katsayısı ve iki fazlı basınç düşüşü korelasyonlarının literatürdeki belirsizliği yüksektir
  (Bell–Delaware ΔP için ±%30–50); kritik ölçütler ve duyarlılık analizi bu nedenle uygulanır.
- DWSIM kabuk tarafını basitleştirilmiş Tinker yöntemiyle ve yoğuşmayı ayrıca modellemeden hesaplar; U ve ΔP
  karşılaştırması bilgi amaçlıdır, asıl doğrulama termodinamik katmandır.
- Titreşimde Strouhal sayısı yaklaşık alınmıştır; ayna kalınlığı ve gerilme analizi imalatçı tarafından kesinleştirilmelidir.
- Girişteki asetonun bir kısmı hidrojenle birlikte gaz olarak çıkar; geri kazanımı bu çalışmanın kapsamı dışındadır.
- Sonuçlar ön tasarım niteliğindedir.
""")


if "res" not in ss:
    gv.yaz(st.info, "m03", "**Nereden başlamalı?** Soldaki **1 · Proses akışı** ve **2 · Soğutma suyu** bölümlerine kendi verinizi girin "
            "(bileşen debileri, sıcaklıklar, basınçlar, su giriş sıcaklığı), sonra **🧭 Benim için tasarla**'ya basın. "
            "Uygulama geri kalan her şeyi standart aralıklardan seçer ve gerekçesini **🧭 Otomatik tasarım** sekmesinde açıklar.")
    tanitim_goster()
    st.stop()

R = ss["res"]
r, curve, m, c = R["r"], R["curve"], R["mech"], R["cost"]
z, ts, sg, lf, w, geo = r["zones"], r["ts"], r["sg"], r["lf"], r["water"], r["geo"]
lim = R["limits"]
fx = c["fx"]


cost_for = cost_fn_from(R)


def opt_fingerprint(R_):
    """Maliyet grafiklerini etkileyen girdiler; değişirse su çıkış taraması yeniden yapılır."""
    c_, l_ = R_["cost_in"], R_["limits"]
    return (round(float(R_["curve"]["Q_cum"].iloc[-1]), 1), round(float(R_["tube"]["T_in"]), 2), l_["v"], l_["dp_t"], l_["dp_s"],
            tuple(l_["od"]), R_["kw"].get("h_method"), R_["kw"].get("shells"), c_["correlation"], c_["water_price_m3"],
            c_["elec_price"], c_["i_rate"], c_["years"], c_["cepci_now"], c_["hours"])


def ensure_opt(R_):
    """Her hesaptan sonra maliyet grafiklerinin verisi güncel olsun: girdiler değiştiyse tarama otomatik yenilenir."""
    fp = opt_fingerprint(R_)
    O_ = ss.get("opt_res")
    if O_ is not None and O_.get("fp") == fp:
        return O_
    T_in_ = float(R_["tube"]["T_in"])
    T_list_ = list(range(int(math.ceil(T_in_ + 5)), 45 + 1))
    bar_ = st.progress(0.0, text="Girdiler değişti: grafikler için su çıkış sıcaklıkları yeniden taranıyor…")
    allc_, best_ = ec.optimize(R_["tube"], R_["curve"], T_list_, cost_fn_from(R_), R_["kw"], R_["limits"], clearance=R_["clearance"],
                               progress=lambda f, T: bar_.progress(f, text=f"Grafikler için su çıkışı {T} °C taranıyor…"))
    bar_.empty()
    top_ = allc_.head(5).copy()
    if not top_.empty:
        top_["asiri_tasarim_h_dusuk"] = [ec.robustness(R_["tube"], R_["curve"], rw, R_["kw"], 0.8) for _, rw in top_.iterrows()]
    O_ = dict(allc=allc_, best=best_, top=top_, T_list=T_list_, h_rob=0.8, fp=fp)
    ss["opt_res"] = O_
    return O_


def make_figs(R_, O_):
    """Rapor grafikleri: seçilen (girdilerdeki) tasarıma göre, her hesaptan sonra otomatik yeniden çizilir."""
    key = (id(R_), id(O_), round(R_["cost"]["TAC"], 2))
    cache = ss.get("fig_cache")
    if cache and cache[0] == key:
        return cache[1]
    r_, cst = R_["r"], R_["cost"]
    g_ = r_["geo"]
    T_out_ = float(r_["water"]["T_out"])
    tube_o = dict(R_["tube"], T_out=T_out_, m=None)
    sweep = ec.tac_vs_tubes_sweep(tube_o, R_["curve"], float(g_["Do"]), float(g_["t"]), int(g_["passes"]), (3048, 4877, 6096),
                                  cost_fn_from(R_), R_["kw"], R_["limits"], clearance=R_["clearance"], N_ref=g_["N"] * g_["L"] / 3048)
    figs = {"1_TAC_tup_sayisi": gr.tac_vs_tubes(sweep, "Seçilen tasarım", float(g_["Do"]), int(g_["passes"]), T_out_,
                                                selected=(g_["N"], cst["TAC"])),
            "2_sicaklik_isi_yuku_profili": gr.temperature_profile(r_["zones"], R_["T_dew"])}
    best_ = O_["best"]
    if not best_.empty:
        figs["3_yillik_toplam_maliyet"] = gr.tac_vs_water_temp(best_, O_["T_list"], selected=(T_out_, cst["TAC"]))
        figs["4_yatirim_isletme_dengesi"] = gr.capex_opex_balance(best_, O_["T_list"])
    pngs = {}
    for k_, f_ in figs.items():
        b_ = io.BytesIO()
        f_.savefig(b_, format="png", dpi=200)
        pngs[k_] = b_.getvalue()
    ss["fig_cache"] = (key, (figs, pngs))
    return figs, pngs


if ss.pop("opt_fp_pending", False) and "opt_res" in ss:
    ss["opt_res"]["fp"] = opt_fingerprint(R)           # otomatik tasarımın taraması bu hesapla aynı girdilere ait
O_cur = ensure_opt(R)
FIGS, PNGS = make_figs(R, O_cur)


def app_values(R_):
    """DWSIM ile karşılaştırılan uygulama değerleri (DWSIM sekmesi, Excel ve dwsim_girdi.json aynı tanımı kullanır)."""
    r_ = R_["r"]
    return dict(Q=r_["Q_kW"], T_dew=R_["T_dew"], m_V_out=float(R_["curve"].m_V.iloc[-1]), product_V_out=R_.get("prod_V"),
                dp_tube=r_["dp_tube_kPa"], U=tf.zonal_U(r_)["U_dirty"], dp_shell=r_["dp_shell_kPa"])


def money(usd):
    return f"${fmt(usd, 0)}" + (f"  (≈ {fmt(usd * fx, 0)} TL)" if fx else "")


T = gv.sekmeleri_kur(GA)

# ---------------------------------------------------------------- TANITIM
for _ in gv.sekme(T, "tanitim"):
    tanitim_goster()

# ---------------------------------------------------------------- OTOMATİK TASARIM
for _ in gv.sekme(T, "otomatik"):
    A = ss.get("auto_report")
    if not A:
        gv.yaz(st.info, "m04", "Otomatik tasarım için soldaki **Benim için tasarla** butonuna basın.")
    else:
        row, g_ = A["row"], A["geo"]
        if ss.pop("auto_done", False):
            gv.yaz(st.success, "m05", "Otomatik tasarım tamamlandı ve girdilere yazıldı. Tüm sekmeler bu tasarıma göre güncellendi.")
        mc = st.columns(4)
        mc[0].metric("Tüp (çap × boy)", f"Ø{g_['Do']} × {g_['L'] / 1000:.2f} m")
        mc[1].metric("Tüp sayısı · geçiş", f"{g_['N']} · {g_['passes']}")
        mc[2].metric("Kabuk çapı · alan", f"{g_['Ds']} mm · {r['A']:.0f} m²")
        mc[3].metric("Yıllık toplam maliyet", f"${fmt(c['TAC'], 0)}", help=f"Aşırı tasarım %{r['over_design']:.1f}")
        n_crit, n_all = len(ec.CRITICAL_CRITERIA), len(ec.TECH_CRITERIA)
        gv.yaz(st.caption, "m93", "Su girişi {su_girisleri} °C × su çıkışı × tüm standart geometriler denendi; kısıtları sağlayan "
               "{aday} tasarım değerlendirildi. Sıralama: önce TEMA zorunlu hükümleri, sonra kritik ölçütler "
               "(**{kritik}**), sonra toplam puan (**{puan}**); eşitlikte en düşük yıllık maliyet.",
               degerler=dict(su_girisleri=", ".join(str(t) for t in A["T_ins"]), aday=A["n_feasible"],
                             kritik=f"{int(row.kritik)}/{n_crit}", puan=f"{int(row.puan)}/{n_all}"))
        gv.yaz(st.subheader, "m06", "Teknik ölçütler (seçilen tasarım)")
        crit_vals = {
            "K1": f"%{row.asiri_tasarim:.1f}", "K2": f"%{row.asiri_tasarim_h_dusuk:.1f}", "K3": f"%{row.asiri_tasarim_su_sicak:.1f}",
            "K4": f"{row.dP_kabuk:.1f} / {A['limits']['dp_s']:.0f} kPa", "K5": f"{row.dP_tup:.1f} / {A['limits']['dp_t']:.0f} kPa",
            "K6": f"{row.v_tup:.2f} m/s", "K7": f"{row.aciklik_mm:.0f} / {ec.TEMA_MAX_SPAN.get(float(row.Do), 0):.0f} mm",
            "K8": f"{row.F:.3f}", "K9": f"{row.dP_kabuk_en_kotu:.1f} / {A['limits']['dp_s']:.0f} kPa",
            "K10": f"V/Vc = {row.V_Vc:.2f} ({row.titresim_kritik_aciklik})", "K11": f"V/Vc = {row.V_Vc:.2f}"}
        st.table(pd.DataFrame([(k, "Kritik" if k in ec.CRITICAL_CRITERIA else "İyi uygulama", lab, crit_vals[k],
                                "✅" if bool(row[k]) else "❌") for k, lab in ec.TECH_CRITERIA],
                              columns=["", "Tür", "Ölçüt", "Değer", "Sonuç"]).set_index(""))
        gv.yaz(st.caption, "m94", "TEMA zorunlu hükümleri: desteksiz açıklık {aciklik} mm (RCB-4.52) · demet giriş/çıkış ρv² "
               "{rhov2} kg/m·s² (RCB-4.62) → {sonuc}",
               degerler=dict(aciklik=f"{row.aciklik_mm:.0f} ≤ {ec.TEMA_MAX_SPAN.get(float(row.Do), 0):.0f}",
                             rhov2=f"{row.rhov2_demet_giris:.0f} / {row.rhov2_demet_cikis:.0f} ≤ {ec.RHOV2_BUNDLE_MAX:.0f}",
                             sonuc="✅ sağlanıyor" if row.tema_zorunlu else "❌ sağlanmıyor"))
        st.caption("Su giriş sıcaklığına göre: " + " · ".join(
            f"{b.T_su_giris:.0f} °C → {int(b.aday)} aday, en iyi kritik {int(b.en_yuksek_kritik)}/{n_crit}, "
            f"en yüksek puan {int(b.en_yuksek_puan)}/{n_all}" for b in A["by_inlet"].itertuples()))
        lm_ = A["limits"]
        gv.yaz(st.subheader, "m07", "Ne seçildi, neden?")
        dec = pd.DataFrame([
            ("Akışkan yerleşimi", "Soğutma suyu → tüp, yoğuşan buhar → kabuk",
             "Su kirletir ve tüp içi mekanik temizlenir; yoğuşma kabukta daha verimli, kabuk ΔP düşük.", "Towler & Sinnott, akışkan yerleşimi kuralları"),
            ("Tüp dış çapı", f"{g_['Do']} mm", "Denenen: 19,05 mm (3/4\") ve 25,4 mm (1\") → teknik ölçütlere göre seçildi.", "TEMA RCB-2.21 standart çap/et tablosu"),
            ("Et kalınlığı", f"{g_['t']} mm", "3/4\" için 14 BWG (2,11 mm), 1\" için 16 BWG (1,65 mm); iç basınç için yeterli (Mekanik sekmesi).", "Karbon çeliği tüp, yaygın uygulama"),
            ("Tüp adımı / dizilim", f"{g_['Pt']} mm ({g_['Pt'] / g_['Do']:.2f}·Do), üçgen", "En küçük adım 1,25·Do; kabuk tarafı mekanik temizlenmeyeceği için üçgen dizilim kullanılabilir.", "TEMA B-2.5 (adım ≥ 1,25·Do), RCB-2.42"),
            ("Tüp boyu", f"{g_['L']} mm", "Denenen standart boylar: 3048 / 4877 / 6096 mm (10 / 16 / 20 ft).", "TEMA RCB-2.1 yaygın tüp boyları"),
            ("Geçiş sayısı", f"{g_['passes']}", f"Denenen: 2 / 4 / 6 → tüp hızı {lm_['v'][0]}–{lm_['v'][1]} m/s aralığında kalanlar.", "Towler & Sinnott"),
            ("Tüp sayısı / kabuk çapı", f"{g_['N']} / {g_['Ds']} mm", "Ds 500–1400 mm (50 mm adım); N, tüplerin aynaya geometrik yerleşimiyle sayıldı (demet–kabuk boşluğu sabit aynada 15 mm, geçiş şeridi dahil).", "Geometrik ayna sayımı; TEMA RCB-4.3 boşlukları"),
            ("Baffle aralığı", f"{g_['B']} mm ({g_['B'] / g_['Ds']:.2f}·Ds)", "Denenen: 0,3–1,0·Ds. Alt sınır max(Ds/5, 51 mm); pencere tüplerinde destek açıklığı (2·B) ≤ TEMA sınırı; titreşimde V/Vc > 0,8 çıkarsa aralık daraltılarak yeniden denenir.", "TEMA RCB-4.51, RCB-4.52; TEMA Bölüm 6 (V-10, V-13.2)"),
            ("Su giriş sıcaklığı", f"{row.T_su_giris:.0f} °C", f"Denenen: {', '.join(str(t) for t in A['T_ins'])} °C (20 °C'den 5'er; giriş + 5 °C < proses çıkışı).", "Teknik ölçütler (K1–K11)"),
            ("Su çıkış sıcaklığı", f"{row.T_su_cikis:.0f} °C", f"Denenen: {A['T_list'][0]}–{A['T_list'][-1]} °C (giriş + 5 °C ile 45 °C arası; üst sınır kireçlenme).", "Teknik ölçütler (K1–K11)"),
            ("TEMA tipi", A["head"], A["tema"]["tema_reason"], "TEMA T-4 (metal sıcaklıkları), A.15 (sabit ayna, J = 1), E-4.12 (temizlik)"),
            ("Kirlenme dirençleri", "Kabuk 0,0002 · su 0,0002 m²K/W", "TEMA: arıtılmış kule suyu (>0,9 m/s, <52 °C) ve solvent buharları 0,001 h·ft²·°F/BTU ≈ 0,00018 m²K/W; 0,0002 ile küçük pay.", "TEMA T-2.4, RGP-T-2.4"),
            ("Tüp hızı", f"{lm_['v'][0]}–{lm_['v'][1]} m/s", "1 m/s altında kirlenme artar; üstte aşınma ve ΔP.", "Towler & Sinnott"),
            ("Tüp ΔP sınırı", f"{lm_['dp_t']:.0f} kPa", "Viskozitesi < 1 mPa·s sıvılar için tipik izin verilen değer.", "Towler & Sinnott"),
            ("Kabuk ΔP sınırı", f"{lm_['dp_s']:.0f} kPa", "Proses verisinden: giriş – çıkış basıncı.", "Sizin akış şemanız"),
            ("Aşırı tasarım", f"%{lm_['od'][0]:.0f}–{lm_['od'][1]:.0f} (tercih %20–35)", "H₂'li yoğuşmada kabuk h belirsiz → kabul için ≥ %15, tercih %20–35 (K1); h %20 düşükken ≥ %10 (K2).", "Tasarım uygulaması"),
            ("Kabuk h yöntemi", "Bell–Delaware (gaz) + Nusselt (film) + Colburn–Hougen (H₂)", "Gaz tarafında kaçak ve dolaşan akışlar Bell–Delaware ile; yoğuşmayan H₂'nin difüzyon direnci Colburn–Hougen ile hesaba katıldı.", "TEMA T-1.1 (Delaware Bülten 5); Colburn & Hougen (1934)"),
            ("Kabuk ΔP yöntemi", "Bölgesel iki faz (C = 8) + en kötü durum kontrolü", "Her bölge yerel gaz/sıvı debisiyle; ayrıca gaz hiç yoğuşmasa bile sınır aşılmamalı (K9).", "Ishihara vd. (1980); TEMA RGP-T-2.4 (belirsizlikler açıkça gösterilmeli)"),
            ("Seçim kuralı", f"Kritik {int(row.kritik)}/{n_crit} · toplam {int(row.puan)}/{n_all}", "Önce TEMA zorunlu hükümleri (RCB-4.52, RCB-4.62), sonra kritik ölçütler (K2, K3, K9, K10 titreşim), sonra toplam puan; eşitlikte en düşük yıllık maliyet.", "TEMA 'shall' hükümleri zorunlu; RGP-T-2.4"),
        ], columns=["Karar", "Seçilen", "Kural / denenen aralık", "Dayanak"])
        st.table(dec.set_index("Karar"))
        ss["decisions_df"] = dec
        gv.yaz(st.subheader, "m08", "En uygun 5 alternatif (teknik puana göre)")
        t5 = A["top"][["kritik", "puan", "T_su_giris", "T_su_cikis", "Do", "L", "gecis", "Ds", "B", "N", "A", "asiri_tasarim",
                       "asiri_tasarim_h_dusuk", "asiri_tasarim_su_sicak", "v_tup", "dP_tup", "dP_kabuk", "dP_kabuk_en_kotu",
                       "V_Vc", "yillik_toplam_USD"]].copy()
        t5.insert(0, "Sıra", range(1, len(t5) + 1))
        t5["Seçilen"] = ["◀" if i == 0 else "" for i in range(len(t5))]
        st.dataframe(t5.round(2).assign(Do=t5.Do.astype(str)).rename(columns={
            "asiri_tasarim_h_dusuk": "aşırı tas. (h×0,8)", "asiri_tasarim_su_sicak": "aşırı tas. (su +5 °C)"}),
            hide_index=True, width="stretch")
        gv.yaz(st.markdown, "m09", """
**Sonraki adımlar**
1. **📋 Özet** ve **🔩 Mekanik** sekmelerinde tüm kriterlerin ✅ olduğunu kontrol edin.
2. **🔁 Bölgeler & doğrulama** sekmesinde hidrolik ve ısıl doğrulama tablolarını inceleyin.
3. **🧪 DWSIM karşılaştırma** sekmesinde bu tasarımı ve diğer 4 alternatifi DWSIM ile karşılaştırın.
4. **📄 TEMA föyü**'nden föyü (PDF), tüm hesapları (Excel) ve grafikleri indirin.
""")

# ---------------------------------------------------------------- ÖZET
for _ in gv.sekme(T, "ozet"):
    cols = st.columns(4)
    cols[0].metric("Isı yükü Q", f"{fmt(r['Q_kW'], 0)} kW")
    cols[1].metric("Mevcut alan", f"{r['A']:.1f} m²")
    cols[2].metric("Gerekli alan", f"{r['A_req']:.1f} m²")
    cols[3].metric("Aşırı tasarım", f"%{r['over_design']:.1f}")
    cols = st.columns(4)
    Uz = tf.zonal_U(r)
    cols[0].metric("U kirli (bölgesel)", f"{Uz['U_dirty']:.0f} W/m²K",
                   help=f"Bölgesel hesaptan: UA/A_gerekli. Bölgesel MTD (F dahil) {Uz['mtd']:.1f} °C. "
                        f"Tek LMTD·F ile eşdeğer değer {r['Ud']:.0f} W/m²K (yoğuşma eğrisi doğrusal olmadığı için farklı).")
    cols[1].metric("Tüp ΔP / Kabuk ΔP", f"{r['dp_tube_kPa']:.1f} / {r['dp_shell_kPa']:.1f} kPa")
    cols[2].metric("Kurulu yatırım", f"${fmt(c['C_installed'] / 1000, 0)} bin")
    cols[3].metric("Yıllık toplam maliyet", f"${fmt(c['TAC'] / 1000, 0)} bin/yıl")

    gv.yaz(st.subheader, "m10", "Kriter kontrolü")
    st.dataframe(ec.checks(r, lim["v"], lim["dp_t"], lim["dp_s"], lim["od"], m, R["clearance"]), hide_index=True,
                 width="stretch")
    if r["over_design"] > lim["od"][1]:
        gv.yaz(st.warning, "m95", "Eşanjör **fazla büyük** (aşırı tasarım %{od}). **Tasarım taraması** sekmesinden daha küçük "
               "bir tasarım seçin.", degerler=dict(od=f"{r['over_design']:.0f}"))
    if r["over_design"] < lim["od"][0]:
        gv.yaz(st.warning, "m96", "Alan **yetersiz**. Tüp sayısı/boyu artırılmalı veya taramadan geometri seçilmeli.")
    if liquid_model == "ideal" and curve_source.startswith("Bileşim"):
        gv.yaz(st.warning, "m97", "Sıvı modeli **ideal**: su–aseton ve IPA–su ideal değildir; çiy noktası ~10 °C yüksek çıkar. "
               "NRTL önerilir.")
    if R.get("prod_V") is not None:
        _feed = dict(zip(R["comps"], R["flows"])).get(product, 0.0)
        gv.yaz(st.info, "m98", "Çıkış gazıyla **{kayip} kg/h {urun}** kaçıyor (girişin %{yuzde}'i). Soğutma suyuyla daha düşük "
               "çıkış sıcaklığına inilemediği için bu gaz bir absorber veya soğutulmuş su ile ayrıca işlenmelidir.",
               degerler=dict(kayip=f"{R['prod_V']:.0f}", urun=product, yuzde=f"{R['prod_V'] / _feed * 100:.1f}" if _feed else "–"))
    gv.yaz(st.caption, "m99", "Kabuk h yöntemi: {yontem} · Bölge sayısı: {bolge} · Seçilen TEMA tipi: {tema}",
           degerler=dict(yontem=h_method, bolge=len(z), tema=m["tema"]))

# ---------------------------------------------------------------- TEMA FÖYÜ
for _ in gv.sekme(T, "tema"):
    D = tf.DEFAULT_CHOICES
    with st.expander("Föy bilgileri ve malzeme/standart seçimleri (hesaplanmayan alanlar)"):
        gv.yaz(st.caption, "m11", "Isıl, hidrolik ve mekanik değerler hesaplardan otomatik gelir. Burada yalnız hesapla bulunmayan bilgileri "
                   "girin; varsayılanlar karbon çeliği kimya prosesi eşanjörü için tipiktir.")
        c1, c2, c3 = st.columns(3)
        chs = dict(
            customer=c1.text_input("Müşteri", key=K("tf_customer", proj["project"])),
            address=c2.text_input("Adres", key=K("tf_address", "")),
            location=c3.text_input("Tesis konumu", key=K("tf_location", proj["plant"])),
            job_no=c1.text_input("İş no", key=K("tf_job_no", "")),
            ref_no=c2.text_input("Referans / PO no", key=K("tf_ref_no", "")),
            proposal_no=c3.text_input("Teklif no", key=K("tf_proposal_no", "")),
            shell_fluid=c1.text_input("Kabuk akışkanı adı", key=K("tf_shell_fluid", D["shell_fluid"])),
            tube_fluid=c2.text_input("Tüp akışkanı adı", key=K("tf_tube_fluid", D["tube_fluid"])),
            mdmt=c3.text_input("MDMT (°C), bilinmiyorsa boş", key=K("tf_mdmt", "")),
        )
        gv.yaz(st.markdown, "m12", "**Kod ve standart**")
        c1, c2, c3, c4 = st.columns(4)
        chs.update(
            code=c1.text_input("Tasarım kodu", key=K("tf_code", D["code"])),
            tema_class=c2.selectbox("TEMA sınıfı", ["B", "C", "R"], key=K("tf_tema_class", "B"),
                                    help="R: petrol rafinerisi, C: genel ticari, B: kimya prosesi"),
            radiography=c3.selectbox("Radyografi", ["Nokta (E = 0,85)", "Tam (E = 1,0)", "Yok"], key=K("tf_radiography", D["radiography"])),
            special_service=c4.text_input("Özel hizmet", key=K("tf_special", D["special_service"])),
            flange_class=c1.selectbox("Flanş sınıfı", ["150", "300", "600"], key=K("tf_flange", "150")),
            facing=c2.selectbox("Flanş yüzeyi", ["RF", "RTJ", "FF"], key=K("tf_facing", "RF")),
            CA_tube=c3.number_input("Tüp tarafı korozyon payı (mm)", key=K("tf_CA_tube", 3.0)),
            pwht=c4.checkbox("PWHT gerekli", key=K("tf_pwht", False)),
            code_stamp=c1.checkbox("Kod damgası", key=K("tf_stamp", True)),
            mawp=c2.checkbox("MAWP hesaplansın", key=K("tf_mawp", True)),
            inspection=c3.checkbox("Müşteri muayenesi", key=K("tf_insp", True)),
        )
        gv.yaz(st.markdown, "m13", "**Konstrüksiyon**")
        c1, c2, c3, c4 = st.columns(4)
        chs.update(
            baffle_cut=c1.number_input("Baffle kesimi (%)", min_value=15, max_value=45, key=K("tf_cut", 25)),
            cut_dir=c2.selectbox("Kesim yönü", ["Dikey", "Yatay"], key=K("tf_cut_dir", "Dikey"),
                                 help="Yoğuşmada dikey kesim: kondensat baffle altından akar."),
            tube_joint=c3.selectbox("Tüp–ayna bağlantısı", ["Genişletme + oluk", "Genişletme", "Kaynak + genişletme"], key=K("tf_joint", D["tube_joint"])),
            n_supports=c4.number_input("Mesnet sayısı", min_value=2, max_value=4, key=K("tf_nsup", 2)),
            vent_shell=c1.text_input("Kabuk havalık", key=K("tf_vs", D["vent_shell"])),
            drain_shell=c2.text_input("Kabuk drenaj", key=K("tf_ds", D["drain_shell"])),
            vent_tube=c3.text_input("Tüp havalık", key=K("tf_vt", D["vent_tube"])),
            drain_tube=c4.text_input("Tüp drenaj", key=K("tf_dt", D["drain_tube"])),
        )
        gv.yaz(st.markdown, "m14", "**Malzemeler**")
        c1, c2, c3, c4 = st.columns(4)
        chs.update(
            mat_tube=c1.text_input("Tüpler", key=K("tf_m_tube", tf.MATERIALS.get(material, tf.MATERIALS["Karbon çeliği"])[0])),
            mat_shell=c2.text_input("Kabuk", key=K("tf_m_shell", tf.MATERIALS.get(material, tf.MATERIALS["Karbon çeliği"])[1])),
            mat_tubesheet=c3.text_input("Ayna", key=K("tf_m_ts", tf.MATERIALS.get(material, tf.MATERIALS["Karbon çeliği"])[1])),
            mat_channel=c4.text_input("Kanal / bonnet", key=K("tf_m_ch", tf.MATERIALS.get(material, tf.MATERIALS["Karbon çeliği"])[1])),
            mat_baffle=c1.text_input("Baffle'lar", key=K("tf_m_baf", D["mat_baffle"])),
            mat_tierod=c2.text_input("Bağlantı çubukları / mesnet", key=K("tf_m_tie", D["mat_tierod"])),
            mat_flange=c3.text_input("Flanşlar", key=K("tf_m_fl", D["mat_flange"])),
            mat_neck=c4.text_input("Nozul boyunları", key=K("tf_m_neck", D["mat_neck"])),
            gasket=c1.text_input("Contalar", key=K("tf_gasket", D["gasket"])),
            bolting=c2.text_input("Cıvatalar", key=K("tf_bolt", D["bolting"])),
            insulation=c3.text_input("İzolasyon", key=K("tf_ins", D["insulation"])),
            paint=c4.text_input("Boya", key=K("tf_paint", D["paint"])),
        )
    inp = dict(proj=proj, head_type=head_type, P_in_bar=P_in_bar, P_w_bar=w_P, dp_shell_max=lim["dp_s"], dp_tube_max=lim["dp_t"],
               P_d_shell=P_d_shell, P_d_tube=P_d_tube, CA_shell=CA_shell,
               thermo_note=("DWSIM eğrisi (tablo)" if R["use_table"] else
                            "gaz Peng–Robinson, sıvı " + ("NRTL (NIST/OSTI seti)" if liquid_model == "nrtl" else "ideal (Raoult)")))
    fd = tf.collect(R, inp, chs)
    item = proj["item"] or "esanjor"

    # Excel'e eklenecek ayrıntı tabloları
    dw_vals = {k: ss.get(f"dw_{s_}", 0.0) for k, s_ in
               (("Q", "Q"), ("T_dew", "Tdew"), ("m_V_out", "mV"), ("product_V_out", "prod"), ("U", "U"), ("dp_tube", "dpt"), ("dp_shell", "dps"))}
    app_vals = app_values(R)
    sheets = {
        "Kriter kontrolü": ec.checks(r, lim["v"], lim["dp_t"], lim["dp_s"], lim["od"], m, R["clearance"]),
        "Ne seçildi, neden": ss.get("decisions_df") if ss.get("auto_report") else None,
        "Yoğuşma eğrisi": ec.curve_to_table(curve),
        "Bölgeler": z,
        "Duyarlılık": R["sens"],
        "Nozullar": m["nozzles"],
        "Optimizasyon ilk 5": ss["opt_res"]["top"] if "opt_res" in ss else None,
        "Su sıcaklığına göre en iyi": ss["opt_res"]["best"] if "opt_res" in ss else None,
        "Doğrulama - tüp sayısı": R["hv"]["tup_sayisi"],
        "Doğrulama - tüp ΔP": R["hv"]["tup_dp"],
        "Doğrulama - Bell-Delaware": R["hv"]["bd_carpanlar"],
        "Doğrulama - kabuk ΔP": R["hv"]["kabuk_dp"],
        "DWSIM karşılaştırma": ec.dwsim_compare(app_vals, dw_vals, ec.dwsim_bounds(r, R["hv"])) if any(dw_vals.values()) else None,
    }
    gv.yaz(st.caption, "m15", "Föy, uygulamanın hesapladığı değerlerle otomatik doldurulur; girdileri değiştirip yeniden hesapladığınızda "
               "föy de güncellenir. Proje başlığı sol paneldeki **0 · Proje bilgileri** bölümünden gelir.")
    pdf_tr_, xls_ = tf.pdf_tr(fd), tf.to_excel(fd, R, sheets)
    c1, c3 = st.columns(2)
    c1.download_button("⬇ Föy — PDF (Türkçe)", pdf_tr_, f"{item}_TEMA_foyu_TR.pdf", mime="application/pdf",
                       type="primary", width="stretch")
    c3.download_button("⬇ Tüm hesaplar — Excel", xls_, f"{item}_hesaplar.xlsx",
                       mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet", width="stretch",
                       help="1. sayfa TEMA föyü, 2. sayfa tüm hesaplanan veriler (kaynağıyla), sonra ayrıntı tabloları.")
    import zipfile
    zb = io.BytesIO()
    with zipfile.ZipFile(zb, "w") as zf:
        zf.writestr(f"{item}_TEMA_foyu_TR.pdf", pdf_tr_)
        zf.writestr(f"{item}_hesaplar.xlsx", xls_)
        for k, b_ in PNGS.items():
            zf.writestr(f"grafikler/{k}.png", b_)
    st.download_button("📦 Rapor paketi — föy + Excel + grafikler (ZIP)", zb.getvalue(), f"{item}_rapor_paketi.zip",
                       mime="application/zip", type="primary", width="stretch",
                       help="Her hesaptan sonra otomatik yeniden hazırlanır.")
    st.markdown(tf.to_html(fd), unsafe_allow_html=True)


# ---------------------------------------------------------------- YOĞUŞMA EĞRİSİ
for _ in gv.sekme(T, "egri"):
    Q = r["Q_kW"]
    Q_dew = float(np.interp(R["T_dew"], curve["T"][::-1], curve["Q_cum"][::-1]))
    c1, c2, c3 = st.columns(3)
    c1.metric("Çiy noktası", f"{R['T_dew']:.1f} °C")
    c2.metric("Kızgın buhar soğuması", f"{Q_dew:.0f} kW (%{Q_dew / Q * 100:.0f})")
    c3.metric("Çıkışta yoğuşan", f"%{curve.m_L.iloc[-1] / (curve.m_V.iloc[0] + curve.m_L.iloc[0]) * 100:.0f}")
    fig, ax = plt.subplots(figsize=(8, 4.2))
    ax.plot(curve.Q_cum, curve["T"], color=BLUE, lw=2, marker="o", ms=4, label="Sıcak akış (yoğuşma eğrisi)")
    ax.plot([0, Q], [curve["T"].iloc[0], curve["T"].iloc[-1]], color=ORANGE, lw=2, ls="--", label="Tek LMTD kabulü (doğrusal)")
    ax.plot([0, Q], [w["T_out"], w["T_in"]], color=GREEN, lw=2, label="Soğutma suyu (ters akış)")
    ax.axhline(R["T_dew"], color=GRAY, lw=1, ls=":")
    ax.annotate(f"çiy noktası {R['T_dew']:.1f} °C", (Q * 0.55, R["T_dew"] + 1.5), color="#555", fontsize=9)
    ax.set(xlabel="Çekilen ısı, Q [kW]", ylabel="Sıcaklık [°C]")
    ax.grid(color="#e6e6e3", lw=0.8)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False)
    st.pyplot(fig)
    gv.yaz(st.markdown, "m16", "Eğrinin doğrudan (turuncu) sapması, **tek LMTD'nin neden yetersiz olduğunu** gösterir: "
                "ısının büyük kısmı çiy noktasının hemen altında çekilir.")
    table = ec.curve_to_table(curve)
    st.dataframe(table.round(4), hide_index=True, width="stretch")
    st.download_button("⬇ Eğriyi CSV indir (tablo girişi için şablon)", table.to_csv(index=False), "yogusma_egrisi.csv")

# ---------------------------------------------------------------- HESAP ADIMLARI
for _ in gv.sekme(T, "adimlar"):
    gv.yaz(st.markdown, "m17", "Tüm hesaplar formülleriyle ve sayılarıyla. Bölge hesaplarında bir bölge seçip ayrıntısını görün.")

    gv.yaz(st.subheader, "m18", "1. Isı yükü ve su debisi")
    st.latex(rf"Q = H(T_{{\text{{giriş}}}}) - H(T_{{\text{{çıkış}}}}) = {r['Q_kW']:.1f}\ \text{{kW}} \quad \text{{(yoğuşma eğrisinden)}}")
    st.latex(rf"\dot m_{{\text{{su}}}} = \frac{{Q}}{{c_p\,(t_2 - t_1)}} = \frac{{{r['Q_kW']:.1f}}}{{{w['cp']}\times({w['T_out']:.1f}-{w['T_in']:.1f})}}"
             rf"\times 3600 = {w['m']:.0f}\ \text{{kg/h}}")

    gv.yaz(st.subheader, "m19", "2. LMTD ve F düzeltme faktörü")
    Th1, Th2 = r["T_hot"]
    st.latex(rf"\Delta T_{{lm}} = \frac{{(T_1 - t_2) - (T_2 - t_1)}}{{\ln\frac{{T_1 - t_2}}{{T_2 - t_1}}}} = "
             rf"\frac{{({Th1:.1f}-{w['T_out']:.1f}) - ({Th2:.1f}-{w['T_in']:.1f})}}{{\ln\frac{{{Th1 - w['T_out']:.1f}}}{{{Th2 - w['T_in']:.1f}}}}} = {lf['LMTD']:.2f}\ \text{{K}}")
    st.latex(rf"R = \frac{{T_1 - T_2}}{{t_2 - t_1}} = {lf['R']:.3f}, \qquad P = \frac{{t_2 - t_1}}{{T_1 - t_1}} = {lf['P']:.4f}")
    st.latex(r"F = \frac{\sqrt{R^2+1}\,\ln\frac{1-P}{1-RP}}{(R-1)\ln\frac{2-P\left(R+1-\sqrt{R^2+1}\right)}{2-P\left(R+1+\sqrt{R^2+1}\right)}}"
             rf" = {lf['F']:.3f}")
    gv.yaz(st.caption, "m20", "F, tüm eşanjör için hesaplanır ve her bölgeye uygulanır (Towler & Sinnott; muhafazakâr).")

    gv.yaz(st.subheader, "m21", "3. Tüp tarafı (soğutma suyu)")
    st.latex(rf"D_i = D_o - 2t = {geo['Do']} - 2\times{geo['t']} = {ts['Di'] * 1000:.2f}\ \text{{mm}}")
    st.latex(rf"A_{{\text{{akış}}}} = \frac{{\pi}}{{4}} D_i^2 \frac{{N}}{{N_p}} = \frac{{\pi}}{{4}}({ts['Di']:.5f})^2 \frac{{{geo['N']}}}{{{geo['passes']}}} = {ts['A_flow']:.4f}\ \text{{m}}^2")
    st.latex(rf"v = \frac{{\dot m}}{{\rho A_{{\text{{akış}}}}}} = {ts['v']:.3f}\ \text{{m/s}}, \quad "
             rf"Re = \frac{{\rho v D_i}}{{\mu}} = {ts['Re']:.0f}, \quad Pr = \frac{{c_p \mu}}{{k}} = {ts['Pr']:.2f}")
    st.latex(rf"f_D = {ts['fd']:.5f}\ \text{{(Colebrook)}}, \qquad Nu = \frac{{(f_D/8)(Re-1000)Pr}}{{1+12.7\sqrt{{f_D/8}}\,(Pr^{{2/3}}-1)}} = {ts['Nu']:.1f}\ \text{{(Gnielinski)}}")
    st.latex(rf"h_i = \frac{{Nu\,k}}{{D_i}} = {ts['h']:.0f}\ \text{{W/m}}^2\text{{K}}")
    st.latex(rf"\Delta P_t = N_p\left(f_D \frac{{L}}{{D_i}} + K\right)\frac{{\rho v^2}}{{2}} = {geo['passes']}\left({ts['fd']:.4f}\frac{{{sg['L']:.3f}}}{{{ts['Di']:.4f}}} + {ts['K_pass']}\right)"
             rf"\frac{{{w['rho']:.0f}\times{ts['v']:.3f}^2}}{{2}} = {ts['dp_kPa']:.1f}\ \text{{kPa}}")

    gv.yaz(st.subheader, "m22", "4. Kabuk geometrisi")
    st.latex(rf"D_b = D_s - \text{{boşluk}} = {geo['Ds']:.0f} - {geo['Ds'] - sg['Db'] * 1000:.0f} = {sg['Db'] * 1000:.0f}\ \text{{mm}}")
    gv.yaz(st.caption, "m86", "Tüp sayısı, tüplerin bu demet çapı içine üçgen adımla tek tek yerleştirilmesiyle (geçiş şeritleri "
                              "boş bırakılarak) sayılır; tahmini bir bağıntı kullanılmaz.")
    st.latex(rf"D_e = {sg['De'] * 1000:.2f}\ \text{{mm}}, \quad A_s = \frac{{D_s B (P_t - D_o)}}{{P_t}} = {sg['As']:.4f}\ \text{{m}}^2, \quad "
             rf"N_r = \tfrac{{2}}{{3}}\frac{{D_b}}{{P_t}} = {sg['Nr']:.1f}, \quad n_{{\text{{baffle}}}} = {sg['n_baffle']}")

    gv.yaz(st.subheader, "m23", "5. Bölge hesabı (seçilen bölge)")
    labels = [f"{i + 1}: {a:.1f} → {b:.1f} °C" for i, (a, b) in enumerate(zip(z.T_a, z.T_b))]
    default = int(np.argmax(z.dQ_kW))
    i = labels.index(st.selectbox("Bölge", labels, index=default))
    q = z.iloc[i]
    st.latex(rf"\Delta Q = {q.dQ_kW:.1f}\ \text{{kW}}, \quad \dot m_V = {q.m_V:.0f}, \ \dot m_L = {q.m_L:.0f}\ \text{{kg/h}}, \quad "
             rf"t_{{\text{{su}}}}: {q.t_su_a:.2f} \to {q.t_su_b:.2f}\ ^\circ\text{{C}}")
    if r.get("shell_method") == "bell_delaware" and r.get("bd"):
        gv.yaz(st.markdown, "m24", "**Gaz (duyulur) — Bell–Delaware**")
        _c0 = curve.iloc[0]
        _bd0 = ec.bell_delaware(_c0.m_V / 3600, _c0.rho_V, _c0.mu_V, _c0.cp_V, _c0.k_V, r["bd"])
        st.latex(rf"h_G = h_{{\text{{ideal}}}}\,J_c\,J_l\,J_b\,J_s\,J_r = {q.h_G:.1f}\ \text{{W/m}}^2\text{{K}}"
                 rf"\quad (\text{{Kern ile: }} {q.h_G_kern:.1f})")
        gv.yaz(st.caption, "m138", "Düzeltme faktörleri (kabuk girişi): Jc = {Jc} (pencere), Jl = {Jl} (baffle boşluklarından "
               "kaçak), Jb = {Jb} (demet çevresinden dolaşan akış), Js = {Js} (uç bölmeler), Jr = {Jr}. Çapraz akış alanı "
               "Sm = {Sm} m².", degerler=dict(**{k_: f"{_bd0[k_]:.3f}" for k_ in ("Jc", "Jl", "Jb", "Js", "Jr")},
                                              Sm=f"{r['bd']['Sm']:.4f}"))
    else:
        gv.yaz(st.markdown, "m87", "**Gaz (duyulur) — Kern**")
        st.latex(rf"Re_V = \frac{{D_e (\dot m_V / A_s)}}{{\mu_V}} = {q.Re_v:.0f}, \quad Pr_V = {q.Pr_v:.3f}, \quad "
                 rf"h_G = 0.36\,Re^{{0.55}}Pr^{{1/3}}\frac{{k_V}}{{D_e}} = {q.h_G:.1f}\ \text{{W/m}}^2\text{{K}}")
    if np.isfinite(q.h_F):
        gv.yaz(st.markdown, "m25", "**Kondensat filmi — Nusselt–Kern**")
        st.latex(rf"\Gamma = \frac{{\dot m_L}}{{L N_t}} = {q.Gamma:.5f}\ \text{{kg/m·s}}, \quad "
                 rf"h_F = 0.95\,k_L\left[\frac{{\rho_L(\rho_L-\rho_V)g}}{{\mu_L\Gamma}}\right]^{{1/3}}N_r^{{-1/6}} = {q.h_F:.0f}\ \text{{W/m}}^2\text{{K}}")
        gv.yaz(st.markdown, "m26", "**SBG**")
        st.latex(rf"Z = \frac{{\dot m_V c_{{p,V}} \Delta T}}{{\Delta Q}} = {q.Z:.3f}, \quad "
                 rf"\frac{{1}}{{h_{{SBG}}}} = \frac{{1}}{{h_F}} + \frac{{Z}}{{h_G}} \Rightarrow h_{{SBG}} = {q.h_SBG:.0f}\ \text{{W/m}}^2\text{{K}}")
        if np.isfinite(q.T_i):
            gv.yaz(st.markdown, "m27", "**Colburn–Hougen**")
            st.latex(rf"Sc = \frac{{\mu_V}}{{\rho_V D_{{AB}}}} = {q.Sc:.3f}, \quad K = \frac{{h_G}}{{\rho_V c_{{p,V}}}}\left(\frac{{Pr}}{{Sc}}\right)^{{2/3}}c = {q.K:.3f}\ \text{{mol/m}}^2\text{{s}}, \quad \lambda = {q.lam / 1000:.0f}\ \text{{kJ/kg}}")
            st.latex(rf"h_G(T_G - T_i) + K\ln\frac{{1-y_i}}{{1-y_G}}\,M_c\lambda = U_w(T_i - t_{{\text{{su}}}}) \Rightarrow "
                     rf"T_i = {q.T_i:.2f}\ ^\circ\text{{C}}, \ y_G = {q.y_G:.3f}, \ y_i = {q.y_i:.3f}")
            st.latex(rf"h_{{\text{{etkin}}}} = {q.get('h_CH', np.nan):.0f}\ \text{{W/m}}^2\text{{K}}")
    else:
        gv.yaz(st.info, "m28", "Bu bölge çiy noktasının üstünde: sıvı yok, kuru duvar kabulüyle h = h_G.")
    gv.yaz(st.markdown, "m29", "**Toplam U, alan ve ε-NTU kontrolü**")
    dr = sg["Do"] / ts["Di"]
    st.latex(rf"\frac{{1}}{{U_d}} = \frac{{1}}{{h_{{\text{{kabuk}}}}}} + R_{{f,k}} + \frac{{D_o\ln(D_o/D_i)}}{{2k_w}} + \frac{{D_o}}{{D_i}}R_{{f,t}} + \frac{{D_o}}{{D_i}}\frac{{1}}{{h_i}}"
             rf"= \frac{{1}}{{{q.h_shell:.0f}}} + {r['Rf_shell']} + {r['R_wall']:.2e} + {dr:.3f}\times{w['Rf']} + \frac{{{dr:.3f}}}{{{ts['h']:.0f}}} \Rightarrow U_d = {q.Ud:.0f}")
    st.latex(rf"A_i = \frac{{\Delta Q}}{{U_d F \Delta T_{{lm,i}}}} = \frac{{{q.dQ_kW * 1000:.0f}}}{{{q.Ud:.0f}\times{lf['F']:.3f}\times{q.LMTD:.2f}}} = {q.A_m2:.2f}\ \text{{m}}^2")
    st.latex(rf"\varepsilon = \frac{{\Delta Q}}{{C_{{min}}(T_{{\text{{h,giriş}}}} - t_{{\text{{c,giriş}}}})}} = {q.eps:.4f}, \ C_r = {q.Cr:.3f} \Rightarrow NTU = {q.NTU:.4f}, \quad "
             rf"A_{{\varepsilon\text{{-}}NTU}} = \frac{{NTU\,C_{{min}}}}{{U_d}} = {q.A_eNTU_m2:.2f}\ \text{{m}}^2 \ (= A_i F = {q.A_m2 * lf['F']:.2f})")

    gv.yaz(st.subheader, "m30", "6. Toplam alan ve aşırı tasarım")
    st.latex(rf"A_{{\text{{gerekli}}}} = \sum A_i = {r['A_req']:.1f}\ \text{{m}}^2, \quad A = \pi D_o L N = \pi\times{sg['Do']}\times{sg['L']}\times{geo['N']} = {r['A']:.1f}\ \text{{m}}^2")
    st.latex(rf"\text{{Aşırı tasarım}} = \frac{{A}}{{A_{{\text{{gerekli}}}}}} - 1 = \%{r['over_design']:.1f}")

    gv.yaz(st.subheader, "m31", "7. Kabuk ΔP (iki faz)")
    _C = r.get("two_phase_C", 8.0)
    _pp = r.get("dp_shell_parts") or {}
    st.latex(rf"X = \sqrt{{\frac{{\Delta P_L}}{{\Delta P_G}}}} = {q.X:.3f}, \quad \phi_G^2 = 1 + C X + X^2 = 1 + {_C:g}\times{q.X:.3f} + {q.X:.3f}^2 = {q.phi2:.3f}")
    st.latex(rf"\Delta P_{{\text{{kabuk}}}} = \underbrace{{\sum_i \phi_{{G,i}}^2\,\Delta P_{{G,i}}\,\frac{{A_i}}{{\sum A}}}}_{{\text{{bölmeler ve pencereler}}}}"
             rf" + \Delta P_{{\text{{uç bölmeler}}}} + \Delta P_{{\text{{nozullar}}}}"
             rf" = {_pp.get('bolmeler_pencereler', float('nan')):.2f} + {_pp.get('giris_cikis_bolmeleri', float('nan')):.2f}"
             rf" + {_pp.get('nozullar', float('nan')):.2f} = {r['dp_shell_kPa']:.2f}\ \text{{kPa}}")
    gv.yaz(st.caption, "m32", "ΔP_G ve ΔP_L: her bölgede ilgili faz tek başına akıyormuş gibi Bell–Delaware yöntemiyle; kaçak ve "
               "dolaşan akış etkileri dahil. C = 8, tüp demeti üzerinden çapraz akış için kullanılan değerdir. Her bölgenin "
               "kaybı gerekli alandaki payı kadar kabuk boyuna uygulanır; giriş/çıkış bölmeleri ilk/son bölge koşullarıyla, "
               "nozul kayıpları 1,0 + 0,5 hız yüküyle eklenir. Yoğuşmadan kaynaklanan küçük basınç kazancı ihtiyatlı olmak "
               "için düşülmez.")

# ---------------------------------------------------------------- BÖLGELER & DOĞRULAMA
for _ in gv.sekme(T, "bolgeler"):
    gv.yaz(st.subheader, "m33", "Bölge tablosu")
    show = z[["T_a", "T_b", "dQ_kW", "t_su_a", "t_su_b", "h_G", "h_F", "Z", "T_i", "h_shell", "Ud", "LMTD",
              "A_m2", "A_eNTU_m2", "dP_2f_kPa", "model"]]
    st.dataframe(show.round(2), hide_index=True, width="stretch")

    gv.yaz(st.subheader, "m34", "ε-NTU doğrulaması")
    gv.yaz(st.markdown, "m110", "Aynı bölgeler ve aynı U değerleriyle iki farklı formülasyon (F = 1, ters akış):  \n"
           "- Bölgesel LMTD: **{A_lmtd} m²**  \n- Bölgesel ε-NTU: **{A_entu} m²**  \n"
           "Aynı çıkması, alan hesabının doğru kurulduğunu gösterir. Tasarım alanı F = {F} ile: **{A_req} m²**.",
           degerler=dict(A_lmtd=f"{r['A_LMTD_cf']:.2f}", A_entu=f"{r['A_eNTU_cf']:.2f}", F=f"{lf['F']:.3f}", A_req=f"{r['A_req']:.1f}"))

    gv.yaz(st.subheader, "m35", "Yöntem karşılaştırması")
    names_ = {"fixed": f"Bölgesel, sabit h = {h_fixed:g}", "sbg": "Bölgesel, SBG", "colburn_hougen": "Bölgesel, Colburn–Hougen"}
    rows = [(f"Tek LMTD, sabit h = {h_fixed:g} (orijinal kod)", h_fixed, R["single"]["A_req"], R["single"]["A_eNTU"],
             R["single"]["over_design"])]
    rows += [(names_[k], rr["h_shell"], rr["A_req"], rr["A_eNTU_cf"] / rr["lf"]["F"], rr["over_design"]) for k, rr in R["methods"]]
    st.dataframe(pd.DataFrame(rows, columns=["Yöntem", "h kabuk (W/m²K)", "A gerekli (m²)", "A ε-NTU (m²)", "Aşırı tasarım (%)"]).round(1),
                 hide_index=True, width="stretch")

    gv.yaz(st.subheader, "m36", "Bölge sayısı kontrolü")
    if R["r2"] is not None:
        d = (R["r2"]["A_req"] / r["A_req"] - 1) * 100
        gv.yaz(st.markdown, "m111", "{n1} bölge: **{A1} m²** · {n2} bölge: **{A2} m²** → fark **%{fark}**. {sonuc}",
               degerler=dict(n1=len(z), A1=f"{r['A_req']:.1f}", n2=len(R["r2"]["zones"]), A2=f"{R['r2']['A_req']:.1f}",
                             fark=f"{d:+.1f}", sonuc="Sonuç bölge sayısından bağımsız ✅" if abs(d) < 3 else "Bölge sayısını artırın ⚠️"))
    else:
        gv.yaz(st.info, "m37", "Tablo girişinde bölge sayısı tablodaki satır sayısıyla sabittir.")

    gv.yaz(st.subheader, "m38", "Duyarlılık analizi")
    st.dataframe(R["sens"].round(2), hide_index=True, width="stretch")
    gv.yaz(st.caption, "m39", "Her satırda tek bir parametre değiştirilmiştir; geometri sabittir.")

    # ---------------- hidrolik doğrulama
    hv = R["hv"]
    gv.yaz(st.subheader, "m100", "Hidrolik doğrulama")
    gv.yaz(st.markdown, "m101", "Tüp sayısı ve basınç düşüşleri, uygulamanın kendi denklemlerinden **bağımsız** kaynaklarla "
           "karşılaştırılır. Amaç iki soruyu ayırmaktır: (1) kod denklemleri doğru uyguluyor mu, (2) yöntem seçimi sonucu ne "
           "kadar değiştiriyor ve tasarım kararı bu değişime dayanıklı mı?")
    gv.yaz(st.markdown, "m102", "**1. Tüp sayısı** — demet çapı Db = {Db} mm için", degerler=dict(Db=f"{sg['Db'] * 1000:.0f}"))
    st.dataframe(hv["tup_sayisi"].round(1), hide_index=True, width="stretch")
    gv.yaz(st.caption, "m103", "Phadke (1984) yöntemi sayı teorisiyle tam sayım verir (ht kütüphanesi). 1 ve 2 geçişte iki sayım "
           "aynı çıkmalıdır. 4 geçiş ve üstünde küçük fark, Phadke'nin geçiş bölme şeritleri için yaklaşık düzeltme kullanmasından "
           "gelir (bu uygulama şerit eksenine Do/2 + 9,5 mm içindeki tüpleri tek tek çıkarır). Towler bağıntısı ampiriktir.")
    gv.yaz(st.markdown, "m104", "**2. Tüp tarafı ΔP**")
    st.dataframe(hv["tup_dp"].round(2), hide_index=True, width="stretch")
    _tp = hv["tup_dp_parca"]
    gv.yaz(st.caption, "m105", "Uygulamanın tüp ΔP'si = sürtünme {surt} kPa + giriş/çıkış/dönüş kayıpları {donus} kPa. "
           "Kern yöntemi dönüş kayıplarını geçiş başına 4 hız yükü alır; Sinnott'un 2,5 hız yükü giriş, çıkış ve dönüşü "
           "birlikte kapsar. Tüp ΔP sınırı nozul kayıpları dahil olarak da kontrol edilmelidir.",
           degerler=dict(surt=f"{_tp['surtunme']:.2f}", donus=f"{_tp['donus']:.2f}"))
    if r.get("shell_method") == "bell_delaware":
        gv.yaz(st.markdown, "m106", "**3. Bell–Delaware düzeltme çarpanları** — kabuk girişi, tek faz gaz (Re = {Re})",
               degerler=dict(Re=f"{hv['bd_giris']['Re']:.0f}"))
        st.dataframe(hv["bd_carpanlar"].round(3), hide_index=True, width="stretch")
        gv.yaz(st.caption, "m107", "Uygulama Taborek'in (HEDH) denklemlerini kullanır; ht kütüphanesi aynı çarpanları Taborek "
               "grafiklerinin sayısallaştırılmasından okur. %5 içinde uyum, grafik okuma hassasiyeti mertebesindedir.")
    gv.yaz(st.markdown, "m108", "**4. Kabuk ΔP — yöntem ve iki faz sabiti duyarlılığı** (sınır {sinir} kPa)",
           degerler=dict(sinir=f"{lim['dp_s']:.0f}"))
    st.dataframe(hv["kabuk_dp"].round(2), hide_index=True, width="stretch")
    _sp = hv["kabuk_dp_parca"]
    gv.yaz(st.caption, "m109", "Tek faz gaz Bell–Delaware ΔP'sinin dağılımı: iç çapraz akış bölmeleri {capraz} kPa, pencereler "
           "{pencere} kPa, giriş/çıkış bölmeleri {uc} kPa. Tasarımın bütün satırlarda sınırın altında kalması, kararın yöntem "
           "seçimine duyarlı olmadığını gösterir. Literatürde Bell–Delaware ΔP belirsizliği ±%30–50 mertebesindedir.",
           degerler=dict(capraz=f"{_sp['B_capraz']:.2f}", pencere=f"{_sp['B_pencere']:.2f}", uc=f"{_sp['B_uc']:.2f}"))

# ---------------------------------------------------------------- MEKANİK
for _ in gv.sekme(T, "mekanik"):
    gv.yaz(st.subheader, "m40", "TEMA tipi ve ısıl genleşme")
    c1, c2, c3 = st.columns(3)
    c1.metric("Kabuk metal T (ort.)", f"{m['T_shell_metal']:.1f} °C")
    c2.metric("Tüp metal T (ort.)", f"{m['T_tube_metal']:.1f} °C")
    c3.metric("Fark / boy farkı", f"{m['dT_metal']:.1f} °C / {m['delta_mm']:.2f} mm")
    st.latex(r"\delta = \alpha\,L\,(T_{\text{kabuk}} - T_{\text{tüp}}), \quad \alpha_{\text{çelik}} = 12\times10^{-6}\ \text{K}^{-1}")
    st.success(f"**Öneri: {m['tema']}** — {m['tema_reason']}")
    gv.yaz(st.caption, "m112", "{not_} Kural: fark < 50 °C ise sabit ayna uygundur.", degerler=dict(not_=m["tema_note"]))

    gv.yaz(st.subheader, "m41", "Tüp yerleşimi (ayna)")
    D_otl_max = geo["Ds"] - R["clearance"]
    c1, c2 = st.columns([1, 1])
    fig, n_fit = draw_tubesheet(geo["Ds"], D_otl_max, geo["Do"], geo["Pt"], geo["layout"], geo["passes"], geo["N"])
    c1.pyplot(fig)
    c2.markdown(f"""
| | Tüp sayısı |
|---|---|
| Girilen N | **{geo['N']}** |
| Bu kabuğa sığan en çok (geometrik sayım) | **{n_fit}** |
| Towler bağıntısı (Db = Ds − boşluk) | {int(sg['K1'] * (D_otl_max / geo['Do']) ** sg['n1'])} |
| Tüp başına geçiş | {geo['N'] / geo['passes']:.0f} |
| Demet çapı (girilen N için) | {geo['Db']:.0f} mm |
""")
    gv.yaz(c2.caption, "m42", "Mavi: kullanılan tüpler (merkezden dışa), gri: boş kalan yerler, turuncu: geçiş bölme şeritleri. "
               "Şerit ekseninden tüp merkezine en az Do/2 + 9,5 mm bırakılır.")
    if geo["N"] > n_fit:
        gv.yaz(c2.error, "m113", "Girilen tüp sayısı bu kabuğa sığmıyor: Ds'yi büyütün veya N'yi azaltın.")

    gv.yaz(st.subheader, "m43", "Nozullar (TEMA ρv²)")
    st.dataframe(m["nozzles"].round(1), hide_index=True, width="stretch")
    gv.yaz(st.markdown, "m114", "**Çarpma plakası:** {karar}", degerler=dict(karar=m["impingement"]))
    gv.yaz(st.caption, "m44", "TEMA RCB-4.61: aşındırıcı olmayan tek faz için ρv² > 2232 kg/m·s² ise çarpma plakası; "
               "doymuş buhar ve iki fazlı girişte her zaman. Tüp nozulu sınırı 8928 kg/m·s².")

    gv.yaz(st.subheader, "m45", "Tüp destek açıklığı")
    gv.yaz(st.markdown, "m115", "Segmental baffle'da pencere tüpleri her iki baffle'da bir desteklenir: açıklık = 2B = "
           "**{aciklik} mm** · TEMA RCB-4.52 sınırı (Ø{Do} karbon çeliği): **{sinir} mm** → {sonuc}",
           degerler=dict(aciklik=f"{m['span_mm']:.0f}", Do=geo["Do"], sinir=f"{m['span_max_mm']:.0f}",
                         sonuc="✅ uygun" if m["span_ok"] else "❌ baffle aralığını azaltın"))

    gv.yaz(st.subheader, "m88", "Akış kaynaklı titreşim (TEMA Bölüm 6)")
    vb = m["vib"]
    c1, c2, c3 = st.columns(3)
    c1.metric("En büyük V/Vc", f"{vb['V_Vc_max']:.2f}", help=vb["kritik_aciklik"])
    c2.metric("En düşük doğal frekans", f"{vb['fn_min']:.1f} Hz")
    c3.metric("Baffle kalınlığı (TEMA R-4.41)", f"{vb['tb_mm']:.1f} mm")
    if not vb["fei_ok"]:
        gv.yaz(st.error, "m116", "Akışkan-elastik kararsızlık riski: bir açıklıkta çapraz akış hızı kritik hızı aşıyor "
               "(TEMA V-10). Baffle aralığını daraltın.")
    elif not vb["margin_ok"]:
        gv.yaz(st.warning, "m117", "V < Vc sağlanıyor, ancak %20 güvenlik payı (V/Vc ≤ 0,8) sağlanmıyor.")
    else:
        gv.yaz(st.success, "m118", "Bütün açıklıklarda çapraz akış hızı kritik hızın %80'inin altında: akışkan-elastik "
               "kararsızlık beklenmez.")
    sp = vb["spans"][["acıklık", "l_mm", "uc_kosulu", "fn_Hz", "V_ms", "Vc_ms", "V_Vc"]].rename(columns={
        "acıklık": "Açıklık", "l_mm": "Uzunluk (mm)", "uc_kosulu": "Uç koşulu", "fn_Hz": "fn (Hz)", "V_ms": "V (m/s)",
        "Vc_ms": "Vc (m/s)", "V_Vc": "V/Vc"})
    st.dataframe(sp.round(2), hide_index=True, width="stretch")
    gv.yaz(st.caption, "m89", "Hesap TEMA Bölüm 6 bağıntılarıyla her tipik desteksiz açıklık için yapılır: doğal frekans (V-5.3), "
           "etkin tüp ağırlığı (V-7), sönüm (V-8), kritik hız Vc = D·fn·do (V-10). Çapraz akış hızında tüm kabuk debisinin "
           "çapraz akış alanından geçtiği kabul edilir (kaçak ve bypass düşülmez, güvenli tarafta). Akustik rezonans: {akustik}.",
           degerler=dict(akustik="olası, gerekirse de-tuning baffle eklenir (V-12.51)" if vb["acoustic_risk"] else "beklenmez"))

    gv.yaz(st.subheader, "m90", "Ayna kalınlığı (TEMA Ek A.13, RCB-7.11)")
    tsh = m["tubesheet"]
    st.latex(rf"T_{{\text{{eğilme}}}} = \frac{{F\,G}}{{3}}\sqrt{{\frac{{P}}{{\eta\,S}}}} = \frac{{{tsh['F']}\times{tsh['G_mm']:.0f}}}{{3}}"
             rf"\sqrt{{\frac{{{tsh['P_kPa']:.0f}}}{{{tsh['eta']:.3f}\times{tsh['S_kPa']:.0f}}}}} = {tsh['T_bend_mm']:.1f}\ \text{{mm}}")
    st.latex(rf"T_{{\text{{kesme}}}} = \frac{{0.31\,D_L}}{{1 - d_o/P_t}}\,\frac{{P}}{{S}} = {tsh['T_shear_mm']:.1f}\ \text{{mm}}, \qquad "
             rf"\eta = 1 - \frac{{0.907}}{{(P_t/d_o)^2}} = {tsh['eta']:.3f}")
    c1, c2, c3 = st.columns(3)
    c1.metric("Hesaplanan + korozyon payı", f"{tsh['T_req_mm']:.1f} mm", help=f"Belirleyici: {tsh['kontrol']}; korozyon payı {tsh['CA_mm']:.0f} mm (iki taraf)")
    c2.metric("TEMA R-7.11 alt sınırı", f"{tsh['T_min_mm']:.1f} mm", help="Mandrelle genişletilen bölgede (T − CA) ≥ tüp dış çapı, toplam ≥ 19,1 mm")
    c3.metric("Seçilen ayna kalınlığı", f"{tsh['T_sel_mm']} mm")
    gv.yaz(st.caption, "m91", "Sabit (kabuğa kaynaklı) ayna için F, Şekil A.131 eğri H'nin üst değeri 1,25 alınmıştır (güvenli tarafta). "
           "G = kabuk iç çapı, P = kabuk ve tüp tasarım basınçlarının büyüğü. Isıl genleşmeden gelen etkin basınç (A.15) ayrıca "
           "hesaplanmamış; genleşme kuvveti yukarıdaki ısıl genleşme kontrolünde sınanmıştır.")

    gv.yaz(st.subheader, "m46", "Et kalınlıkları (ASME VIII Div.1, UG-27)")
    st.latex(rf"t_{{\text{{kabuk}}}} = \frac{{P R}}{{S E - 0.6 P}} = \frac{{{P_d_shell / 10:.2f}\times{geo['Ds'] / 2:.0f}}}{{{S_shell:.0f}\times{E_joint} - 0.6\times{P_d_shell / 10:.2f}}}"
             rf" = {m['t_shell_calc']:.2f}\ \text{{mm}} + CA\ {CA_shell:.1f} = {m['t_shell_req']:.2f}\ \text{{mm}}")
    st.latex(rf"t_{{\text{{tüp}}}} = {m['t_tube_calc']:.3f}\ \text{{mm}} \le {geo['t']}\ \text{{mm}} \Rightarrow " + (r"\text{yeterli}" if m["t_tube_ok"] else r"\text{yetersiz}"))
    gv.yaz(st.caption, "m47", "Basınç düşük olduğundan hesaplanan kalınlık küçüktür; TEMA minimum kabuk kalınlığı (karbon çeliği, bu çap için) "
               "genellikle belirleyicidir — TEMA tablosundan kontrol edin.")

# ---------------------------------------------------------------- MALİYET
for _ in gv.sekme(T, "maliyet"):
    ci = R["cost_in"]
    gv.yaz(st.subheader, "m48", "Yatırım maliyeti")
    st.latex(rf"C_{{e,2010}} = a + b\,S^n = {c['a']:.0f} + {c['b']:g}\times{c['S']:.1f}^{{{c['n']}}} = \${c['Ce_2010']:.0f}")
    st.latex(rf"C_e = C_{{e,2010}}\times\frac{{CEPCI_{{\text{{güncel}}}}}}{{CEPCI_{{2010}}}} = {c['Ce_2010']:.0f}\times\frac{{{ci['cepci_now']:.1f}}}{{{ci['cepci_base']}}} = \${c['Ce']:.0f}")
    f = ci["install"]
    st.latex(rf"C_{{\text{{kurulu}}}} = C_e\left[(1+f_p)f_m + (f_{{er}}+f_{{el}}+f_i+f_c+f_s+f_l)\right] = {c['Ce']:.0f}\times[(1+{f['f_p']})\times{ci['f_m']} + "
             rf"{f['f_er'] + f['f_el'] + f['f_i'] + f['f_c'] + f['f_s'] + f['f_l']:.2f}] = \${c['C_installed']:.0f}")
    gv.yaz(st.caption, "m119", "Bağıntı: {bag} — Towler & Sinnott (2. bs.) Tablo 7.2, ABD Körfez Kıyısı, Ocak 2010, CEPCI = 532,9. "
           "Geçerlilik 10–1000 m²{aralik}. Sabit aynalı tip tabloda yoktur; U-tüp bağıntısı yaklaşık olarak kullanılır.",
           degerler=dict(bag=ci["correlation"], aralik="" if c["in_range"] else " ⚠️ alan aralık dışında!"))

    gv.yaz(st.subheader, "m49", "İşletme maliyeti")
    st.latex(rf"P_{{\text{{pompa}}}} = \frac{{\Delta P_t\,\dot V}}{{\eta}} = \frac{{{r['dp_tube_kPa']:.1f}\times{c['V_w_m3h'] / 3600:.4f}}}{{{ci['pump_eff']}}} = {c['P_pump_kW']:.2f}\ \text{{kW}}")
    st.latex(rf"C_{{\text{{elektrik}}}} = {c['P_pump_kW']:.2f}\times{ci['hours']:.0f}\times{ci['elec_price']} = \${c['E_cost']:.0f}/\text{{yıl}}, \quad "
             rf"C_{{\text{{su}}}} = {c['V_w_m3h']:.0f}\ \text{{m}}^3/\text{{h}}\times{ci['hours']:.0f}\times{ci['water_price_m3']} = \${c['W_cost']:.0f}/\text{{yıl}}")

    gv.yaz(st.subheader, "m50", "Yıllık toplam maliyet")
    st.latex(rf"CRF = \frac{{i(1+i)^n}}{{(1+i)^n - 1}} = {c['crf']:.4f}, \quad TAC = C_{{\text{{kurulu}}}}\,CRF + C_{{\text{{elektrik}}}} + C_{{\text{{su}}}} = \${c['TAC']:.0f}/\text{{yıl}}")
    st.dataframe(pd.DataFrame([
        ("Satın alma (güncel)", money(c["Ce"])), ("Kurulu yatırım", money(c["C_installed"])),
        ("Yıllık sermaye payı", money(c["annual_capital"])), ("Pompa elektriği / yıl", money(c["E_cost"])),
        ("Soğutma suyu / yıl", money(c["W_cost"])), ("Yıllık toplam maliyet", money(c["TAC"]))],
        columns=["Kalem", "Tutar"]), hide_index=True, width="stretch")
    gv.yaz(st.caption, "m51", "Elektrik ve su fiyatları örnek değerlerdir; kendi değerlerinizi girin.")

# ---------------------------------------------------------------- DWSIM
for _ in gv.sekme(T, "dwsim"):
    gv.yaz(st.markdown, "m60", """
DWSIM ile iki katmanda karşılaştırın:

| Katman | Ne karşılaştırılır | Beklenti |
|---|---|---|
| **1. Termodinamik** | Q, çiy noktası, çıkış gaz ve sıvı debisi, gazdaki aseton | Aynı NRTL parametreleriyle **çok yakın** olmalı (±%3, ±2 °C). Asıl doğrulama budur. |
| **2. Isı transferi ve hidrolik** | U, tüp ΔP, kabuk ΔP | DWSIM kabuk tarafını basitleştirilmiş **Tinker** yöntemiyle ve yoğuşmayı ayrıca modellemeden hesaplar; bu uygulama **Bell–Delaware + Nusselt + Colburn–Hougen** ve iki faz çarpanı kullanır. Fark **beklenir**; hidroliğin bağımsız doğrulaması **Bölgeler & doğrulama** sekmesindedir. |
""")
    adaylar_ = dwsim_adaylari()
    dwsim_girdi = dwsim_girdisi(R, adaylar_)
    with st.expander("DWSIM'de elle çalıştırma (Python betiği)"):
        gv.yaz(st.markdown, "m74", """
Uygulama klasöründeki **dwsim_dogrulama.py** betiği, aşağıdaki dosyadaki değerlerle DWSIM'de akış şemasını
kendisi kurar, çözer ve sonuçları karşılaştırır. Yalnızca DWSIM kurulu bir Windows bilgisayarda çalışır.

1. Bir kez: `pip install pythonnet pywin32`
2. Aşağıdaki **dwsim_girdi.json** dosyasını indirip uygulama klasörüne koyun.
3. Komut satırında: `python dwsim_dogrulama.py dwsim_girdi.json`
4. Ekrana çıkan tablodaki DWSIM değerlerini aşağıdaki kutulara girin. Betik ayrıca DWSIM'de açılabilen
   **E-102_dwsim.dwxmz** dosyasını oluşturur.
""")
        st.download_button("dwsim_girdi.json indir", json.dumps(dwsim_girdi, ensure_ascii=False, indent=2),
                           file_name="dwsim_girdi.json", mime="application/json")
    gv.yaz(st.subheader, "m71", "DWSIM ile otomatik karşılaştırma")
    if _platform.system() != "Windows":
        gv.yaz(st.info, "m72", "Otomatik karşılaştırma yalnızca DWSIM kurulu bir Windows bilgisayarda çalışır. Başka bir "
               "bilgisayarda yukarıdaki betiği kullanın veya DWSIM değerlerini aşağıya elle girin.")
    elif (dw_klasor := dwv.dwsim_klasoru()) is None:
        gv.yaz(st.warning, "m73", "Bu bilgisayarda DWSIM bulunamadı. DWSIM'i kurun ya da kurulum klasörünü DWSIM_BIN "
                                  "ortam değişkeniyle belirtip uygulamayı yeniden başlatın.")
    else:
        gv.yaz(st.caption, "m121", "DWSIM bulundu: {klasor}", degerler=dict(klasor=dw_klasor))
        dw_coklu = st.checkbox(f"Tasarım taramasındaki en uygun {0 if adaylar_ is None else len(adaylar_)} tasarımı da dene "
                               "(daha uzun sürer)", value=adaylar_ is not None, disabled=adaylar_ is None)
        if st.button("DWSIM ile karşılaştır", type="primary"):
            girdi = dict(dwsim_girdi)
            if not dw_coklu:
                girdi.pop("tasarimlar", None)
            with st.spinner("DWSIM çalışıyor, 1–3 dakika sürebilir…"):
                ss["dwsim_otomatik"] = dict(dwsim_calistir(girdi), adaylar=adaylar_ if dw_coklu else None)
            sat = ss["dwsim_otomatik"]["satirlar"]
            if sat:                                    # elle giriş kutularını da doldur
                def _bul(test, ad):
                    return next((x["dwsim"] for x in sat if x["test"] == test and x["buyukluk"] == ad), None)
                ilk = next((x["test"] for x in sat if x["test"] not in ("Termodinamik", "U·A verilerek")), None)
                for key, test, ad in (("dw_Q", "Termodinamik", "Isı yükü"), ("dw_Tdew", "Termodinamik", "Çiy noktası"),
                                      ("dw_mV", "Termodinamik", "Çıkış gaz debisi"), ("dw_prod", "Termodinamik", "Gazdaki aseton"),
                                      ("dw_U", ilk, "U (kirli)"), ("dw_dpt", ilk, "Tüp tarafı ΔP"), ("dw_dps", ilk, "Kabuk ΔP")):
                    v_ = _bul(test, ad)
                    if v_ is not None and v_ == v_:
                        ss[key] = float(v_)
            st.rerun()

        sonuc_dw = ss.get("dwsim_otomatik")
        if sonuc_dw:
            sat = sonuc_dw["satirlar"]
            if not sat:
                gv.yaz(st.error, "m122", "DWSIM karşılaştırması tamamlanamadı. Ayrıntılar aşağıdaki kayıtta.")
            else:
                df_dw = pd.DataFrame(sat)
                ilk = next((t_ for t_ in df_dw.test if t_ not in ("Termodinamik", "U·A verilerek")), None)

                def _ort(maske):
                    v_ = df_dw.loc[maske, "uyum"].dropna()
                    return None if v_.empty else float(v_.mean())
                termo = _ort(df_dw.test.eq("Termodinamik") & df_dw.buyukluk.isin(
                    ["Isı yükü", "Çiy noktası", "Çıkış gaz debisi", "Çıkış sıvı debisi"]))
                esj = _ort(df_dw.test.eq(ilk) & df_dw.buyukluk.eq("U (kirli)"))
                sinir = df_dw[(df_dw.test == ilk) & (df_dw.buyukluk == "Kabuk ΔP (sınır)")]
                mc_ = st.columns(3)
                mc_[0].metric("Termodinamik uyum", "–" if termo is None else f"%{termo:.1f}",
                              help="Isı yükü, çiy noktası, çıkış gaz ve sıvı debilerinin ortalama benzerliği (100 − bağıl fark).")
                mc_[1].metric("U uyumu (bilgi)", "–" if esj is None else f"%{esj:.1f}",
                              help="İki program farklı kabuk tarafı yöntemleri kullandığı için tam uyum beklenmez.")
                if not sinir.empty:
                    mc_[2].metric("DWSIM'e göre kabuk ΔP", f"{sinir.dwsim.iloc[0]:.1f} kPa",
                                  delta=sinir.sonuc.iloc[0], delta_color="normal" if "İÇİNDE" in sinir.sonuc.iloc[0] else "inverse")
                if termo is not None and termo >= 95:
                    gv.yaz(st.success, "m123", "Termodinamik hesap DWSIM ile uyumlu.")
                elif termo is not None:
                    gv.yaz(st.warning, "m124", "Termodinamik hesapta beklenenden büyük fark var; bileşenleri ve NRTL "
                           "parametrelerini kontrol edin.")
                if sonuc_dw.get("adaylar") is not None:
                    gv.yaz(st.markdown, "m125", "**Tasarımların DWSIM'deki sonucu**")
                    tb_ = dwsim_adaylari_degerlendir(sat, sonuc_dw["adaylar"], R)
                    st.dataframe(tb_.drop(columns="gecti"), hide_index=True, width="stretch")
                    gv.yaz(st.caption, "m78", "Geçme koşulu: DWSIM'e göre ısı yükünün en az %98'i aktarılmalı ve kabuk/tüp "
                           "basınç düşüşleri izin verilen sınırların altında kalmalı. Bir tasarım geçmezse sıradaki "
                           "tasarımı **Tasarım taraması** sekmesinden girdilere aktarabilirsiniz.")
                st.dataframe(df_dw.rename(columns={"test": "Test", "buyukluk": "Büyüklük", "birim": "Birim",
                                                   "uygulama": "Uygulama", "dwsim": "DWSIM", "fark": "Fark",
                                                   "sonuc": "Sonuç", "uyum": "Uyum (%)"}),
                             hide_index=True, width="stretch")
                st.download_button("Sonuçları indir (CSV)", df_dw.to_csv(index=False, sep=";").encode("utf-8-sig"),
                                   file_name="dwsim_sonuclar.csv", mime="text/csv")
            with st.expander("DWSIM çalışma kaydı"):
                st.code(sonuc_dw["log"][-6000:] or "(boş)")

    gv.yaz(st.subheader, "m61", "DWSIM sonuçları (elle giriş)")
    app_vals = app_values(R)
    for _k in ("dw_Q", "dw_Tdew", "dw_mV", "dw_prod", "dw_U", "dw_dpt", "dw_dps"):
        ss.setdefault(_k, 0.0)
    cc = st.columns(4)
    dw = dict(Q=cc[0].number_input("Q (kW)", key="dw_Q"),
              T_dew=cc[1].number_input("Çiy noktası (°C)", key="dw_Tdew"),
              m_V_out=cc[2].number_input("Çıkış buharı (kg/h)", key="dw_mV"),
              product_V_out=cc[3].number_input(f"Buhardaki {product or 'ürün'} (kg/h)", key="dw_prod"))
    cc = st.columns(4)
    dw.update(U=cc[0].number_input("U (W/m²K)", key="dw_U"),
              dp_tube=cc[1].number_input("Tüp ΔP (kPa)", key="dw_dpt"),
              dp_shell=cc[2].number_input("Kabuk ΔP (kPa)", key="dw_dps"))
    cmp_df = ec.dwsim_compare(app_vals, dw, ec.dwsim_bounds(r, R["hv"]))
    if cmp_df.empty:
        gv.yaz(st.info, "m62", "DWSIM değerlerini girdiğinizde karşılaştırma tablosu burada oluşur. Boş (0) bırakılan satırlar atlanır.")
        st.dataframe(pd.DataFrame([{"Büyüklük": n, "Birim": u, "Uygulama": app_vals[k]}
                                   for k, n, u, *_ in ec.DWSIM_ITEMS if app_vals.get(k) is not None]).round(2),
                     hide_index=True, width="stretch")
    else:
        st.dataframe(cmp_df.round(2), hide_index=True, width="stretch")
        gv.yaz(st.caption, "m139", "Termodinamik satırlar toleransla karşılaştırılır. U ve ΔP satırlarında DWSIM yoğuşmayı ve "
               "H₂ direncini ayrıca modellemediği için tek bir değer beklenmez; DWSIM sonucu uygulamanın iki uç durumu "
               "arasında (±%10) kalmalıdır. U: yoğuşmasız gaz ↔ yoğuşmalı · tüp ΔP: yalnız sürtünme ↔ Kern (4 hız yükü) · "
               "kabuk ΔP: iki fazlı ↔ gaz hiç yoğuşmuyor.")
        thermo_ok = all(cmp_df[cmp_df["Büyüklük"].isin(["Isı yükü", "Çiy noktası", "Çıkış buhar debisi"])]["Durum"].str.startswith("✅"))
        if thermo_ok:
            gv.yaz(st.success, "m63", "Termodinamik katman uyumlu: yoğuşma eğrisi DWSIM ile tutarlı.")
        else:
            gv.yaz(st.warning, "m64", "Termodinamik katmanda fark var: iki programda aynı paketi seçtiğinizden emin olun "
                       "ya da DWSIM eğrisini tablo olarak içe aktarın.")

# ---------------------------------------------------------------- TASARIM TARAMASI (OPTİMİZASYON)
for _ in gv.sekme(T, "optimizasyon"):
    gv.yaz(st.markdown, "m65", """
**Amaç:** tüm kriterleri sağlayan tasarımlar arasında **yıllık toplam maliyeti (TAC)** en düşük olanı bulmak.

**Karar değişkenleri:** soğutma suyu çıkış sıcaklığı + geometri (tüp çapı, boy, geçiş sayısı, kabuk çapı, baffle aralığı).
Her su sıcaklığı için tüm standart geometriler denenir; her birinde bölgesel ısıl hesap, tüp sayısı (ayna sayımı) ve
basınç düşüşleri yeniden hesaplanır. Su çıkış sıcaklığı yükseldikçe su debisi azalır (**işletme maliyeti düşer**) ama
sürücü sıcaklık farkı ve F küçülür (**alan ve yatırım artar**). En uygun nokta bu ikisinin dengesidir.
""")
    st.latex(r"TAC = \underbrace{C_{kurulu}\cdot CRF}_{\text{sabit (yıllık sermaye)}} + "
             r"\underbrace{C_{elektrik} + C_{su}}_{\text{işletme}}")
    with st.form("opt"):
        c1, c2, c3, c4 = st.columns(4)
        T_rng = c1.slider("Su çıkış sıcaklığı aralığı (°C)", 27, 50, (30, 44))
        T_step = c2.selectbox("Adım (°C)", [1, 2], index=0)
        n_top = c3.slider("Tabloda model sayısı", 3, 10, 5)
        h_rob = c4.number_input("Sağlamlık testi: kabuk h ×", value=0.8, step=0.05)
        sizes = st.multiselect("Tüp çapları", [19.05, 25.4], default=[19.05, 25.4])
        go_opt = st.form_submit_button("🏆 Taramayı çalıştır (her sıcaklık ~2–3 s)")
    if go_opt:
        T_list = list(range(T_rng[0], T_rng[1] + 1, T_step))
        bar = st.progress(0.0, text="Başlıyor…")
        allc, best = ec.optimize(R["tube"], curve, T_list, cost_for, R["kw"], lim, clearance=R["clearance"],
                                 tube_sizes=tuple((d, TUBE_SIZES[d]) for d in sizes),
                                 progress=lambda f, T: bar.progress(f, text=f"Su çıkışı {T} °C tarandı"))
        bar.empty()
        top = allc.head(n_top).copy() if not allc.empty else allc
        if not top.empty:
            top["asiri_tasarim_h_dusuk"] = [ec.robustness(R["tube"], curve, row, R["kw"], h_rob) for _, row in top.iterrows()]
        ss["opt_res"] = dict(allc=allc, best=best, top=top, T_list=T_list, h_rob=h_rob, fp=opt_fingerprint(R))

    if "opt_res" in ss:
        O = ss["opt_res"]
        allc, best, top = O["allc"], O["best"], O["top"]
        if allc.empty:
            gv.yaz(st.warning, "m66", "Kriterleri sağlayan tasarım bulunamadı; sınırları veya sıcaklık aralığını değiştirin.")
        else:
            opt = allc.sort_values("yillik_toplam_USD").iloc[0]
            mc = st.columns(4)
            mc[0].metric("Denenen uygun tasarım", f"{len(allc)}")
            mc[1].metric("En düşük TAC", f"${fmt(opt.yillik_toplam_USD, 0)}/yıl", help=f"Su çıkışı {opt.T_su_cikis:.0f} °C")
            mc[2].metric("Seçilen tasarımın TAC'ı", f"${fmt(c['TAC'], 0)}/yıl")
            mc[3].metric("Seçilen: su çıkışı / alan", f"{w['T_out']:.0f} °C / {r['A']:.0f} m²")
            gv.yaz(st.caption, "m67", "Grafikler her hesaptan sonra seçilen tasarıma göre otomatik yeniden çizilir.")
            figs, pngs = make_figs(R, O)
            yazilar = {
                "1": (lambda: gv.yaz(st.subheader, "m126", "Grafik 1 — TAC – tüp sayısı"),
                      lambda: gv.yaz(st.caption, "m127", "Seçilen tasarımın çapı, geçişi ve su sıcaklığıyla tüp sayısı "
                                     "değiştirilmiştir (◆ seçilen). Kesikli kısım kriterleri sağlamaz.")),
                "2": (lambda: gv.yaz(st.subheader, "m128", "Grafik 2 — Sıcaklık – ısı yükü profili"),
                      lambda: gv.yaz(st.caption, "m129", "Her bölgede proses, su ve arayüz sıcaklığı. Gri alan sürücü sıcaklık farkıdır.")),
                "3": (lambda: gv.yaz(st.subheader, "m130", "Grafik 3 — Yıllık toplam maliyet"),
                      lambda: gv.yaz(st.caption, "m131", "Her su çıkış sıcaklığının en ucuz uygun geometrisi; ◆ seçilen tasarım.")),
                "4": (lambda: gv.yaz(st.subheader, "m132", "Grafik 4 — Yatırım – işletme dengesi"),
                      lambda: gv.yaz(st.caption, "m133", "Yatırım (yıllık sermaye payı) ile işletme maliyetinin (elektrik + su) "
                                     "karşılıklı değişimi.")),
            }
            for k, f in figs.items():
                yazilar[k[0]][0]()
                st.pyplot(f)
                yazilar[k[0]][1]()
            import zipfile
            zbuf = io.BytesIO()
            with zipfile.ZipFile(zbuf, "w") as zf:
                for k, b_ in pngs.items():
                    zf.writestr(f"{k}.png", b_)
            st.download_button("⬇ Grafikleri indir (PNG, rapor için)", zbuf.getvalue(), "grafikler.zip")

            gv.yaz(st.subheader, "m68", "En uygun modeller")
            cols_ = ["T_su_cikis", "Do", "L", "gecis", "Ds", "B", "N", "A", "asiri_tasarim", "asiri_tasarim_h_dusuk",
                     "v_tup", "dP_tup", "dP_kabuk", "yatirim_USD", "yillik_sermaye_USD", "elektrik_USD", "su_USD", "yillik_toplam_USD"]
            show = top[[c_ for c_ in cols_ if c_ in top]].copy()
            show.insert(0, "Sıra", range(1, len(show) + 1))
            if "asiri_tasarim_h_dusuk" in show:
                show["Sağlam"] = np.where(show.asiri_tasarim_h_dusuk >= 0, "✅", "⚠️")
            st.dataframe(show.round(1).assign(Do=show.Do.astype(str)).rename(columns={"asiri_tasarim_h_dusuk": f"aşırı tas. (h×{O['h_rob']:g})"}),
                         hide_index=True, width="stretch")
            gv.yaz(st.caption, "m134", "Sağlamlık: kabuk h'si ×{h} çıkarsa aşırı tasarım hâlâ ≥ 0 mı? H₂'li yoğuşmada kabuk "
                   "h'si en belirsiz büyüklüktür.", degerler=dict(h=f"{O['h_rob']:g}"))

            robust = top[top.asiri_tasarim_h_dusuk >= 0] if "asiri_tasarim_h_dusuk" in top else top.iloc[0:0]
            if robust.empty:
                gv.yaz(st.warning, "m69", "İlk modellerin hiçbiri sağlamlık testini geçmiyor; aşırı tasarım alt sınırını (%10) artırmayı düşünün.")
            else:
                rb = robust.iloc[0]
                if rb.name == top.index[0]:
                    gv.yaz(st.success, "m135", "Önerilen tasarım: **#1** — hem en ucuz hem sağlam ({tac}/yıl).",
                           degerler=dict(tac=f"${fmt(rb.yillik_toplam_USD, 0)}"))
                else:
                    gv.yaz(st.info, "m136", "En ucuz model (#1) sağlamlık testini geçmiyor. Sağlam olan en ucuz model **#{sira}**: "
                           "{tac}/yıl (+%{fark}).", degerler=dict(sira=list(top.index).index(rb.name) + 1,
                                                                  tac=f"${fmt(rb.yillik_toplam_USD, 0)}",
                                                                  fark=f"{(rb.yillik_toplam_USD / top.iloc[0].yillik_toplam_USD - 1) * 100:.1f}"))

            pick = st.selectbox("Girdilere aktarılacak model", list(range(len(top))),
                                format_func=lambda i: f"#{i + 1}: Ø{top.iloc[i].Do} · {top.iloc[i].L:.0f} mm · "
                                                      f"{int(top.iloc[i].gecis)} geçiş · Ds {top.iloc[i].Ds:.0f} · "
                                                      f"su {top.iloc[i].T_su_cikis:.0f} °C · ${top.iloc[i].yillik_toplam_USD:,.0f}/yıl")

            def apply_opt(i=pick, top=top):
                row = top.iloc[i]
                upd = dict(geo_Do=float(row.Do), geo_t=TUBE_SIZES[float(row.Do)], geo_N=int(row.N),
                           geo_passes=int(row.gecis), geo_L=int(row.L), geo_B=int(row.B), geo_Ds=int(row.Ds),
                           geo_Pt=round(1.25 * float(row.Do), 2), geo_layout="triangular", w_T_out=float(row.T_su_cikis))
                if "T_su_giris" in row:
                    upd["w_T_in"] = float(row.T_su_giris)
                ss.update(upd)
            st.button("Bu modeli girdilere aktar (geometri + su sıcaklıkları)", on_click=apply_opt)
            gv.yaz(st.caption, "m70", "Aktardıktan sonra soldan **▶ Hesapla**'ya basın; tüm sekmeler bu tasarıma göre güncellenir.")
