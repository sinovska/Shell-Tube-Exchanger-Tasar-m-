"""
Kabuk–tüp eşanjör tasarımı: kısmi yoğuşturucu (kabukta yoğuşan buhar + H2, tüpte soğutma suyu).
Hesap motoru — arayüzden bağımsızdır; app.py bu modülü kullanır.

Birimler (girdi/çıktı): debi kg/h · sıcaklık °C · basınç bar · ρ kg/m³ · Cp kJ/kg·K · μ Pa·s
k W/m·K · kirlenme m²K/W · geometri mm · ΔP kPa · ısı yükü kW · alan m² · maliyet USD
"""
import math
import numpy as np
import pandas as pd
import ht
import fluids
from scipy.optimize import brentq

G = 9.80665          # m/s²
R_GAS = 8.314        # J/mol·K

# =====================================================================
# 1. GEOMETRİ, ISI DENGESİ, TÜP TARAFI
# =====================================================================
SQUARE_LAYOUTS = ("square", "rotated square")
# Towler & Sinnott demet çapı katsayıları: Nt = K1·(Db/Do)^n1
K1N1 = {"triangular": {1: (0.319, 2.142), 2: (0.249, 2.207), 4: (0.175, 2.285), 6: (0.0743, 2.499), 8: (0.0365, 2.675)},
        "square":     {1: (0.215, 2.207), 2: (0.156, 2.291), 4: (0.158, 2.263), 6: (0.0402, 2.617), 8: (0.0331, 2.643)}}


def equivalent_diameter(layout, Do, Pt):
    """Kabuk tarafı eşdeğer çap [m] (Kern)."""
    if layout in SQUARE_LAYOUTS:
        return 4 * (Pt**2 - math.pi * Do**2 / 4) / (math.pi * Do)
    return 8 * (0.43301 * Pt**2 - math.pi * Do**2 / 8) / (math.pi * Do)


def shell_geometry(geo):
    """Kabuk tarafı geometrisi (SI)."""
    Do, Pt, B, Ds, L = (geo[k] / 1000 for k in ("Do", "Pt", "B", "Ds", "L"))
    K1, n1 = K1N1["triangular" if geo["layout"] == "triangular" else "square"][geo["passes"]]
    Db = geo["Db"] / 1000 if geo.get("Db") else Do * (geo["N"] / K1) ** (1 / n1)   # sayımdan veya Towler'dan
    return dict(Do=Do, Pt=Pt, B=B, Ds=Ds, L=L, Db=Db, K1=K1, n1=n1,
                De=equivalent_diameter(geo["layout"], Do, Pt),
                As=Ds * B * (Pt - Do) / Pt,
                Nr=max((2 / 3) * Db / Pt, 1.0),
                n_baffle=max(int(L / B - 1), 0),
                A=math.pi * Do * L * geo["N"])


def lmtd_F(Th1, Th2, Tc1, Tc2, shells=1):
    if not (Th1 > Tc2 and Th2 > Tc1):
        raise ValueError("Sıcaklık kesişimi: sıcak çıkış > soğuk giriş ve sıcak giriş > soğuk çıkış olmalı.")
    return dict(LMTD=ht.LMTD(Th1, Th2, Tc1, Tc2, counterflow=True),
                F=ht.F_LMTD_Fakheri(Th1, Th2, Tc1, Tc2, shells=shells),
                R=(Th1 - Th2) / (Tc2 - Tc1), P=(Tc2 - Tc1) / (Th1 - Tc1))


def water_side_balance(tube, Q_kW):
    """Soğutma suyunda tek bilinmeyeni (debi veya bir sıcaklık) Q'dan bulur."""
    w = dict(tube)
    missing = [k for k in ("m", "T_in", "T_out") if w.get(k) is None]
    if len(missing) > 1:
        raise ValueError("Soğutma suyunda debi ve iki sıcaklıktan en fazla biri boş bırakılabilir.")
    if w.get("m") is None:
        w["m"] = Q_kW / (w["cp"] * (w["T_out"] - w["T_in"])) * 3600
    elif w.get("T_out") is None:
        w["T_out"] = w["T_in"] + Q_kW / (w["m"] / 3600 * w["cp"])
    elif w.get("T_in") is None:
        w["T_in"] = w["T_out"] - Q_kW / (w["m"] / 3600 * w["cp"])
    return w


def tube_side(tube, geo, shells=1, K_pass=2.5):
    """Tüp içi hız, Re, Pr, Nu, h ve ΔP (Darcy sürtünme + geçiş başına K_pass hız yükü)."""
    N, Np = geo["N"], geo["passes"]
    Di, L = (geo["Do"] - 2 * geo["t"]) / 1000, geo["L"] / 1000
    A_flow = math.pi / 4 * Di**2 * N / Np
    v = tube["m"] / 3600 / (tube["rho"] * A_flow)
    Re = tube["rho"] * v * Di / tube["mu"]
    Pr = tube["cp"] * 1000 * tube["mu"] / tube["k"]
    eD = 0.000045 / Di
    fd = fluids.friction_factor(Re, eD=eD)                         # Darcy, Colebrook
    if Re > 2300:                                                  # Gnielinski (2300 < Re < 5·10⁶)
        Nu = (fd / 8) * (Re - 1000) * Pr / (1 + 12.7 * (fd / 8) ** 0.5 * (Pr ** (2 / 3) - 1))
    else:                                                          # laminer, sabit duvar sıcaklığı
        Nu = 3.66
    dp = shells * Np * (fd * L / Di + K_pass) * tube["rho"] * v**2 / 2
    return dict(Di=Di, A_flow=A_flow, v=v, Re=Re, Pr=Pr, fd=fd, Nu=Nu, h=Nu * tube["k"] / Di,
                dp_kPa=dp / 1000, K_pass=K_pass)


def overall_U(h_shell, h_tube, Rf_shell, Rf_tube, Do, Di, k_wall):
    """Dış yüzeye göre temiz ve kirli U [W/m²K]."""
    dr = Do / Di
    R_clean = 1 / h_shell + Do * math.log(dr) / (2 * k_wall) + dr / h_tube
    return 1 / R_clean, 1 / (R_clean + Rf_shell + dr * Rf_tube)


def kern_dp(m, rho, mu, sg):
    """Kern kabuk ΔP [Pa], bir kabuk. m kg/s."""
    return ht.conv_tube_bank.dP_Kern(m, rho, mu, sg["Ds"], sg["B"], sg["Pt"], sg["Do"], sg["n_baffle"], mu_w=None)


# ---------------------------------------------------------------------
# Bell–Delaware (Taborek, HEDH) — tek fazlı kabuk tarafı h ve ΔP; kaçak ve bypass düzeltmeleri
# ---------------------------------------------------------------------
# İdeal demet j ve f katsayıları (Taborek): Re aralığı alt sınırı -> (a1, a2, b1, b2); ayrıca (a3, a4, b3, b4)
BD_COEF = {
    30: ([(1e4, .321, -.388, .372, -.123), (1e3, .321, -.388, .486, -.152), (1e2, .593, -.477, 4.570, -.476),
          (10, 1.360, -.657, 45.10, -.973), (0, 1.400, -.667, 48.00, -1.000)], (1.450, .519, 7.00, .500)),
    45: ([(1e4, .370, -.396, .303, -.126), (1e3, .370, -.396, .333, -.136), (1e2, .730, -.500, 3.500, -.476),
          (10, .498, -.656, 26.20, -.913), (0, 1.550, -.667, 32.00, -1.000)], (1.930, .500, 6.59, .520)),
    90: ([(1e4, .370, -.395, .391, -.148), (1e3, .107, -.266, .0815, .022), (1e2, .408, -.460, 6.0900, -.602),
          (10, .900, -.631, 32.10, -.963), (0, .970, -.667, 35.00, -1.000)], (1.187, .370, 6.30, .378)),
}
LAYOUT_ANGLE = {"triangular": 30, "rotated square": 45, "square": 90}
# TEMA RCB-4.3: kabuk iç çapı − baffle dış çapı (mm), iç çap üst sınırına göre
TEMA_BAFFLE_CLEARANCE = [(432, 3.2), (991, 4.8), (1372, 6.4), (1753, 7.9), (2134, 9.5), (2540, 11.1)]
TWO_PHASE_C = 8.0   # Ishihara, Palen & Taborek (1980): kabuk tarafı çapraz akış için Chisholm sabiti


def bd_geometry(geo, Bc=0.25, n_ss=None, lane_clear_mm=19.0):
    """Bell–Delaware akış alanları ve sayıları (SI). geo['Db'] = dış tüp sınırı (OTL), mm.
    Bc: baffle kesim oranı. n_ss: sızdırmazlık şeridi çifti (None: demet boşluğu > 20 mm ise her 6 sıraya bir).
    lane_clear_mm: geçiş bölmesi şeridinin tüpler arası açıklığı (akışa paralel şerit bypass'ı için)."""
    Do, Pt, B, Ds, L = (geo[k] / 1000 for k in ("Do", "Pt", "B", "Ds", "L"))
    N, Np = geo["N"], geo["passes"]
    ang = LAYOUT_ANGLE.get(geo["layout"], 30)
    Dotl = (geo.get("Db") or otl_for_tubes(N, geo["Do"], geo["Pt"], geo["layout"], Np)) / 1000
    Lbb = Ds - Dotl
    Dctl = Dotl - Do
    Lpp = {30: 0.866 * Pt, 45: 0.707 * Pt, 90: Pt}[ang]                 # akış yönünde sıra aralığı
    Lpe = {30: Pt, 45: 0.707 * Pt, 90: Pt}[ang]                         # akışa dik etkin adım
    Sm = B * (Lbb + Dctl / Lpe * (Pt - Do))                              # çapraz akış alanı (merkez çizgisi)
    th_ds = 2 * math.acos(1 - 2 * Bc)
    th_ctl = 2 * math.acos(max(-1.0, min(1.0, Ds * (1 - 2 * Bc) / Dctl)))
    Fw = th_ctl / (2 * math.pi) - math.sin(th_ctl) / (2 * math.pi)       # pencerelerdeki tüp oranı
    Fc = 1 - 2 * Fw
    Ntcc = Ds * (1 - 2 * Bc) / Lpp                                       # çapraz akışta geçilen sıra
    Ntcw = 0.8 / Lpp * (Ds * Bc - (Ds - Dctl) / 2)                       # pencerede etkin sıra
    Swg = Ds**2 / 8 * (th_ds - math.sin(th_ds))
    Sw = Swg - N * Fw * math.pi * Do**2 / 4                              # pencere akış alanı
    Dw = 4 * Sw / (math.pi * Do * N * Fw + th_ds * Ds)
    Lsb = next(c for d, c in TEMA_BAFFLE_CLEARANCE if Ds * 1000 <= d) / 1000
    span = 2 * B
    Ltb = (0.4 if (span > 0.914 and geo["Do"] <= 31.8) else 0.8) / 1000  # TEMA RCB-4.2
    Ssb = math.pi * Ds * Lsb / 2 * (1 - th_ds / (2 * math.pi))           # kabuk–baffle kaçağı
    Stb = math.pi / 4 * ((Do + Ltb) ** 2 - Do**2) * N * (1 - Fw)         # tüp–baffle deliği kaçağı
    # 2 geçişli, dikey kesimde (yan yana akış) yatay geçiş şeridi akışa paraleldir → bypass
    Lpl = 0.5 * lane_clear_mm / 1000 if Np >= 2 else 0.0
    Sb = B * (Lbb + Lpl)
    Nb = max(int(L / B - 1), 1)
    B_end = (L - (Nb - 1) * B) / 2                                       # giriş/çıkış bölmesi aralığı
    if n_ss is None:
        n_ss = 0 if Lbb * 1000 <= 20 else math.ceil(Ntcc / 6)
    return dict(Do=Do, Pt=Pt, B=B, Ds=Ds, L=L, N=N, ang=ang, Dotl=Dotl, Lbb=Lbb, Dctl=Dctl, Lpp=Lpp, Sm=Sm,
                Fw=Fw, Fc=Fc, Ntcc=Ntcc, Ntcw=Ntcw, Sw=Sw, Dw=Dw, Lsb=Lsb, Ltb=Ltb, Ssb=Ssb, Stb=Stb, Lpl=Lpl,
                Sb=Sb, Fsbp=Sb / Sm, Nb=Nb, B_in=B_end, B_out=B_end, n_ss=n_ss, rss=n_ss / Ntcc, th_ds=th_ds, Bc=Bc)


def _bd_jf(Re, ang, pt_do):
    rows, (a3, a4, b3, b4) = BD_COEF[ang]
    _, a1, a2, b1, b2 = next(r for r in rows if Re >= r[0])
    a = a3 / (1 + 0.14 * Re**a4)
    b = b3 / (1 + 0.14 * Re**b4)
    return a1 * (1.33 / pt_do) ** a * Re**a2, b1 * (1.33 / pt_do) ** b * Re**b2


def bell_delaware(m, rho, mu, cp_kJ, k, bd):
    """Tek fazlı akış, Bell–Delaware. m kg/s. Döndürür: h (W/m²K), ΔP bileşenleri (Pa), J ve R çarpanları."""
    Re = bd["Do"] * m / (mu * bd["Sm"])
    Pr = cp_kJ * 1000 * mu / k
    G = m / bd["Sm"]
    j, f = _bd_jf(Re, bd["ang"], bd["Pt"] / bd["Do"])
    h_id = j * cp_kJ * 1000 * G * Pr ** (-2 / 3)
    lam = Re < 100
    Jc = 0.55 + 0.72 * bd["Fc"]
    rs = bd["Ssb"] / (bd["Ssb"] + bd["Stb"])
    rlm = (bd["Ssb"] + bd["Stb"]) / bd["Sm"]
    Jl = 0.44 * (1 - rs) + (1 - 0.44 * (1 - rs)) * math.exp(-2.2 * rlm)
    bypass = 1 - (2 * bd["rss"]) ** (1 / 3) if bd["rss"] < 0.5 else 0.0
    Jb = math.exp(-(1.25 if lam else 1.35) * bd["Fsbp"] * bypass)
    n = 1 / 3 if lam else 0.6
    Lin, Lout, Nb = bd["B_in"] / bd["B"], bd["B_out"] / bd["B"], bd["Nb"]
    Js = (Nb - 1 + Lin ** (1 - n) + Lout ** (1 - n)) / (Nb - 1 + Lin + Lout)
    Nc = (bd["Ntcc"] + bd["Ntcw"]) * (Nb + 1)
    Jr = 1.0 if Re >= 100 else ((10 / Nc) ** 0.18 if Re <= 20 else
                                (10 / Nc) ** 0.18 + (Re - 20) / 80 * (1 - (10 / Nc) ** 0.18))
    h = h_id * Jc * Jl * Jb * Js * Jr
    # basınç düşüşü
    dPbi = 2 * f * bd["Ntcc"] * G**2 / rho                                # ideal bir çapraz akış bölmesi
    p = -0.15 * (1 + rs) + 0.8
    Rl = math.exp(-1.33 * (1 + rs) * rlm**p)
    Rb = math.exp(-(4.5 if lam else 3.7) * bd["Fsbp"] * bypass)
    n2 = 1.0 if lam else 0.2
    Rs = (1 / Lin) ** (2 - n2) + (1 / Lout) ** (2 - n2)
    mw = m / math.sqrt(bd["Sm"] * bd["Sw"])
    if lam:
        dPwi = 26 * mu * mw / rho * (bd["Ntcw"] / (bd["Pt"] - bd["Do"]) + bd["B"] / bd["Dw"] ** 2) + 2 * mw**2 / (2 * rho)
    else:
        dPwi = (2 + 0.6 * bd["Ntcw"]) * mw**2 / (2 * rho)
    dPc = (Nb - 1) * dPbi * Rb * Rl                                         # iç çapraz akış bölmeleri
    dPw = Nb * dPwi * Rl                                                     # pencereler
    dPe = 2 * dPbi * (1 + bd["Ntcw"] / bd["Ntcc"]) * Rb * Rs                 # giriş + çıkış bölmeleri
    return dict(Re=Re, Pr=Pr, j=j, f=f, h_ideal=h_id, h=h, Jc=Jc, Jl=Jl, Jb=Jb, Js=Js, Jr=Jr,
                Rl=Rl, Rb=Rb, Rs=Rs, dPc=dPc, dPw=dPw, dPe=dPe, dP=dPc + dPw + dPe)


# =====================================================================
# 2. YOĞUŞMA EĞRİSİ (thermo ile veya tablodan)
# =====================================================================
FULLER_ATOM = {"C": 15.9, "H": 2.31, "O": 6.11, "N": 4.54, "S": 22.9, "Cl": 21.0}
FULLER_MOLECULE = {"H2": 6.12, "H2O": 13.1, "N2": 18.5, "O2": 16.3, "CO": 18.0, "CO2": 26.9}
LIQ_COLS = ["rho_L", "mu_L", "k_L", "cp_L"]


def fuller_D(T_C, P_bar, M_a, M_b, V_a, V_b):
    """Fuller ikili difüzyon katsayısı [m²/s]."""
    M_ab = 2 / (1 / M_a + 1 / M_b)
    return 1.43e-7 * (T_C + 273.15) ** 1.75 / (P_bar * M_ab**0.5 * (V_a ** (1 / 3) + V_b ** (1 / 3)) ** 2)


# Su–aseton–2-propanol için tek ve tutarlı NRTL seti (NIST ölçümleri + literatür; OSTI 799235, Tablo III):
# τ_ij = (g_ij − g_jj)/RT, α = 0,3. Burada b_ij = (g_ij − g_jj)/R [K]. CAS: su 7732-18-5, aseton 67-64-1, IPA 67-63-0.
# Doğrulama: IPA–su azeotropu 1 atm'de x_IPA = 0,677 (%87,5 kütle), 80,0 °C hesaplanır (literatür %87,4–87,7, 80,3–80,4 °C).
_R = 8.314462618
NRTL_DEFAULT = {("7732-18-5", "67-64-1"): (4648 / _R, 2095 / _R, 0.3),
                ("7732-18-5", "67-63-0"): (6768 / _R, 46.4 / _R, 0.3),
                ("67-64-1", "67-63-0"): (1756 / _R, -86.6 / _R, 0.3)}


def _nrtl_model(consts, zs, user_params):
    """Varsayılan NRTL seti (NRTL_DEFAULT) + ChemSep + kullanıcı parametreleri. Eksik ikili varsa durur."""
    from thermo import NRTL
    from thermo.interaction_parameters import IPDB
    names = list(consts.names)
    cas = list(consts.CASs)
    defaults = {}
    for (ca, cb), (b_ab, b_ba, al) in NRTL_DEFAULT.items():
        if ca in cas and cb in cas:
            defaults[(names[cas.index(ca)], names[cas.index(cb)])] = (b_ab, b_ba, al)
    user_params = {**defaults, **(user_params or {})}
    b = np.array(IPDB.get_ip_asymmetric_matrix("ChemSep NRTL", consts.CASs, "bij"), float)
    alpha = np.array(IPDB.get_ip_symmetric_matrix("ChemSep NRTL", consts.CASs, "alphaij"), float)
    for (a, c), (b_ac, b_ca, al) in (user_params or {}).items():
        i, j = names.index(a), names.index(c)
        b[i][j], b[j][i], alpha[i][j], alpha[j][i] = b_ac, b_ca, al, al
    missing = [f"{names[i]}–{names[j]}" for i in range(len(names)) for j in range(i + 1, len(names))
               if consts.Tcs[i] > 200 and consts.Tcs[j] > 200 and b[i][j] == 0 and b[j][i] == 0]
    if missing:
        raise ValueError("NRTL için eksik ikili parametre: " + ", ".join(missing) +
                         ". Doğrulanmış değer girin veya ideal sıvı modelini seçin.")
    return NRTL(T=300.0, xs=zs, tau_bs=b.tolist(), alpha_cs=alpha.tolist())


class Mixture:
    """Bileşenler + debiler için termodinamik paket (gaz: Peng–Robinson, sıvı: ideal/NRTL)."""

    def __init__(self, components, flows_kgh, liquid_model="ideal", nrtl_params=None):
        from thermo import ChemicalConstantsPackage, FlashVL, CEOSGas, GibbsExcessLiquid, PRMIX
        self.components = list(components)
        self.consts, self.props = ChemicalConstantsPackage.from_IDs(self.components)
        m = np.asarray(flows_kgh, float)
        if len(m) != len(self.components) or (m < 0).any() or m.sum() <= 0:
            raise ValueError("Bileşen debileri bileşen sayısıyla aynı uzunlukta, ≥ 0 ve toplamı > 0 olmalı.")
        self.flows = m
        self.MW = np.array(self.consts.MWs)
        n = m / self.MW / 3.6                                   # mol/s
        self.n_tot, self.zs = n.sum(), list(n / n.sum())
        self.noncond = np.array(self.consts.Tcs) < 200.0        # H2, N2 ... yoğuşmayan
        self.V_fuller = np.array([FULLER_MOLECULE.get(f, sum(FULLER_ATOM.get(el, 15.0) * k for el, k in at.items()))
                                  for f, at in zip(self.consts.formulas, self.consts.atomss)])
        c = self.consts
        gas = CEOSGas(PRMIX, eos_kwargs=dict(Tcs=c.Tcs, Pcs=c.Pcs, omegas=c.omegas),
                      HeatCapacityGases=self.props.HeatCapacityGases)
        kw = dict(VaporPressures=self.props.VaporPressures, HeatCapacityGases=self.props.HeatCapacityGases,
                  VolumeLiquids=self.props.VolumeLiquids, EnthalpyVaporizations=self.props.EnthalpyVaporizations,
                  equilibrium_basis="Psat", caloric_basis="Psat", T=300.0, P=1e5, zs=self.zs)
        if liquid_model == "nrtl":
            kw["GibbsExcessModel"] = _nrtl_model(c, self.zs, nrtl_params)
        self.flasher = FlashVL(c, self.props, liquid=GibbsExcessLiquid(**kw), gas=gas)

    def dew_point(self, P_bar):
        return self.flasher.flash(P=P_bar * 1e5, VF=1.0, zs=self.zs).T - 273.15

    def flash(self, T_C, P_bar):
        return self.flasher.flash(T=T_C + 273.15, P=P_bar * 1e5, zs=self.zs)

    def vapor_flows(self, T_C, P_bar):
        """T, P'de buhar fazındaki bileşen debileri [kg/h]."""
        s = self.flash(T_C, P_bar)
        if s.gas is None:
            return dict.fromkeys(self.components, 0.0)
        m_V = s.VF * self.n_tot * s.gas.MW() * 3.6
        return {c: m_V * w for c, w in zip(self.components, s.gas.ws())}


MW_FEED = {"acetone": 58.080, "isopropanol": 60.096, "hydrogen": 2.016, "water": 18.015}


def reactor_effluent(feed_tpy, hours=8000, w_ipa=0.870, conversion=0.90, basis="mixture"):
    """IPA dehidrojenasyonu kütle denkliği: (CH3)2CHOH -> (CH3)2CO + H2. Su tepkimeye girmez.
    feed_tpy: yıllık besleme [t/yıl]; basis 'mixture' = IPA–su karışımının toplamı, 'ipa' = yalnız IPA.
    w_ipa: beslemede IPA kütle kesri. Döndürür: reaktör çıkışı tablosu (kg/h, kmol/h, mol %, kütle %)."""
    kg = feed_tpy * 1000 / hours
    ipa = kg * w_ipa if basis == "mixture" else kg
    water = ipa * (1 - w_ipa) / w_ipa
    n_ipa = ipa / MW_FEED["isopropanol"]
    n = {"acetone": conversion * n_ipa, "isopropanol": (1 - conversion) * n_ipa,
         "hydrogen": conversion * n_ipa, "water": water / MW_FEED["water"]}
    df = pd.DataFrame({"Bileşen (thermo adı)": list(n), "Debi (kg/h)": [n[k] * MW_FEED[k] for k in n],
                       "kmol/h": list(n.values())})
    df["Mol %"] = 100 * df["kmol/h"] / df["kmol/h"].sum()
    df["Kütle %"] = 100 * df["Debi (kg/h)"] / df["Debi (kg/h)"].sum()
    return df


def condensation_curve(mix, T_in, T_out, P_in_bar, P_out_bar, n_zones=15):
    """Yoğuşma eğrisi: her noktada Q, buhar/sıvı debisi, iki fazın özellikleri ve
    Colburn–Hougen için y_NC (buhardaki yoğuşmayan mol kesri), M_c, D_AB."""
    T_dew = mix.dew_point(P_in_bar)
    if T_out < T_dew < T_in:
        # Bölgeler yoğuşma aralığında sıklaştırılır: kuru gaz soğuması Q–T'de neredeyse doğrusaldır (az bölge yeter),
        # ısının çoğu ve eğrinin bükülmesi çiy noktasının altındadır. Kuru bölgeye en az 2, en çok %20 bölge verilir.
        n_dry = max(2, round(0.2 * n_zones))
        T_list = list(np.linspace(T_in, T_dew, n_dry + 1)) + list(np.linspace(T_dew, T_out, n_zones - n_dry + 1))
    else:
        T_list = list(np.linspace(T_in, T_out, n_zones + 1))
    T_list = sorted(set(round(T, 3) for T in T_list), reverse=True)
    rows, H0 = [], None
    for T in T_list:
        P = P_in_bar + (P_out_bar - P_in_bar) * (T_in - T) / (T_in - T_out)
        s = mix.flash(T, P)
        H = s.H() * mix.n_tot
        H0 = H if H0 is None else H0
        gas, liq = s.gas, getattr(s, "liquid0", None)
        row = dict(T=T, P=P, Q_cum=(H0 - H) / 1000,
                   m_V=s.VF * mix.n_tot * gas.MW() * 3.6,
                   m_L=(1 - s.VF) * mix.n_tot * liq.MW() * 3.6 if liq else 0.0,
                   rho_V=gas.rho_mass(), mu_V=gas.mu(), k_V=gas.k(), cp_V=gas.Cp_mass() / 1000)
        if liq:
            row.update(rho_L=liq.rho_mass(), mu_L=liq.mu(), k_L=liq.k(), cp_L=liq.Cp_mass() / 1000)
        y = np.array(gas.zs)
        nc = mix.noncond
        y_nc = y[nc].sum()
        if nc.any() and 0 < y_nc < 1:
            yc, yg = y[~nc] / y[~nc].sum(), y[nc] / y_nc
            Mc, Mg = (yc * mix.MW[~nc]).sum(), (yg * mix.MW[nc]).sum()
            Vc, Vg = (yc * mix.V_fuller[~nc]).sum(), (yg * mix.V_fuller[nc]).sum()
            row.update(y_NC=y_nc, M_c=Mc, D_AB=fuller_D(T, P, Mc, Mg, Vc, Vg))
        else:
            row.update(y_NC=y_nc, M_c=np.nan, D_AB=np.nan)
        rows.append(row)
    curve = pd.DataFrame(rows)
    for c in LIQ_COLS:
        if c not in curve:
            curve[c] = np.nan
    curve[LIQ_COLS] = curve[LIQ_COLS].bfill()
    return curve, T_dew


# Tablo ile eğri girişi: kullanıcı dostu kolon adları -> iç adlar (μ cP olarak girilir)
TABLE_COLUMNS = {"T_C": "T", "P_bar": "P", "Q_kW": "Q_cum", "mV_kgh": "m_V", "mL_kgh": "m_L",
                 "rhoV": "rho_V", "muV_cP": "mu_V", "kV": "k_V", "cpV_kJkgK": "cp_V",
                 "rhoL": "rho_L", "muL_cP": "mu_L", "kL": "k_L", "cpL_kJkgK": "cp_L",
                 "yNC": "y_NC", "Mc": "M_c", "DAB_m2s": "D_AB"}
TABLE_REQUIRED = list(TABLE_COLUMNS)[:13]


def curve_from_table(df):
    """Simülatörden alınan eğri tablosu. yNC, Mc, DAB_m2s yoksa Colburn–Hougen yerine SBG kullanılır."""
    missing = [c for c in TABLE_REQUIRED if c not in df.columns]
    if missing:
        raise ValueError(f"Eğri tablosunda eksik kolon: {missing}")
    c = df.rename(columns=TABLE_COLUMNS)[[v for k, v in TABLE_COLUMNS.items() if k in df.columns]].astype(float)
    c = c.sort_values("T", ascending=False).reset_index(drop=True)
    c["mu_V"] /= 1000
    c["mu_L"] /= 1000
    for col in ("y_NC", "M_c", "D_AB"):
        if col not in c:
            c[col] = np.nan
    c[LIQ_COLS] = c[LIQ_COLS].bfill()
    T_dew = float(c.loc[c.m_L <= 1e-6, "T"].min()) if (c.m_L <= 1e-6).any() else float(c["T"].iloc[0])
    return c, T_dew


def fill_diffusion(curve, mix):
    """Tablo eğrisinde y_NC var ama M_c/D_AB yoksa: besleme bileşiminden M_c ve Fuller D_AB hesaplanır."""
    c = curve.copy()
    if c["y_NC"].notna().any() and c["D_AB"].isna().all() and mix.noncond.any():
        z = np.array(mix.zs)
        nc = mix.noncond
        zc, zg = z[~nc] / z[~nc].sum(), z[nc] / z[nc].sum()
        Mc, Mg = (zc * mix.MW[~nc]).sum(), (zg * mix.MW[nc]).sum()
        Vc, Vg = (zc * mix.V_fuller[~nc]).sum(), (zg * mix.V_fuller[nc]).sum()
        c["M_c"] = Mc
        c["D_AB"] = [fuller_D(T, P, Mc, Mg, Vc, Vg) for T, P in zip(c["T"], c["P"])]
    return c


def curve_to_table(curve):
    """İç eğri -> kullanıcı tablosu (şablon olarak indirilebilir)."""
    t = curve.copy()
    t["mu_V"] *= 1000
    t["mu_L"] *= 1000
    return t.rename(columns={v: k for k, v in TABLE_COLUMNS.items()})[list(TABLE_COLUMNS)]


# =====================================================================
# 3. KABUK TARAFI ISI TRANSFER KATSAYILARI
# =====================================================================
def h_gas_kern(m_V, mu_V, cp_V, k_V, sg):
    """Gaz (duyulur), Kern: Nu = 0,36 Re^0,55 Pr^(1/3). m_V kg/s."""
    Re = sg["De"] * (m_V / sg["As"]) / mu_V
    Pr = cp_V * 1000 * mu_V / k_V
    return dict(h=0.36 * Re**0.55 * Pr ** (1 / 3) * k_V / sg["De"], Re=Re, Pr=Pr)


def h_film_nusselt(m_L, rho_L, rho_V, mu_L, k_L, sg, N):
    """Yatay demette kondensat filmi, Nusselt–Kern. Γ = ṁ_L/(L·N_t)."""
    gamma = m_L / (sg["L"] * N)
    h = 0.95 * k_L * (rho_L * (rho_L - rho_V) * G / (mu_L * gamma)) ** (1 / 3) * sg["Nr"] ** (-1 / 6)
    return dict(h=h, Gamma=gamma)


def h_colburn_hougen(z, hG, Pr, hF, R_other, t_water, y_cond_at):
    """Colburn–Hougen: arayüz sıcaklığı T_i enerji dengesinden çözülür.
        hG·(T_G − T_i) + N·M_c·λ = U_w·(T_i − t_su),   N = K·ln[(1 − y_i)/(1 − y_G)]
    T_i eğrinin en düşük sıcaklığının altına inerse y_i son değerde sabit kalır (muhafazakâr).
    Çözüm yoksa None döner (çağıran SBG'ye düşer)."""
    T_G = z["T"]
    if not np.isfinite(z.get("D_AB", np.nan)) or z["lam"] <= 0 or T_G <= t_water:
        return None
    Sc = z["mu_V"] / (z["rho_V"] * z["D_AB"])
    c = z["P"] * 1e5 / (R_GAS * (T_G + 273.15))
    K = hG / (z["rho_V"] * z["cp_V"] * 1000) * (Pr / Sc) ** (2 / 3) * c
    U_w = 1 / (1 / hF + R_other)
    y_G = y_cond_at(T_G)

    def balance(T_i):
        N = K * math.log((1 - y_cond_at(T_i)) / (1 - y_G))
        return hG * (T_G - T_i) + N * z["M_c"] / 1000 * z["lam"] - U_w * (T_i - t_water)

    try:
        T_i = brentq(balance, t_water, T_G)
    except ValueError:
        return None
    q = U_w * (T_i - t_water)
    inv_h = (T_G - t_water) / q - R_other
    if inv_h <= 0:
        return None
    return dict(h=1 / inv_h, T_i=T_i, Sc=Sc, K=K, y_G=y_G, y_i=y_cond_at(T_i), q=q, U_w=U_w)


# =====================================================================
# 4. BÖLGESEL RATING (alan, kabuk ΔP, ε-NTU, duvar sıcaklıkları)
# =====================================================================


def rate_condensing(tube, geo, curve, h_method="colburn_hougen", h_fixed=450.0, h_mult=1.0,
                    Rf_shell=0.00023, k_wall=60.0, K_pass=2.5, shells=1, shell_method="bell_delaware",
                    n_ss=None, baffle_cut=0.25, two_phase_C=TWO_PHASE_C):
    """Yoğuşan kabuk için bölgesel rating. Sonuç sözlüğü + bölge tablosu döner.
    shell_method: 'bell_delaware' (kaçak/bypass düzeltmeli, varsayılan) veya 'kern' (karşılaştırma).
    Kabuk ΔP: her bölge kabuk boyunun A_i/ΣA_i kadarını kaplar; bölge ΔP'si = φ_G²(yerel x)·ΔP_gaz,
    giriş ve çıkış bölmeleri ilk ve son bölgenin koşullarıyla; nozul kayıpları eklenir."""
    Q_kW = float(curve["Q_cum"].iloc[-1])
    w = water_side_balance(tube, Q_kW)
    T_h_in, T_h_out = float(curve["T"].iloc[0]), float(curve["T"].iloc[-1])
    lf = lmtd_F(T_h_in, T_h_out, w["T_in"], w["T_out"], shells)
    ts = tube_side(w, geo, shells, K_pass)
    sg = shell_geometry(geo)
    use_bd = shell_method == "bell_delaware"
    bd = bd_geometry(dict(geo, Db=geo.get("Db") or sg["Db"] * 1000), Bc=baffle_cut, n_ss=n_ss) if use_bd else None
    dr = sg["Do"] / ts["Di"]
    R_wall = sg["Do"] * math.log(dr) / (2 * k_wall)
    R_other = Rf_shell + R_wall + dr * w["Rf"] + dr / ts["h"]          # kabuk filmi dışındaki dirençler
    R_tube_side = dr * w["Rf"] + dr / ts["h"]
    C_c = w["m"] / 3600 * w["cp"]                                     # kW/K
    t_water = lambda Q: w["T_out"] - Q / C_c                          # ters akış su profili
    T_up, y_up = curve["T"].values[::-1], (1 - curve["y_NC"]).values[::-1]
    y_cond_at = lambda T: float(np.interp(T, T_up, y_up))
    ch_possible = curve["D_AB"].notna().any()

    rows = []
    for i in range(len(curve) - 1):
        a, b = curve.iloc[i], curve.iloc[i + 1]
        dQ, dT = b.Q_cum - a.Q_cum, a["T"] - b["T"]
        z = {k: 0.5 * (a[k] + b[k]) for k in curve.columns}
        m_V, m_L = z["m_V"] / 3600, z["m_L"] / 3600
        gas = h_gas_kern(m_V, z["mu_V"], z["cp_V"], z["k_V"], sg)
        h_kern = gas["h"]
        if use_bd:                                                    # tek fazlı gaz referansı: Bell–Delaware
            bdg = bell_delaware(m_V, z["rho_V"], z["mu_V"], z["cp_V"], z["k_V"], bd)
            gas = dict(h=bdg["h"], Re=bdg["Re"], Pr=bdg["Pr"])
        hG = gas["h"]
        rec = dict(T_a=a["T"], T_b=b["T"], dQ_kW=dQ, t_su_a=t_water(a.Q_cum), t_su_b=t_water(b.Q_cum),
                   m_V=z["m_V"], m_L=z["m_L"], Re_v=gas["Re"], Pr_v=gas["Pr"], h_G=hG,
                   Gamma=np.nan, h_F=np.nan, Z=1.0, h_SBG=hG, T_i=np.nan, Sc=np.nan, K=np.nan,
                   y_G=np.nan, y_i=np.nan, lam=np.nan, model="kuru duvar (gaz)")
        if b.m_L > 1e-6:                                              # bölgede yoğuşma var
            film = h_film_nusselt(m_L, z["rho_L"], z["rho_V"], z["mu_L"], z["k_L"], sg, geo["N"])
            Z = min(m_V * z["cp_V"] * dT / dQ, 1.0)
            h_sbg = 1 / (1 / film["h"] + Z / hG)
            rec.update(Gamma=film["Gamma"], h_F=film["h"], Z=Z, h_SBG=h_sbg, model="SBG")
            if h_method == "colburn_hougen" and ch_possible:
                z["lam"] = (dQ - (m_V * z["cp_V"] + m_L * z["cp_L"]) * dT) * 1000 / ((b.m_L - a.m_L) / 3600)
                rec["lam"] = z["lam"]
                ch = h_colburn_hougen(z, hG, gas["Pr"], film["h"], R_other, t_water(z["Q_cum"]), y_cond_at)
                if ch:
                    rec.update(T_i=ch["T_i"], Sc=ch["Sc"], K=ch["K"], y_G=ch["y_G"], y_i=ch["y_i"],
                               h_CH=ch["h"], model="Colburn–Hougen")
        h = {"fixed": h_fixed}.get(h_method, rec.get("h_CH", rec["h_SBG"]) * h_mult)
        Uc, Ud = overall_U(h, ts["h"], Rf_shell, w["Rf"], sg["Do"], ts["Di"], k_wall)

        d1, d2 = a["T"] - t_water(a.Q_cum), b["T"] - t_water(b.Q_cum)
        if min(d1, d2) <= 0:
            raise ValueError(f"{a['T']:.1f}→{b['T']:.1f} °C bölgesinde sıcaklık kesişimi.")
        lm = d1 if abs(d1 - d2) < 1e-9 else (d1 - d2) / math.log(d1 / d2)
        A_i = dQ * 1000 / (Ud * lf["F"] * lm)

        # ε-NTU (ters akış) — bağımsız kontrol
        C_h = dQ / max(dT, 1e-9)
        C_min, C_max = min(C_h, C_c), max(C_h, C_c)
        eps = dQ / (C_min * (a["T"] - t_water(b.Q_cum)))
        NTU = ht.NTU_from_effectiveness(eps, C_min / C_max, subtype="counterflow")

        # kabuk ΔP (tüm kabuk için, bölgenin yerel gaz/sıvı debileriyle): gaz tek başına × iki faz çarpanı
        if use_bd:
            g_ = bell_delaware(m_V, z["rho_V"], z["mu_V"], z["cp_V"], z["k_V"], bd)
            dp_G, dp_G_core, dp_G_end = g_["dP"], g_["dPc"] + g_["dPw"], g_["dPe"]
            dp_L = (bell_delaware(m_L, z["rho_L"], z["mu_L"], z["cp_L"], z["k_L"], bd)["dP"] if m_L > 1e-6 else 0.0)
        else:
            dp_G = kern_dp(m_V, z["rho_V"], z["mu_V"], sg)
            nb = sg["n_baffle"]
            dp_G_core, dp_G_end = dp_G * (nb - 1) / (nb + 1), dp_G * 2 / (nb + 1)
            dp_L = kern_dp(m_L, z["rho_L"], z["mu_L"], sg) if m_L > 1e-6 else 0.0
        X = math.sqrt(dp_L / dp_G) if m_L > 1e-6 else 0.0
        phi2 = 1 + two_phase_C * X + X**2

        # tüp metal sıcaklığı (ısıl genleşme kontrolü için)
        q = Ud * lm * lf["F"]                                         # W/m² (dış yüzey)
        t_avg = 0.5 * (t_water(a.Q_cum) + t_water(b.Q_cum))
        rec.update(h_shell=h, Uc=Uc, Ud=Ud, LMTD=lm, A_m2=A_i, eps=eps, NTU=NTU, Cr=C_min / C_max,
                   A_eNTU_m2=NTU * C_min * 1000 / Ud, X=X, phi2=phi2, dP_G_kPa=dp_G / 1000,
                   dP_2f_kPa=phi2 * dp_G / 1000, dP_core_2f_kPa=phi2 * dp_G_core / 1000,
                   dP_end_2f_kPa=phi2 * dp_G_end / 1000, h_G_kern=h_kern,
                   T_tube_metal=t_avg + q * (R_tube_side + R_wall / 2), T_hot_avg=0.5 * (a["T"] + b["T"]))
        rows.append(rec)

    zones = pd.DataFrame(rows)
    A_req = zones.A_m2.sum()
    A = sg["A"] * shells
    wgt = zones.A_m2 / A_req
    # kabuk ΔP, fiziksel toplama: bölge i kabuk boyunun wgt_i kadarını kaplar (iç bölmeler + pencereler);
    # giriş bölmesi ilk bölgenin, çıkış bölmesi son bölgenin koşullarında; nozullar: giriş 1,0, çıkış 0,5 hız yükü
    zones["dP_seg_kPa"] = zones.dP_core_2f_kPa * wgt
    dp_core = zones.dP_seg_kPa.sum()
    dp_ends = 0.5 * zones.dP_end_2f_kPa.iloc[0] + 0.5 * zones.dP_end_2f_kPa.iloc[-1]
    inl, out = curve.iloc[0], curve.iloc[-1]
    m_tot = inl.m_V + inl.m_L
    x_out = out.m_V / (out.m_V + out.m_L)
    rho_out = 1 / (x_out / out.rho_V + (1 - x_out) / out.rho_L) if out.m_L > 0 else out.rho_V
    n_in = pick_nozzle(m_tot, inl.rho_V, RHOV2_SHELL_INLET)
    n_out = pick_nozzle(m_tot, rho_out, RHOV2_SHELL_INLET)
    dp_noz = (1.0 * n_in["rhov2"] / 2 + 0.5 * n_out["rhov2"] / 2) / 1000
    # momentum (yoğuşmayla yavaşlama) basınç geri kazanımı, homojen model — bilgi amaçlı, toplamdan düşülmez
    Gm = m_tot / 3600 / (bd["Sm"] if use_bd else sg["As"])
    dp_mom = Gm**2 * (1 / rho_out - 1 / inl.rho_V) / 1000
    dp_shell = shells * (dp_core + dp_ends) + dp_noz
    return dict(
        Q_kW=Q_kW, water=w, T_hot=(T_h_in, T_h_out), lf=lf, ts=ts, sg=sg, geo=geo, shells=shells,
        h_method=h_method, R_wall=R_wall, R_other=R_other, Rf_shell=Rf_shell, k_wall=k_wall,
        A=A, A_req=A_req, over_design=(A / A_req - 1) * 100,
        A_LMTD_cf=(zones.A_m2 * lf["F"]).sum(), A_eNTU_cf=zones.A_eNTU_m2.sum(),
        h_shell=(zones.h_shell * wgt).sum(), h_shell_min=zones.h_shell.min(),
        Ud=Q_kW * 1000 / (lf["F"] * lf["LMTD"] * A_req),
        Uc=Q_kW * 1000 / (lf["F"] * lf["LMTD"] * (zones.dQ_kW * 1000 / (zones.Uc * lf["F"] * zones.LMTD)).sum()),
        dp_tube_kPa=ts["dp_kPa"], dp_shell_kPa=dp_shell,
        dp_shell_gas_kPa=shells * (zones.dP_G_kPa * wgt).sum(),
        dp_shell_parts=dict(bolmeler_pencereler=shells * dp_core, giris_cikis_bolmeleri=shells * dp_ends,
                            nozullar=dp_noz, momentum_geri_kazanim=dp_mom),
        shell_method=shell_method, bd=bd, two_phase_C=two_phase_C,
        T_shell_metal=(zones.T_hot_avg * wgt).sum(), T_tube_metal=(zones.T_tube_metal * wgt).sum(),
        zones=zones)


def rate_single_lmtd(r, h_shell):
    """Karşılaştırma: aynı Q ve su tarafıyla tek LMTD + sabit kabuk h (orijinal kod yaklaşımı)."""
    Uc, Ud = overall_U(h_shell, r["ts"]["h"], r["Rf_shell"], r["water"]["Rf"], r["sg"]["Do"], r["ts"]["Di"], r["k_wall"])
    A_req = r["Q_kW"] * 1000 / (Ud * r["lf"]["F"] * r["lf"]["LMTD"])
    (Th1, Th2), Tc1 = r["T_hot"], r["water"]["T_in"]
    C_h, C_c = r["Q_kW"] / (Th1 - Th2), r["water"]["m"] / 3600 * r["water"]["cp"]
    C_min, C_max = min(C_h, C_c), max(C_h, C_c)
    eps = r["Q_kW"] / (C_min * (Th1 - Tc1))
    A_eNTU = ht.NTU_from_effectiveness(eps, C_min / C_max, subtype="S&T", n_shell_tube=r["shells"]) * C_min * 1000 / Ud
    return dict(A_req=A_req, A_eNTU=A_eNTU, Ud=Ud, over_design=(r["A"] / A_req - 1) * 100)


# =====================================================================
# 4b. TÜP YERLEŞİMİ (geometrik sayım) VE BOYUTLANDIRMA
# =====================================================================
# Tipik demet–kabuk çap boşlukları (mm). Towler & Sinnott Şekil 19.12'den okunur; kendi değerinizi girin.
BUNDLE_CLEARANCE = {"Sabit ayna (BEM)": 15.0, "U-tüp (BEU)": 15.0,
                    "Yüzer kafa, split-ring (AES)": 65.0, "Yüzer kafa, pull-through (AET)": 95.0}


def tube_layout(D_otl, Do, Pt, layout="triangular", passes=1):
    """Ayna üzerine tüp yerleşimi. D_otl: dış tüp sınırı çapı (= demet çapı Db), mm.
    Geçiş bölme levhası şeritlerinde tüp konmaz: tüp merkezi şerit ekseninden en az Do/2 + 9,5 mm uzakta
    (Ds < 559 mm ise 7,9 mm). En çok tüp sığan ızgara kaydırması seçilir.
    Döndürür: (tüp sayısı, merkez koordinatları [N×2], şeritler)."""
    R = (D_otl - Do) / 2
    gap = Do / 2 + (7.9 if D_otl < 559 else 9.5)
    if layout == "triangular":
        dx, dy, shift = Pt, Pt * math.sqrt(3) / 2, Pt / 2
    elif layout == "rotated square":
        dx, dy, shift = Pt * math.sqrt(2), Pt * math.sqrt(2) / 2, Pt * math.sqrt(2) / 2
    else:
        dx, dy, shift = Pt, Pt, 0.0
    lanes_h = {1: [], 2: [0.0], 4: [0.0], 6: [-R / 3, R / 3], 8: [-R / 2, 0.0, R / 2]}[passes]
    lanes_v = [0.0] if passes in (4, 6, 8) else []
    best = None
    for ox in (0.0, dx / 2):
        for oy in (0.0, dy / 2):
            rows = np.arange(-math.ceil(R / dy) - 1, math.ceil(R / dy) + 2)
            pts = []
            for j in rows:
                y = j * dy + oy
                x0 = ox + (shift if j % 2 else 0.0)
                xs = np.arange(-math.ceil(R / dx) - 1, math.ceil(R / dx) + 2) * dx + x0
                pts.append(np.column_stack([xs, np.full_like(xs, y)]))
            P = np.vstack(pts)
            keep = np.hypot(P[:, 0], P[:, 1]) <= R
            for yl in lanes_h:
                keep &= np.abs(P[:, 1] - yl) >= gap
            for xl in lanes_v:
                keep &= np.abs(P[:, 0] - xl) >= gap
            if best is None or keep.sum() > len(best):
                best = P[keep]
    return len(best), best, dict(h=lanes_h, v=lanes_v)


def otl_for_tubes(N_target, Do, Pt, layout, passes):
    """N_target tüpün sığdığı en küçük dış tüp sınırı çapı (mm, 5 mm adım)."""
    lo, hi = 2 * Do, 4000.0
    while hi - lo > 5:
        mid = 0.5 * (lo + hi)
        if tube_layout(mid, Do, Pt, layout, passes)[0] >= N_target:
            hi = mid
        else:
            lo = mid
    return math.ceil(hi / 5) * 5


# =====================================================================
# 5. MEKANİK: TEMA TİPİ, NOZULLAR, DESTEK AÇIKLIĞI, ET KALINLIĞI
# =====================================================================
# Standart boru nominal çapları (inç) ve STD et iç çapları (mm)
NOZZLE_SIZES = [(2, 52.5), (3, 77.9), (4, 102.3), (6, 154.1), (8, 202.7), (10, 254.5), (12, 304.8),
                (14, 336.6), (16, 387.4), (18, 438.2), (20, 489.0), (24, 590.6)]
RHOV2_SHELL_INLET = 2232.0     # kg/m·s², TEMA RCB-4.61: üstünde çarpma plakası
RHOV2_TUBE_NOZZLE = 8928.0     # kg/m·s², TEMA: tüp tarafı nozul sınırı
# TEMA RCB-4.52: desteksiz en büyük düz tüp açıklığı (karbon çeliği), Do mm -> mm
TEMA_MAX_SPAN = {15.88: 1321, 19.05: 1524, 25.4: 1880, 31.75: 2235, 38.1: 2540}


def pick_nozzle(m_kgh, rho, rhov2_max):
    """ρv² sınırını sağlayan en küçük standart nozul."""
    for nps, d_mm in NOZZLE_SIZES:
        v = m_kgh / 3600 / (rho * math.pi / 4 * (d_mm / 1000) ** 2)
        if rho * v**2 <= rhov2_max:
            return dict(NPS=nps, D_mm=d_mm, v=v, rhov2=rho * v**2)
    nps, d_mm = NOZZLE_SIZES[-1]
    v = m_kgh / 3600 / (rho * math.pi / 4 * (d_mm / 1000) ** 2)
    return dict(NPS=nps, D_mm=d_mm, v=v, rhov2=rho * v**2)


PLATES_MM = [5, 6, 8, 10, 12, 14, 16, 18, 20, 25, 30]
TEMA_MIN_SHELL_CB = [(305, 6.4), (737, 7.9), (991, 9.5), (1524, 11.1), (2540, 12.7)]   # TEMA CB-3.13, karbon çeliği levha


def shell_plate_thickness(Ds_mm, P_barg, S_MPa, E_joint, CA_mm):
    """ASME UG-27 + korozyon payı ile TEMA CB-3.13 asgari kalınlığının büyüğü, standart levhaya yuvarlanmış [mm]."""
    P = P_barg / 10
    t_req = P * Ds_mm / 2 / (S_MPa * E_joint - 0.6 * P) + CA_mm
    t_min = next(t for d, t in TEMA_MIN_SHELL_CB if Ds_mm <= d)
    return next((p for p in PLATES_MM if p >= max(t_req, t_min)), PLATES_MM[-1])


def mechanical(r, curve, T_dew, mech):
    """TEMA tipi önerisi, nozullar, tüp destek açıklığı, et kalınlıkları.
    mech: P_design_shell_barg, P_design_tube_barg, S_shell_MPa, S_tube_MPa, E_joint, CA_shell_mm, alpha"""
    g, sg, w = r["geo"], r["sg"], r["water"]
    out = {}
    # --- ısıl genleşme ve TEMA tipi (TEMA T-4 ortalama metal sıcaklıkları; A.15 mantığı, J = 1) ---
    # Sabit aynalı yapıda kabuk ve tüpler aynalara bağlıdır; serbest uzama farkı ΔL eksenel kuvvete dönüşür.
    # Basitleştirme: aynalar rijit kabul edilir (üst sınır). Kuvvet: F = ΔL·E / [L·(1/A_kabuk + 1/A_tüp)].
    dT_metal = r["T_shell_metal"] - r["T_tube_metal"]
    delta_mm = mech["alpha"] * sg["L"] * abs(dT_metal) * 1000
    E_mod = mech.get("E_MPa", 200e3) * 1e6
    t_sh = shell_plate_thickness(g["Ds"], mech["P_design_shell_barg"], mech["S_shell_MPa"], mech["E_joint"], mech["CA_shell_mm"]) / 1000
    A_sh = math.pi * (sg["Ds"] + t_sh) * t_sh
    Di_ = sg["Do"] - 2 * g["t"] / 1000
    a_t1 = math.pi / 4 * (sg["Do"]**2 - Di_**2)
    A_tb = g["N"] * a_t1
    F_ax = delta_mm / 1000 * E_mod / (sg["L"] * (1 / A_sh + 1 / A_tb))
    sig_sh, sig_tb = F_ax / A_sh / 1e6, F_ax / A_tb / 1e6                  # MPa (sıcak taraf basmada, soğuk taraf çekmede)
    joint_kN = F_ax / A_tb * a_t1 / 1000
    ok_exp = sig_sh <= mech["S_shell_MPa"] and sig_tb <= mech["S_tube_MPa"]
    hot = "kabuk" if dT_metal > 0 else "tüpler"
    if ok_exp:
        tema = "BEM (sabit aynalı)"
        why = (f"Ortalama metal sıcaklıkları (TEMA T-4) kabuk {r['T_shell_metal']:.0f} °C, tüp {r['T_tube_metal']:.0f} °C; "
               f"serbest uzama farkı {delta_mm:.1f} mm. Aynalar rijit kabul edilse bile kabuk {sig_sh:.0f} MPa, tüp {sig_tb:.0f} MPa "
               f"eksenel gerilme (izin verilen {mech['S_shell_MPa']:.0f} / {mech['S_tube_MPa']:.0f} MPa): genleşme elemanı gerekmez. "
               "Kabuk akışkanı temiz, mekanik temizlik gerektirmez; bu nedenle sökülebilir demet gereksiz.")
    else:
        tema = "AES (yüzer kafa) veya BEM + genleşme körüğü"
        why = (f"Serbest uzama farkı {delta_mm:.1f} mm; kısıtlı durumda kabuk {sig_sh:.0f} MPa, tüp {sig_tb:.0f} MPa — izin verilen "
               "gerilme aşılıyor: sabit ayna yalnız esnek kabuk elemanıyla (TEMA RCB-8) kullanılabilir veya yüzer kafa seçilir.")
    out.update(T_shell_metal=r["T_shell_metal"], T_tube_metal=r["T_tube_metal"], dT_metal=dT_metal,
               delta_mm=delta_mm, F_axial_kN=F_ax / 1000, sigma_shell_MPa=sig_sh, sigma_tube_MPa=sig_tb,
               joint_load_kN=joint_kN, hot_side=hot, expansion_ok=ok_exp, t_shell_used_mm=t_sh * 1000,
               tema=tema, tema_reason=why,
               tema_note="Soğutma suyu tüpte: tüp içi bonnet sökülerek temizlenir (TEMA E-4.12).")
    # --- nozullar ---
    inlet, outlet = curve.iloc[0], curve.iloc[-1]
    m_tot = inlet.m_V + inlet.m_L
    x_out = outlet.m_V / m_tot
    rho_mix_out = 1 / (x_out / outlet.rho_V + (1 - x_out) / outlet.rho_L)          # homojen karışım
    saturated_inlet = r["T_hot"][0] <= T_dew + 1.0
    n_in = pick_nozzle(m_tot, inlet.rho_V, RHOV2_SHELL_INLET)
    n_out = pick_nozzle(m_tot, rho_mix_out, RHOV2_SHELL_INLET)
    n_w = pick_nozzle(w["m"], w["rho"], RHOV2_TUBE_NOZZLE)
    out["nozzles"] = pd.DataFrame([
        ("Kabuk girişi (buhar)", n_in["NPS"], n_in["D_mm"], n_in["v"], n_in["rhov2"], RHOV2_SHELL_INLET),
        ("Kabuk çıkışı (iki faz, homojen)", n_out["NPS"], n_out["D_mm"], n_out["v"], n_out["rhov2"], RHOV2_SHELL_INLET),
        ("Su giriş / çıkış", n_w["NPS"], n_w["D_mm"], n_w["v"], n_w["rhov2"], RHOV2_TUBE_NOZZLE),
    ], columns=["Nozul", "NPS (inç)", "İç çap (mm)", "Hız (m/s)", "ρv² (kg/m·s²)", "Sınır ρv²"])
    out["impingement"] = ("GEREKLİ — giriş doymuş/iki fazlı (TEMA: doymuş buhar ve iki fazda her zaman)"
                          if saturated_inlet else
                          "Gerekmez — kızgın tek faz buhar ve seçilen nozulda ρv² ≤ 2232 kg/m·s²")
    # --- destek açıklığı (segmental baffle: pencere tüpleri her iki baffle'da bir desteklenir) ---
    span = 2 * g["B"]
    span_max = TEMA_MAX_SPAN.get(g["Do"], np.nan)
    out.update(span_mm=span, span_max_mm=span_max, span_ok=bool(span <= span_max))
    # --- et kalınlıkları (ASME VIII Div.1 UG-27, çevresel gerilme): t = P·R / (S·E − 0,6·P) ---
    P_s = mech["P_design_shell_barg"] / 10          # MPa
    R_s = g["Ds"] / 2
    t_shell = P_s * R_s / (mech["S_shell_MPa"] * mech["E_joint"] - 0.6 * P_s)
    P_t = mech["P_design_tube_barg"] / 10
    R_t = (g["Do"] - 2 * g["t"]) / 2
    t_tube = P_t * R_t / (mech["S_tube_MPa"] * 1.0 - 0.6 * P_t)                   # dikişsiz tüp: E = 1
    out.update(t_shell_calc=t_shell, t_shell_req=t_shell + mech["CA_shell_mm"], t_tube_calc=t_tube,
               t_tube_ok=bool(t_tube <= g["t"]))
    # --- akış kaynaklı titreşim (TEMA Bölüm 6) ve ayna kalınlığı (TEMA Ek A.13, RCB-7.11) ---
    out["vib"] = vibration_check(r, curve, dict(mech, sigma_tube_MPa=sig_tb))
    out["tubesheet"] = tubesheet_thickness(r, mech)
    return out


# =====================================================================
# 5b. AKIŞ KAYNAKLI TİTREŞİM (TEMA 9. bs., Bölüm 6) ve AYNA KALINLIĞI (RCB-7, Ek A.13)
# =====================================================================
# TEMA Bölüm 6 bağıntıları İngiliz birimleriyle verilir (in, ft/s, lb/ft, lb/ft³, psi); hesap bu birimlerde yapılır.
_IN, _FT, _LBFT, _LBFT3, _PSI = 1 / 25.4, 3.28084, 0.671969, 0.062428, 145.0377e-6   # mm→in, m→ft, kg/m→lb/ft, kg/m³→lb/ft³, Pa→psi
TEMA_C_SPAN = {"ss-ss": 9.9, "ankastre-ss": 15.42, "ankastre-ankastre": 22.37}       # Tablo V-5.3
TEMA_K_AXIAL = {"ss-ss": math.pi, "ankastre-ss": 4.49, "ankastre-ankastre": 2 * math.pi}  # V-6.1
TEMA_CM_TRI = [(1.20, 2.25), (1.25, 2.03), (1.33, 1.78), (1.50, 1.47)]             # Tablo V-8, üçgen
TEMA_CM_SQ = [(1.20, 1.87), (1.25, 1.72), (1.33, 1.56), (1.50, 1.35)]              # Tablo V-8, kare
TEMA_CL_30 = [(1.20, 0.090), (1.25, 0.091), (1.33, 0.065), (1.50, 0.025)]          # Tablo V-11.2, 30°
VIB_MARGIN = 0.8        # V/Vc için önerilen tasarım payı (TEMA yalnız V < Vc ister; %20 pay mühendislik tercihi)


def _interp(tab, x):
    xs, ys = zip(*tab)
    return float(np.interp(x, xs, ys))


def baffle_thickness_mm(Ds_mm, span_mm):
    """TEMA R-4.41 en küçük baffle kalınlığı (merkez baffle aralığına göre)."""
    cols = [356, 711, 965, 1524, 2540]
    rows = [(610, [3.2, 4.8, 6.4, 6.4, 9.5]), (914, [4.8, 6.4, 7.9, 9.5, 12.7]), (1219, [6.4, 9.5, 9.5, 12.7, 15.9]),
            (1524, [9.5, 9.5, 12.7, 15.9, 19.1]), (1e9, [9.5, 12.7, 15.9, 15.9, 19.1])]
    j = next((i for i, c in enumerate(cols) if Ds_mm <= c), len(cols) - 1)
    return next(v[j] for lim, v in rows if span_mm <= lim)


def tema_D_factor(X, angle, pitch_ratio):
    """Tablo V-10.1 kritik hız katsayısı D. (30° satırı: 8,86·(P/do − 0,9)·X^0,34 / X^0,5.)"""
    if angle == 30:
        return 8.86 * (pitch_ratio - 0.9) * (X**0.34 if X <= 1 else X**0.5)
    if angle == 60:
        return 2.80 * (X**0.17 if X <= 1 else X**0.5)
    if angle == 90:
        return 2.10 * X**0.15 if X <= 0.7 else 2.35 * X**0.5
    return 4.13 * (pitch_ratio - 0.5) * X**0.5          # 45°


def _local_shell_state(r, curve, x_frac):
    """Kabuk boyunca x_frac (0 = giriş, 1 = çıkış) konumundaki gaz/sıvı debisi ve yoğunlukları.
    Konum, bölge alanlarının birikimli payından Q'ya, Q'dan yoğuşma eğrisine eşlenir."""
    z = r["zones"]
    A_cum = np.concatenate([[0], np.cumsum(z.A_m2.values)])
    Q_cum = np.concatenate([[0], np.cumsum(z.dQ_kW.values)])
    Q = float(np.interp(x_frac * A_cum[-1], A_cum, Q_cum))
    c = curve.sort_values("Q_cum")
    g = {k: float(np.interp(Q, c.Q_cum, c[k])) for k in ("m_V", "m_L", "rho_V", "rho_L", "mu_V", "mu_L", "T", "P")}
    W = g["m_V"] + g["m_L"]
    vV, vL = g["m_V"] / g["rho_V"], g["m_L"] / g["rho_L"]
    eps = vV / (vV + vL)
    g.update(W=W, eps_g=eps, rho_TP=g["rho_L"] * (1 - eps) + g["rho_V"] * eps)
    return g


def vibration_check(r, curve, mech=None, strouhal=None, gamma=1.15):
    """TEMA Bölüm 6 akış kaynaklı titreşim taraması, her tipik desteksiz açıklık için:
    doğal frekans (V-5.3, V-6, V-7), sönüm (V-8, buhar ifadesi), akışkan-elastik parametre (V-4.2),
    kritik hız (V-10), girdap kopması ve türbülans çırpınması genliği (V-11), akustik rezonans (V-12).
    Referans çapraz akış hızı: tüm kabuk debisi Bell–Delaware çapraz akış alanından geçer kabulü
    (kaçak/bypass düşülmez → V-9.21'den büyük, güvenli tarafta)."""
    g, sg, bd = r["geo"], r["sg"], r["bd"]
    mech = mech or {}
    E = mech.get("E_MPa", 200e3) * 1e6 * _PSI                    # psi
    do_m, di_m = g["Do"] / 1000, (g["Do"] - 2 * g["t"]) / 1000
    do, di = do_m * 1000 * _IN, di_m * 1000 * _IN
    I = math.pi / 64 * (do**4 - di**4)                           # in⁴
    At = math.pi / 4 * (do**2 - di**2)                           # in²
    wt = 7850 * math.pi / 4 * (do_m**2 - di_m**2) * _LBFT        # lb/ft
    wfi = 0.00545 * r["water"]["rho"] * _LBFT3 * di**2
    pr = g["Pt"] / g["Do"]
    ang = bd["ang"]
    Cm = _interp(TEMA_CM_TRI if ang in (30, 60) else TEMA_CM_SQ, pr)
    B, B_end = g["B"], bd["B_in"] * 1000
    Nb = bd["Nb"]
    tb = baffle_thickness_mm(g["Ds"], B)
    Ltot = g["L"]
    # Tipik açıklıklar: (ad, uzunluk mm, uç koşulu, konum x/L)
    spans = [("Giriş açıklığı (ayna–1. baffle)", B_end, "ankastre-ss", 0.0),
             ("Pencere tüpü, girişte (ayna–2. baffle)", B_end + B, "ankastre-ss", 0.0),
             ("Merkez açıklık (baffle–baffle)", B, "ss-ss", 0.5),
             ("Pencere tüpü, ortada (2·B)", 2 * B, "ss-ss", 0.5) if Nb >= 3 else None,
             ("Pencere tüpü, çıkışta (2. son baffle–ayna)", B_end + B, "ankastre-ss", 1.0),
             ("Çıkış açıklığı (son baffle–ayna)", B_end, "ankastre-ss", 1.0)]
    spans = [s for s in spans if s]
    # Eksenel kuvvet (V-6): kabuk tüpten sıcaksa tüpler çekmede (frekansı artırır) → güvenli tarafta A = 1.
    sig_t = mech.get("sigma_tube_MPa", 0.0) * (1 if r["T_shell_metal"] >= r["T_tube_metal"] else -1)
    xt, xl = pr, (0.866 if ang == 30 else (0.5 if ang == 60 else (0.707 if ang == 45 else 1.0))) * pr
    if ang == 60:
        xt, xl = 1.732 * pr, 0.5 * pr
    S_st = strouhal if strouhal else 0.5                         # Şekil V-12.2B'den yaklaşık okuma (30°, P/do 1,25)
    CL = _interp(TEMA_CL_30, pr)
    rows = []
    for name, l_mm, cond, x in spans:
        st = _local_shell_state(r, curve, x)
        rho_o = st["rho_TP"] * _LBFT3
        Hm = Cm * 0.00545 * rho_o * do**2
        wo = wt + wfi + Hm
        l = l_mm * _IN
        C, K = TEMA_C_SPAN[cond], TEMA_K_AXIAL[cond]
        Fcr = K**2 * E * I / l**2
        F_ax = min(0.0, sig_t * 1e6 * _PSI * At)                 # yalnız basma hesaba katılır
        A_ax = math.sqrt(max(1 + F_ax / Fcr, 0.0))
        fn = 10.838 * A_ax * C / l**2 * math.sqrt(E * I / wo)
        Nsp = Nb + 1
        delta = 0.314 * (Nsp - 1) / Nsp * math.sqrt(tb * _IN / l)
        X = 144 * wo * delta / (rho_o * do**2)
        D = tema_D_factor(X, ang, pr)
        Vc = D * fn * do / 12                                    # ft/s
        L_cross = (B_end if x in (0.0, 1.0) else B) / 1000
        Sm = bd["Sm"] * L_cross / bd["B"]
        V = st["W"] / 3600 / (st["rho_TP"] * Sm) * _FT           # ft/s
        f_vs = 12 * S_st * V / do
        f_tb = 12 * V / (do * xl * xt) * (3.05 * (1 - 1 / xt) ** 2 + 0.28)
        y_vs = CL * rho_o * do * V**2 / (20.6 * fn**2 * delta * wo) if fn < 2 * f_vs else 0.0
        entrance = x in (0.0, 1.0)
        CF = (0.022 if fn <= 40 else max(-0.00045 * fn + 0.04, 0)) if entrance else (0.012 if fn <= 40 else max(-0.00025 * fn + 0.022, 0))
        y_tb = CF * rho_o * do * V**2 / (8.67 * delta**0.5 * fn**1.5 * wo)
        rows.append(dict(acıklık=name, l_mm=l_mm, uc_kosulu=cond, fn_Hz=fn, wo_lbft=wo, delta=delta, X=X, D=D,
                         V_ms=V / _FT, Vc_ms=Vc / _FT, V_Vc=V / Vc, f_vs_Hz=f_vs, f_tb_Hz=f_tb,
                         y_vs_mm=y_vs * 25.4, y_tb_mm=y_tb * 25.4, y_sinir_mm=0.02 * g["Do"],
                         rho_o=st["rho_TP"], eps_g=st["eps_g"]))
    df = pd.DataFrame(rows)
    # Akustik rezonans (V-12): yalnız gaz/buhar varken; yansıtıcı duvarlar arası mesafe ≈ kabuk iç çapı.
    inl = _local_shell_state(r, curve, 0.0)
    w_in = g["Ds"] * _IN
    Ps_psia = inl["P"] * 1e5 * _PSI
    rho_g = inl["rho_V"] * _LBFT3
    fa1 = 409 / w_in * math.sqrt(Ps_psia * gamma / (rho_g * (1 + 0.5 / (xl * xt))))
    V_in = df.V_ms.iloc[0] * _FT
    fvs_in, ftb_in = df.f_vs_Hz.iloc[0], df.f_tb_Hz.iloc[0]
    xo = 2 * xl if ang in (30, 45, 60) else xl
    Re = 124.13 * do * V_in * rho_g / (inl["mu_V"] * 1000)
    acu = []
    for n in (1, 2):
        fa = n * fa1
        A_ = (0.8 * fvs_in < fa < 1.2 * fvs_in) or (0.8 * ftb_in < fa < 1.2 * ftb_in)
        B_ = V_in > fa * do * (xl - 0.5) / 6
        C_ = V_in > fa * do / (12 * S_st) and Re / (S_st * xo) * (1 - 1 / xt) ** 2 / xt > 2000
        acu.append(dict(mod=n, fa_Hz=fa, A=A_, B=B_, C=C_))
    acu = pd.DataFrame(acu)
    worst = df.loc[df.V_Vc.idxmax()]
    ok_fei = bool((df.V_Vc < 1).all())
    ok_amp = bool(((df.y_vs_mm <= df.y_sinir_mm) & (df.y_tb_mm <= df.y_sinir_mm)).all())
    return dict(spans=df, acoustic=acu, Cm=Cm, tb_mm=tb, strouhal=S_st, V_Vc_max=float(worst.V_Vc),
                kritik_aciklik=worst["acıklık"], fn_min=float(df.fn_Hz.min()), fei_ok=ok_fei, amp_ok=ok_amp,
                margin_ok=bool(df.V_Vc.max() <= VIB_MARGIN), acoustic_risk=bool((acu.A | acu.B | acu.C).any()))


def tubesheet_thickness(r, mech, F=1.25, CA_tube_mm=3.0, S_ts_MPa=None):
    """Sabit ayna kalınlığı, TEMA Ek A.131 (eğilme) ve A.132 (kesme), RCB-7.11 alt sınırı.
    F = 1,25: sabit (her iki tarafa kaynaklı) ayna için Şekil A.131 eğri H'nin üst değeri (güvenli tarafta).
    G = kabuk iç çapı (aynanın kabuğa kaynaklı olduğu bölge). P = kabuk ve tüp tasarım basınçlarının büyüğü;
    A.15'teki ısıl genleşme kaynaklı etkin basınç ayrıca hesaplanmamış, genleşme gerilmesi mechanical()'da sınanır."""
    g, bd = r["geo"], r["bd"]
    S = (S_ts_MPa or mech.get("S_shell_MPa", 138.0)) * 1000            # kPa
    P = max(mech["P_design_shell_barg"], mech["P_design_tube_barg"]) * 100   # kPa
    pr = g["Pt"] / g["Do"]
    eta = 1 - (0.907 if bd["ang"] in (30, 60) else 0.785) / pr**2
    G = g["Ds"]
    T_bend = F * G / 3 * math.sqrt(P / (eta * S))
    DL = bd["Dctl"] * 1000                                              # 4A/C ≈ tüp merkez sınırı çapı (dairesel kabul)
    T_shear = 0.31 * DL / (1 - g["Do"] / g["Pt"]) * (P / S)
    CA = mech.get("CA_shell_mm", 3.0) + CA_tube_mm
    T_min = max(g["Do"] + CA, 19.1)                                     # R-7.11: (T − CA) ≥ do, toplam ≥ 19,1 mm
    T_req = max(max(T_bend, T_shear) + CA, T_min)
    T_sel = next(t for t in [20, 22, 25, 28, 30, 32, 35, 38, 40, 45, 50, 55, 60, 65, 70, 75, 80, 90, 100] if t >= T_req - 1e-9)
    return dict(F=F, G_mm=G, P_kPa=P, S_kPa=S, eta=eta, T_bend_mm=T_bend, T_shear_mm=T_shear, DL_mm=DL,
                CA_mm=CA, T_min_mm=T_min, T_req_mm=T_req, T_sel_mm=T_sel,
                kontrol="eğilme" if T_bend >= T_shear else "kesme")


# =====================================================================
# 6. MALİYET (Towler & Sinnott, 2. bs., Tablo 7.2; ABD Körfez Kıyısı, Ocak 2010, CEPCI = 532,9)
# =====================================================================
COST_CORRELATIONS = {           # Ce = a + b·S^n,  S = alan m², geçerli 10–1000 m²
    "U-tüp (sabit ayna için de yaklaşık)": (28000.0, 54.0, 1.2),
    "Yüzer kafa": (32000.0, 70.0, 1.2),
}
MATERIAL_FACTORS = {"Karbon çeliği": 1.0, "Paslanmaz 304": 1.3, "Paslanmaz 316": 1.3, "Paslanmaz 321": 1.5,
                    "Hastelloy C": 1.55, "Monel": 1.65, "Nikel / Inconel": 1.7}
INSTALL_FACTORS = dict(f_er=0.3, f_p=0.8, f_i=0.3, f_el=0.2, f_c=0.3, f_s=0.2, f_l=0.1)   # akışkan prosesleri


def capital_recovery_factor(i, n):
    return i * (1 + i) ** n / ((1 + i) ** n - 1)


def cost_estimate(r, cost):
    """Yatırım (satın alma + kurulu) ve işletme maliyeti, yıllık toplam maliyet.
    cost: correlation, f_m, cepci_now, cepci_base, install(dict), i_rate, years, hours,
          elec_price, pump_eff, water_price_m3, fx_rate"""
    a, b, n = COST_CORRELATIONS[cost["correlation"]]
    shells = r.get("shells", 1)
    S = r["A"] / shells                                               # korelasyon tek kabuk alanı için
    Ce_2010 = shells * (a + b * S**n)
    Ce = Ce_2010 * cost["cepci_now"] / cost["cepci_base"]
    f = cost["install"]
    install_mult = (1 + f["f_p"]) * cost["f_m"] + (f["f_er"] + f["f_el"] + f["f_i"] + f["f_c"] + f["f_s"] + f["f_l"])
    C_installed = Ce * install_mult
    crf = capital_recovery_factor(cost["i_rate"], cost["years"])
    w = r["water"]
    V_w = w["m"] / 3600 / w["rho"]                                   # m³/s
    P_pump_kW = r["dp_tube_kPa"] * V_w / cost["pump_eff"]             # kPa·m³/s = kW
    E_cost = P_pump_kW * cost["hours"] * cost["elec_price"]
    W_cost = V_w * 3600 * cost["hours"] * cost["water_price_m3"]
    annual_capital = C_installed * crf
    tac = annual_capital + E_cost + W_cost
    return dict(S=S, a=a, b=b, n=n, Ce_2010=Ce_2010, Ce=Ce, install_mult=install_mult, C_installed=C_installed,
                crf=crf, annual_capital=annual_capital, P_pump_kW=P_pump_kW, V_w_m3h=V_w * 3600,
                E_cost=E_cost, W_cost=W_cost, TAC=tac, in_range=10 <= S <= 1000, fx=cost.get("fx_rate"))


# =====================================================================
# 7. DUYARLILIK, DOĞRULAMA, TARAMA, KRİTER KONTROLÜ
# =====================================================================
def sensitivity(tube, geo, curve, base, **kw):
    """Tek parametre değişimlerinin gerekli alan ve aşırı tasarıma etkisi."""
    cases = [("Temel durum", {}, {}),
             ("Kabuk h × 0,8", dict(h_mult=0.8 * kw.get("h_mult", 1.0)), {}),
             ("Kabuk h × 1,2", dict(h_mult=1.2 * kw.get("h_mult", 1.0)), {}),
             ("Kabuk kirlenme × 2", dict(Rf_shell=2 * kw["Rf_shell"]), {}),
             ("Tüp (su) kirlenme × 2", {}, dict(Rf=2 * tube["Rf"])),
             ("Su çıkışı −2 °C (daha çok su)", {}, dict(T_out=tube["T_out"] - 2, m=None)),
             ("Su çıkışı +2 °C (daha az su)", {}, dict(T_out=tube["T_out"] + 2, m=None))]
    rows = []
    for name, kw_mod, tube_mod in cases:
        try:
            r = rate_condensing(dict(tube, **tube_mod), geo, curve, **dict(kw, **kw_mod))
            rows.append(dict(Durum=name, A_gerekli_m2=r["A_req"], degisim_yuzde=(r["A_req"] / base["A_req"] - 1) * 100,
                             asiri_tasarim_yuzde=r["over_design"], su_debisi_kgh=r["water"]["m"],
                             tup_hizi_ms=r["ts"]["v"], tup_dP_kPa=r["dp_tube_kPa"]))
        except Exception as e:
            rows.append(dict(Durum=name, A_gerekli_m2=np.nan, degisim_yuzde=np.nan, asiri_tasarim_yuzde=np.nan,
                             su_debisi_kgh=np.nan, tup_hizi_ms=np.nan, tup_dP_kPa=np.nan, hata=str(e)))
    return pd.DataFrame(rows)


HT_ANGLE = {"triangular": 30, "rotated square": 45, "square": 90}


def hydraulic_verification(tube, geo, curve, r, rate_kw):
    """Hidrolik sonuçların bağımsız kontrolü (uygulamanın kendi denklemlerinden ayrı kaynaklarla):
      tüp sayısı   : geometrik ayna sayımı ↔ Phadke (1984, tam sayım), Perry, VDI ve Towler K1–n1 bağıntıları
      tüp ΔP       : aynı girdilerle el hesabı; Blasius (pürüzsüz) ve Kern (geçiş başına 4 hız yükü) alternatifleri
      Bell–Delaware: Jc, Jl, Jb, Js — Taborek denklemleri ↔ ht kütüphanesindeki Taborek grafiklerinin sayısallaştırması
      kabuk ΔP     : Bell–Delaware ↔ Kern (tek faz gaz); iki faz sabiti C ve kabuk yöntemi duyarlılığı
    Döndürür: {tablo adı: DataFrame}."""
    from ht.hx import Ntubes_Phadkeb, Ntubes_Perrys, Ntubes_VDI
    import ht.conv_tube_bank as tbk
    out = {}
    sg, ts, w = r["sg"], r["ts"], r["water"]
    Do, Pt, Db, ang = sg["Do"], sg["Pt"], sg["Db"], HT_ANGLE.get(geo["layout"], 30)

    # --- 1) tüp sayısı
    def _try(f, *a):
        try:
            return int(f(*a))
        except Exception:
            return np.nan
    rows = []
    for p in sorted({1, 2, 4, 6, 8} & set(K1N1["triangular"])):
        app = tube_layout(Db * 1000, geo["Do"], geo["Pt"], geo["layout"], p)[0]
        K1, n1 = K1N1["triangular" if geo["layout"] == "triangular" else "square"][p]
        ph = _try(Ntubes_Phadkeb, Db, Do, Pt, p, ang)
        rows.append({"Geçiş": p, "Uygulama (ayna sayımı)": app, "Phadke (tam sayım)": ph,
                     "Perry": _try(Ntubes_Perrys, Db, Do, p, ang), "VDI": _try(Ntubes_VDI, Db, p, Do, Pt, ang),
                     "Towler K1·(Db/Do)^n1": int(K1 * (Db / Do) ** n1),
                     "Fark Phadke'ye göre (%)": (app / ph - 1) * 100 if ph == ph else np.nan,
                     "": "◀ bu tasarım" if p == geo["passes"] else ""})
    out["tup_sayisi"] = pd.DataFrame(rows)

    # --- 2) tüp tarafı ΔP (el hesabı, yalnız temel denklemlerle)
    Di, L, Np, N = ts["Di"], sg["L"], geo["passes"], geo["N"]
    v = w["m"] / 3600 / (w["rho"] * math.pi / 4 * Di**2 * N / Np)
    Re = w["rho"] * v * Di / w["mu"]
    hv = w["rho"] * v**2 / 2
    f_col = fluids.friction_factor(Re, eD=0.000045 / Di)
    f_bla = 0.316 * Re**-0.25
    sh = r["shells"]
    noz = pick_nozzle(w["m"], w["rho"], RHOV2_TUBE_NOZZLE)
    dp_noz = (1.0 + 0.5) * noz["rhov2"] / 2 / 1000
    tube_rows = [
        ("Uygulama", "Colebrook, geçiş başına K hız yükü", ts["dp_kPa"]),
        ("El hesabı (aynı yöntem)", f"v = {v:.3f} m/s, Re = {Re:.0f}, f = {f_col:.4f}, K = {ts['K_pass']:g}",
         sh * Np * (f_col * L / Di + ts["K_pass"]) * hv / 1000),
        ("Pürüzsüz boru (Blasius)", f"f = 0,316·Re^-0,25 = {f_bla:.4f}", sh * Np * (f_bla * L / Di + ts["K_pass"]) * hv / 1000),
        ("Kern (geçiş başına 4 hız yükü)", "Kern (1950): dönüş kaybı 4·ρv²/2", sh * Np * (f_col * L / Di + 4.0) * hv / 1000),
        ("Uygulama + nozullar", f"NPS {noz['NPS']}\", giriş 1,0 + çıkış 0,5 hız yükü = {dp_noz:.2f} kPa", ts["dp_kPa"] + dp_noz),
    ]
    t = pd.DataFrame(tube_rows, columns=["Yöntem", "Ayrıntı", "Tüp ΔP (kPa)"])
    t["Uygulamaya göre fark (%)"] = (t["Tüp ΔP (kPa)"] / ts["dp_kPa"] - 1) * 100
    out["tup_dp"] = t
    out["tup_dp_parca"] = dict(surtunme=sh * Np * f_col * L / Di * hv / 1000, donus=sh * Np * ts["K_pass"] * hv / 1000)

    # --- 3) Bell–Delaware düzeltme çarpanları (kabuk girişi, tek faz gaz)
    a = curve.iloc[0]
    mV = a.m_V / 3600
    bd = r.get("bd") or bd_geometry(dict(geo, Db=Db * 1000))
    B = bell_delaware(mV, a.rho_V, a.mu_V, a.cp_V, a.k_V, bd)

    def _ht(f, *a_, **k_):
        try:
            return float(f(*a_, **k_))
        except Exception:
            return np.nan
    jr = [("Jc", "Pencere (baffle kesimi)", B["Jc"], _ht(tbk.baffle_correction_Bell, bd["Fc"])),
          ("Jl", "Baffle kaçakları", B["Jl"], _ht(tbk.baffle_leakage_Bell, bd["Ssb"], bd["Stb"], bd["Sm"])),
          ("Jb", "Demet bypass'ı", B["Jb"], _ht(tbk.bundle_bypassing_Bell, bd["Fsbp"], bd["n_ss"], bd["Ntcc"], laminar=B["Re"] < 100)),
          ("Js", "Uç bölme aralıkları", B["Js"], _ht(tbk.unequal_baffle_spacing_Bell, bd["Nb"], bd["B"], bd["B_in"], bd["B_out"],
                                                    laminar=B["Re"] < 100))]
    j = pd.DataFrame(jr, columns=["Çarpan", "Etki", "Uygulama (Taborek denklemi)", "ht (Taborek grafiği)"])
    j["Fark (%)"] = (j["Uygulama (Taborek denklemi)"] / j["ht (Taborek grafiği)"] - 1) * 100
    out["bd_carpanlar"] = j
    out["bd_giris"] = dict(Re=B["Re"], Sm=bd["Sm"], Ntcc=bd["Ntcc"], Nb=bd["Nb"])

    # --- 4) kabuk ΔP: yöntem ve iki faz sabiti duyarlılığı
    kern1 = kern_dp(mV, a.rho_V, a.mu_V, sg)
    rows = [("Tek faz gaz (giriş koşulu), Bell–Delaware", B["dP"] / 1000),
            ("Tek faz gaz (giriş koşulu), Kern", kern1 / 1000)]
    for C in (0.0, 5.0, 8.0, 12.0, 20.0):
        rr = rate_condensing(tube, geo, curve, **dict(rate_kw, two_phase_C=C, shell_method="bell_delaware"))
        rows.append((f"Tasarım, Bell–Delaware, C = {C:g}" + (" (uygulama)" if C == TWO_PHASE_C else ""), rr["dp_shell_kPa"]))
    rk = rate_condensing(tube, geo, curve, **dict(rate_kw, shell_method="kern"))
    rows.append(("Tasarım, Kern, C = 8", rk["dp_shell_kPa"]))
    rows.append(("En kötü durum: gaz hiç yoğuşmuyor", worst_case_shell_dp(r)))
    out["kabuk_dp"] = pd.DataFrame(rows, columns=["Durum", "Kabuk ΔP (kPa)"])
    out["kabuk_dp_parca"] = dict(B_capraz=B["dPc"] / 1000, B_pencere=B["dPw"] / 1000, B_uc=B["dPe"] / 1000)
    return out


def size_scan(rate_fn, tube_sizes=((19.05, 2.11), (25.4, 1.65)), lengths=(3048, 4877, 6096),
              passes_list=(2, 4, 6), layout="triangular", shell_diameters=range(500, 1401, 50),
              baffle_fractions=(0.3, 0.4, 0.5, 0.6, 0.8, 1.0), bundle_clearance=15.0,
              v_range=(1.0, 2.5), dp_tube_max=35.0, dp_shell_max=20.0, od_range=(10.0, 30.0),
              cost_fn=None, top=15, water=None):
    """Geometri taraması. rate_fn(geo) -> sonuç. cost_fn(r) verilirse maliyetler de hesaplanır.
    water=(debi kg/h, ρ) verilirse tüp hızı sınır dışı olan geometriler rating yapılmadan elenir (hız)."""
    rows, counts = [], {}
    for (Do, t), L, Np, Ds in ((a, b, c, d) for a in tube_sizes for b in lengths
                               for c in passes_list for d in shell_diameters):
        key = (Do, Np, Ds)
        if key not in counts:                                         # kabuğa sığan tüp: geometrik ayna sayımı
            counts[key] = tube_layout(Ds - bundle_clearance, Do, 1.25 * Do, layout, Np)[0] // Np * Np
        N = counts[key]
        if N < 10:
            continue
        if water:                                                    # hız ön-elemesi
            Di = (Do - 2 * t) / 1000
            v = water[0] / 3600 / (water[1] * math.pi / 4 * Di**2 * N / Np)
            if not v_range[0] <= v <= v_range[1]:
                continue
        for frac in baffle_fractions:
            g = dict(N=N, passes=Np, Do=Do, t=t, layout=layout, Pt=1.25 * Do, L=L, B=round(frac * Ds), Ds=Ds,
                     Db=Ds - bundle_clearance)
            span_max = TEMA_MAX_SPAN.get(Do, np.inf)
            if N < 10 or 2 * g["B"] > span_max:
                continue
            try:
                r = rate_fn(g)
            except Exception:
                continue
            if (v_range[0] <= r["ts"]["v"] <= v_range[1] and r["lf"]["F"] >= 0.75
                    and od_range[0] <= r["over_design"] <= od_range[1]
                    and r["dp_tube_kPa"] <= dp_tube_max and r["dp_shell_kPa"] <= dp_shell_max):
                row = dict(Do=Do, L=L, gecis=Np, Ds=Ds, B=g["B"], N=N, Db=g["Db"], A=r["A"], A_gerekli=r["A_req"],
                           asiri_tasarim=r["over_design"], v_tup=r["ts"]["v"], dP_tup=r["dp_tube_kPa"],
                           dP_kabuk=r["dp_shell_kPa"], h_kabuk=r["h_shell"])
                if cost_fn:
                    c = cost_fn(r)
                    row.update(yatirim_USD=c["C_installed"], yillik_sermaye_USD=c["annual_capital"],
                               elektrik_USD=c["E_cost"], su_USD=c["W_cost"], yillik_toplam_USD=c["TAC"])
                rows.append(row)
    if not rows:
        return pd.DataFrame()
    df = pd.DataFrame(rows).sort_values("yillik_toplam_USD" if cost_fn else "A").reset_index(drop=True)
    return df if top is None else df.head(top)


def checks(r, v_range, dp_tube_max, dp_shell_max, od_range, mech_res=None, clearance=None):
    g = r["geo"]
    items = [
        ("F ≥ 0,75", f"{r['lf']['F']:.3f}", r["lf"]["F"] >= 0.75),
        (f"Tüp hızı {v_range[0]}–{v_range[1]} m/s", f"{r['ts']['v']:.2f}", v_range[0] <= r["ts"]["v"] <= v_range[1]),
        ("Tüp Re > 10 000 (türbülanslı)", f"{r['ts']['Re']:.0f}", r["ts"]["Re"] > 10000),
        ("Tüp sayısı geçişe tam bölünür", f"{g['N']} / {g['passes']}", g["N"] % g["passes"] == 0),
        ("Baffle aralığı 0,2·Ds – 1,0·Ds", f"{g['B']:.0f} mm", 0.2 * g["Ds"] <= g["B"] <= g["Ds"]),
        (f"Aşırı tasarım %{od_range[0]:.0f}–%{od_range[1]:.0f}", f"%{r['over_design']:.1f}",
         od_range[0] <= r["over_design"] <= od_range[1]),
        (f"Tüp ΔP ≤ {dp_tube_max:g} kPa", f"{r['dp_tube_kPa']:.1f}", r["dp_tube_kPa"] <= dp_tube_max),
        (f"Kabuk ΔP ≤ {dp_shell_max:g} kPa", f"{r['dp_shell_kPa']:.1f}", r["dp_shell_kPa"] <= dp_shell_max),
    ]
    if mech_res:
        items += [("Tüp destek açıklığı ≤ TEMA sınırı", f"{mech_res['span_mm']:.0f} / {mech_res['span_max_mm']:.0f} mm",
                   mech_res["span_ok"]),
                  ("Tüp et kalınlığı yeterli (iç basınç)", f"{mech_res['t_tube_calc']:.2f} ≤ {g['t']} mm",
                   mech_res["t_tube_ok"])]
        if "vib" in mech_res:
            v = mech_res["vib"]
            items += [("Titreşim: V < Vc her açıklıkta (TEMA Bölüm 6, V-10)", f"en büyük V/Vc = {v['V_Vc_max']:.2f}", v["fei_ok"]),
                      (f"Titreşim payı V/Vc ≤ {VIB_MARGIN}", f"{v['V_Vc_max']:.2f} ({v['kritik_aciklik']})", v["margin_ok"])]
        if "tubesheet" in mech_res:
            t_ = mech_res["tubesheet"]
            items.append(("Ayna kalınlığı (TEMA A.131/A.132, R-7.11)", f"gerekli {t_['T_req_mm']:.1f} → seçilen {t_['T_sel_mm']} mm", True))
    if clearance is not None:
        n_fit = tube_layout(g["Ds"] - clearance, g["Do"], g["Pt"], g["layout"], g["passes"])[0]
        items.append(("Tüpler kabuğa sığıyor (geometrik sayım)", f"{g['N']} / en çok {n_fit}", g["N"] <= n_fit))
    return pd.DataFrame([(a, b, "✅ Uygun" if ok else "❌ Uygun değil") for a, b, ok in items],
                        columns=["Kriter", "Değer", "Durum"])


# =====================================================================
# 8. DWSIM KARŞILAŞTIRMASI
# =====================================================================
# (büyüklük, birim, tolerans, tolerans tipi, fark büyükse olası neden)
DWSIM_ITEMS = [
    ("Q", "Isı yükü", "kW", 3.0, "%", "Farklı termodinamik paket (ideal ↔ NRTL) veya farklı ideal gaz Cp / buharlaşma entalpisi verisi."),
    ("T_dew", "Çiy noktası", "°C", 2.0, "abs", "Sıvı modeli/ikili parametreler farklı; aynı paketi seçin."),
    ("m_V_out", "Çıkış buhar debisi", "kg/h", 5.0, "%", "Faz dengesi farklı (aktivite modeli, H₂ çözünürlüğü)."),
    ("product_V_out", "Çıkış buharındaki ürün", "kg/h", 10.0, "%", "Faz dengesi farklı; ürün kaybı VLE'ye çok duyarlıdır."),
    ("dp_tube", "Tüp ΔP", "kPa", 30.0, "%", "Dönüş kaybı katsayısı (2,5 ↔ 4 hız yükü) ve sürtünme bağıntısı farkı."),
    ("U", "U (kirli)", "W/m²K", 30.0, "%", "Beklenen fark: DWSIM basitleştirilmiş Tinker (ortalama özellik), "
                                            "bu uygulama bölgesel Nusselt + Colburn–Hougen kullanır."),
    ("dp_shell", "Kabuk ΔP", "kPa", 50.0, "%", "Beklenen fark: DWSIM Tinker akım analizi (tek faz) ↔ bu uygulama Bell–Delaware + iki faz çarpanı (C = 8)."),
]


def dwsim_bounds(r, hv):
    """DWSIM'in U ve ΔP sonuçları için fiziksel aralık. DWSIM kabukta yoğuşmayı ve H₂ direncini ayrıca modellemez,
    dönüş kayıplarını farklı sayar; bu yüzden doğru bir DWSIM sonucu uygulamanın iki uç durumu arasında kalmalıdır:
      U       : yalnız gaz konveksiyonu (yoğuşmasız)  ↔  bölgesel yoğuşmalı U (uygulama)
      tüp ΔP  : yalnız sürtünme                      ↔  sürtünme + geçiş başına 4 hız yükü (Kern)
      kabuk ΔP: iki fazlı (uygulama)                  ↔  gaz hiç yoğuşmuyor (en kötü durum)"""
    z, ts, w = r["zones"], r["ts"], r["water"]
    Ug = [overall_U(h, ts["h"], r["Rf_shell"], w["Rf"], r["sg"]["Do"], ts["Di"], r["k_wall"])[1] for h in z.h_G]
    U_gas = float((np.array(Ug) * z.A_m2).sum() / z.A_m2.sum())
    U_app = float((z.Ud * z.A_m2).sum() / r["A_req"])
    kern = float(hv["tup_dp"].loc[hv["tup_dp"]["Yöntem"].str.startswith("Kern"), "Tüp ΔP (kPa)"].iloc[0])
    return {"U": (U_gas, U_app), "dp_tube": (hv["tup_dp_parca"]["surtunme"], kern),
            "dp_shell": (r["dp_shell_kPa"], worst_case_shell_dp(r))}


BOUND_MARGIN = 0.10     # aralık uçlarında ±%10 pay (iki programın fiziksel özellik verisi farkı)


def dwsim_compare(app_vals, dwsim_vals, bounds=None):
    """Termodinamik satırlar toleransla; U ve ΔP satırları (bounds verilirse) fiziksel aralıkla karşılaştırılır.
    Boş (None/0) DWSIM değerleri atlanır."""
    rows = []
    for key, name, unit, tol, kind, reason in DWSIM_ITEMS:
        a, d = app_vals.get(key), dwsim_vals.get(key)
        if a is None or not d:
            continue
        diff = (a - d) if kind == "abs" else (a / d - 1) * 100
        if bounds and key in bounds:
            lo, hi = sorted(bounds[key])
            ok = lo * (1 - BOUND_MARGIN) <= d <= hi * (1 + BOUND_MARGIN)
            crit = f"aralık {lo:.1f}–{hi:.1f}"
            why = ("DWSIM değeri uygulamanın iki uç durumunun dışında: girdileri (geometri, debi, kirlenme) iki programda "
                   "karşılaştırın.")
        else:
            ok = abs(diff) <= tol
            crit = f"±{tol:g} {'°C' if kind == 'abs' else '%'}"
            why = reason
        rows.append({"Büyüklük": name, "Birim": unit, "Uygulama": a, "DWSIM": d,
                     "Fark": f"{diff:+.1f} {'°C' if kind == 'abs' else '%'}", "Kabul ölçütü": crit,
                     "Durum": "✅ Uyumlu" if ok else "❌ Uyumsuz", "Açıklama": "" if ok else why})
    return pd.DataFrame(rows)


# =====================================================================
# 9. OPTİMİZASYON (en düşük yıllık toplam maliyet)
# =====================================================================
def optimize(tube, curve, water_T_outs, cost_fn, rate_kw, limits, tube_sizes=((19.05, 2.11), (25.4, 1.65)),
             progress=None, clearance=15.0):
    """Karar değişkenleri: soğutma suyu çıkış sıcaklığı + geometri (Do, L, geçiş, Ds, B).
    Her su çıkış sıcaklığı için tüm geometriler taranır; kısıtları sağlayanlar maliyetleriyle döner.
    Döndürür: (tüm uygun adaylar, her sıcaklıktaki en iyi aday)."""
    Q = float(curve["Q_cum"].iloc[-1])
    frames = []
    for k, T_out in enumerate(water_T_outs):
        tb = dict(tube, T_out=T_out, m=None)
        try:
            w = water_side_balance(tb, Q)
        except Exception:
            continue
        df = size_scan(lambda g: rate_condensing(tb, g, curve, **rate_kw), tube_sizes=tube_sizes,
                       v_range=limits["v"], dp_tube_max=limits["dp_t"], dp_shell_max=limits["dp_s"],
                       od_range=limits["od"], cost_fn=cost_fn, top=None, water=(w["m"], w["rho"]),
                       bundle_clearance=clearance)
        if not df.empty:
            df.insert(0, "T_su_cikis", T_out)
            df["su_debisi_kgh"] = w["m"]
            frames.append(df)
        if progress:
            progress((k + 1) / len(water_T_outs), T_out)
    if not frames:
        return pd.DataFrame(), pd.DataFrame()
    allc = pd.concat(frames, ignore_index=True)
    # aynı eşanjörün yalnız baffle aralığı farklı kopyaları: kabuk ΔP'si en düşük olan tutulur
    allc = (allc.sort_values(["yillik_toplam_USD", "dP_kabuk"])
                .drop_duplicates(["T_su_cikis", "Do", "L", "gecis", "Ds", "N"]).reset_index(drop=True))
    best = allc.loc[allc.groupby("T_su_cikis")["yillik_toplam_USD"].idxmin()].sort_values("T_su_cikis")
    return allc, best.reset_index(drop=True)


def robustness(tube, curve, row, rate_kw, h_factor=0.8):
    """Aday geometride kabuk h'si h_factor kadar düşük çıkarsa aşırı tasarım ne olur?"""
    g = dict(N=int(row.N), passes=int(row.gecis), Do=float(row.Do), t=2.11 if float(row.Do) < 20 else 1.65,
             layout="triangular", Pt=round(1.25 * float(row.Do), 2), L=float(row.L), B=float(row.B), Ds=float(row.Ds),
             Db=float(row.Db) if "Db" in row and pd.notna(row.Db) else None)
    tb = dict(tube, T_out=float(row.T_su_cikis), m=None)
    kw = dict(rate_kw, h_mult=rate_kw.get("h_mult", 1.0) * h_factor)
    return rate_condensing(tb, g, curve, **kw)["over_design"]


TECH_CRITERIA = [
    ("K1", "Aşırı tasarım %20–35"),
    ("K2", "Kabuk h %20 düşükken pay ≥ %10"),
    ("K3", "Su 5 °C sıcak gelirse (aynı debi) pay ≥ %0"),
    ("K4", "Kabuk ΔP ≤ sınırın %85'i"),
    ("K5", "Tüp ΔP ≤ sınırın %85'i"),
    ("K6", "Tüp hızı 1,2–2,0 m/s"),
    ("K7", "Desteksiz açıklık ≤ TEMA sınırının %85'i"),
    ("K8", "F ≥ 0,80"),
    ("K9", "Kabuk ΔP, gaz hiç yoğuşmasa bile ≤ sınır"),
    ("K10", "Titreşim: her açıklıkta V < Vc (akışkan-elastik kararsızlık yok)"),
    ("K11", "Titreşimde %20 pay: V/Vc ≤ 0,8"),
]
# Kritik ölçütler: hesabın belirsiz olduğu noktalarda tasarımın yine çalışmasını sınar. TEMA RGP-T-2.4,
# belirsizliklerin kirlenme direncine gizlenmeden tasarım hesabında açıkça gösterilmesini ister.
# K10: TEMA Bölüm 6 titreşim değerlendirmesini zorunlu tutmaz (V-1.1: "istendiğinde"; V-10 "should"), ancak
# akışkan-elastik kararsızlık tüp hasarının başlıca nedenidir; sektör şartnameleri ve tasarım programları bu kontrolü yapar.
CRITICAL_CRITERIA = ("K2", "K3", "K9", "K10")
RHOV2_BUNDLE_MAX = 5953.0       # kg/m·s², TEMA RCB-4.62: demet/kabuk giriş-çıkış alanında hiçbir durumda aşılamaz


def worst_case_shell_dp(r):
    """Kabuk ΔP üst sınırı: gaz hiç yoğuşmuyor ve tüm boy boyunca giriş debisi ve özellikleriyle akıyor
    (çekirdek + iki uç bölmesi + nozullar)."""
    z = r["zones"]
    return float(z.dP_core_2f_kPa.iloc[0] + 2 * z.dP_end_2f_kPa.iloc[0] + r["dp_shell_parts"]["nozullar"])


def bundle_rhov2(r, curve):
    """TEMA RCB-4.624: darbe plakası yokken demet giriş/çıkış alanı = baffle bölmesindeki çapraz akış alanı.
    Çıkışta iki faz homojen yoğunluk. Döndürür (giriş, çıkış) ρv², kg/(m·s²)."""
    inl, out = curve.iloc[0], curve.iloc[-1]
    G = (inl.m_V + inl.m_L) / 3600 / r["sg"]["As"]
    x = out.m_V / (out.m_V + out.m_L)
    rho_out = 1 / (x / out.rho_V + (1 - x) / out.rho_L) if out.m_L > 0 else out.rho_V
    return G**2 / inl.rho_V, G**2 / rho_out


def technical_evaluation(tube, curve, cands, rate_kw, limits, progress=None):
    """Teknik uygunluk değerlendirmesi, üç kademede:
    1) TEMA zorunlu hükümleri: desteksiz açıklık ≤ RCB-4.52 sınırı, demet giriş/çıkış ρv² ≤ 5953 (RCB-4.62).
       Titreşimde V/Vc > 0,8 çıkarsa baffle aralığı daraltılarak yeniden denenir (TEMA V-13.2).
    2) Kritik ölçütler (K2, K3, K9, K10): belirsizlik ve işletme değişimlerine dayanıklılık.
    3) Tüm ölçütlerin toplam puanı (K1–K11).
    Sıralama: TEMA zorunlu → kritik puan → toplam puan → en düşük yıllık maliyet."""
    T_h1, T_h2 = float(curve["T"].iloc[0]), float(curve["T"].iloc[-1])
    out = []
    for k, (_, row) in enumerate(cands.iterrows()):
        Do = float(row.Do)
        g = dict(N=int(row.N), passes=int(row.gecis), Do=Do, t=2.11 if Do < 20 else 1.65, layout="triangular",
                 Pt=round(1.25 * Do, 2), L=float(row.L), B=float(row.B), Ds=float(row.Ds),
                 Db=float(row.Db) if "Db" in row and pd.notna(row.Db) else None)
        T_in, T_out = float(row.T_su_giris), float(row.T_su_cikis)
        tb = dict(tube, T_in=T_in, T_out=T_out, m=None)
        r0 = rate_condensing(tb, g, curve, **rate_kw)
        vib = vibration_check(r0, curve)
        B_orig = g["B"]
        if vib["V_Vc_max"] > VIB_MARGIN:
            # TEMA V-13.2: desteksiz açıklık titreşimi en çok etkileyen etkendir → baffle aralığı daraltılarak denenir.
            # Önce V/Vc ≤ 0,8 sağlayan en geniş aralık, yoksa V/Vc < 1 sağlayan; ΔP ve aşırı tasarım sınırları korunur.
            fallback = None
            for frac in (0.8, 0.6, 0.5, 0.4, 0.3):
                Bn = round(frac * g["Ds"])
                if Bn >= g["B"] or Bn < 0.2 * g["Ds"] or Bn < 50:
                    continue
                try:
                    gn = dict(g, B=Bn)
                    rn = rate_condensing(tb, gn, curve, **rate_kw)
                except Exception:
                    continue
                if not (rn["dp_shell_kPa"] <= limits["dp_s"] and limits["od"][0] <= rn["over_design"] <= limits["od"][1]):
                    continue
                vn = vibration_check(rn, curve)
                if vn["V_Vc_max"] <= VIB_MARGIN:
                    fallback = (gn, rn, vn)
                    break
                if vn["V_Vc_max"] < 1 and fallback is None:
                    fallback = (gn, rn, vn)
            if fallback:
                g, r0, vib = fallback
                row = row.copy()
                row["B"], row["dP_kabuk"], row["asiri_tasarim"], row["h_kabuk"] = g["B"], r0["dp_shell_kPa"], r0["over_design"], r0["h_shell"]
        dp_worst = worst_case_shell_dp(r0)
        rv_in, rv_out = bundle_rhov2(r0, curve)
        od_low = rate_condensing(tb, g, curve, **dict(rate_kw, h_mult=rate_kw.get("h_mult", 1.0) * 0.8))["over_design"]
        m_w = water_side_balance(tb, float(curve["Q_cum"].iloc[-1]))["m"]
        try:
            od_warm = rate_condensing(dict(tube, T_in=T_in + 5, T_out=None, m=m_w), g, curve, **rate_kw)["over_design"]
        except Exception:
            od_warm = -100.0
        F = lmtd_F(T_h1, T_h2, T_in, T_out, rate_kw.get("shells", 1))["F"]
        span = 2 * float(row.B)
        span_max = TEMA_MAX_SPAN.get(Do, np.inf)
        tema_ok = span <= span_max and max(rv_in, rv_out) <= RHOV2_BUNDLE_MAX
        crit = dict(K1=20 <= row.asiri_tasarim <= 35, K2=od_low >= 10, K3=od_warm >= 0,
                    K4=row.dP_kabuk <= 0.85 * limits["dp_s"], K5=row.dP_tup <= 0.85 * limits["dp_t"],
                    K6=1.2 <= row.v_tup <= 2.0, K7=span <= 0.85 * span_max, K8=F >= 0.80,
                    K9=dp_worst <= limits["dp_s"], K10=vib["fei_ok"], K11=vib["margin_ok"])
        out.append(dict(B_yeni=row.B, dP_kabuk_yeni=row.dP_kabuk, asiri_tasarim_yeni=row.asiri_tasarim, h_kabuk_yeni=row.h_kabuk,
                        B_ilk=B_orig, V_Vc=vib["V_Vc_max"], titresim_kritik_aciklik=vib["kritik_aciklik"],
                        fn_min_Hz=vib["fn_min"], asiri_tasarim_h_dusuk=od_low, asiri_tasarim_su_sicak=od_warm, F=F, aciklik_mm=span,
                        dP_kabuk_en_kotu=dp_worst, rhov2_demet_giris=rv_in, rhov2_demet_cikis=rv_out, tema_zorunlu=tema_ok,
                        **crit, kritik=sum(crit[c] for c in CRITICAL_CRITERIA), puan=sum(crit.values())))
        if progress:
            progress((k + 1) / len(cands))
    ev = pd.concat([cands.reset_index(drop=True), pd.DataFrame(out)], axis=1)
    for c in ("B", "dP_kabuk", "asiri_tasarim", "h_kabuk"):          # titreşim için baffle aralığı daraltıldıysa
        ev[c] = ev.pop(f"{c}_yeni")
    return ev.sort_values(["tema_zorunlu", "kritik", "puan", "yillik_toplam_USD"],
                          ascending=[False, False, False, True]).reset_index(drop=True)


def tac_vs_tubes_sweep(tube, curve, Do, t, passes, lengths, cost_fn, rate_kw, limits, clearance=15.0,
                       B_ratio=0.5, n_points=22, N_span=(0.45, 1.9), N_ref=None, layout="triangular"):
    """Tüp sayısı taraması: çap, geçiş ve su sıcaklığı sabit; her tüp boyu için N değiştirilir.
    Kabuk çapı geometrik sayımdan (OTL + boşluk), baffle aralığı B_ratio·Ds. Her noktada TAC ve uygunluk."""
    Pt = round(1.25 * Do, 2)
    span_max = TEMA_MAX_SPAN.get(Do, np.inf)
    rows = []
    for L in lengths:
        A_tube = math.pi * Do / 1000 * L / 1000
        N_mid = N_ref * 3048 / L if N_ref else 150 / A_tube
        grid = np.append(np.linspace(*N_span, n_points) * N_mid, N_mid + np.arange(-3, 4) * passes * 2)  # optimum çevresi sık
        for N in np.unique((grid // passes * passes).astype(int)):
            if N < 20:
                continue
            otl = otl_for_tubes(N, Do, Pt, layout, passes)
            Ds0 = math.ceil((otl + clearance) / 50) * 50                 # optimizasyonla aynı 50 mm kabuk adımı
            r = None
            for Ds in (Ds0, Ds0 + 50, Ds0 + 100):                         # ΔP sağlanmazsa kabuk bir adım büyütülür
                for frac in (0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0):    # kabuk ΔP'yi sağlayan en sık baffle
                    B = int(min(max(frac * Ds, 0.2 * Ds, 50), span_max / 2))
                    g = dict(N=int(N), passes=passes, Do=Do, t=t, layout=layout, Pt=Pt, L=L, B=B, Ds=Ds, Db=otl)
                    try:
                        r = rate_condensing(tube, g, curve, **rate_kw)
                    except Exception:
                        r = None
                        break
                    if r["dp_shell_kPa"] <= limits["dp_s"]:
                        break
                if r is None or r["dp_shell_kPa"] <= limits["dp_s"]:
                    break
            if r is None:
                continue
            c = cost_fn(r)
            why = []
            if r["over_design"] < limits["od"][0]:
                why.append("alan yetersiz")
            if not limits["v"][0] <= r["ts"]["v"] <= limits["v"][1]:
                why.append("tüp hızı")
            if r["dp_tube_kPa"] > limits["dp_t"]:
                why.append("tüp ΔP")
            if r["dp_shell_kPa"] > limits["dp_s"]:
                why.append("kabuk ΔP")
            rows.append(dict(L=L, N=int(N), Ds=Ds, B=B, A=r["A"], asiri_tasarim=r["over_design"], v_tup=r["ts"]["v"],
                             TAC=c["TAC"], yatirim=c["C_installed"], isletme=c["E_cost"] + c["W_cost"],
                             uygun=not why, neden=", ".join(why)))
    return pd.DataFrame(rows)
