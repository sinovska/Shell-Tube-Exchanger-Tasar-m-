"""
Optimizasyon grafikleri — tümü çizgi grafik, ortak stil.
Her fonksiyon bir matplotlib Figure döndürür (Streamlit'te st.pyplot ile gösterilir, PNG olarak da kaydedilebilir).
"""
import numpy as np
import matplotlib.pyplot as plt

# Renkler (renk körlüğü ve kontrast kontrolünden geçmiş sıra)
BLUE, GREEN, ORANGE = "#2a78d6", "#1baf7a", "#eb6834"
INK, MUTED, GRID, SHADE = "#2b2b29", "#6b6b66", "#e6e6e3", "#f0f0ee"

plt.rcParams.update({
    "font.size": 10, "axes.titlesize": 12, "axes.labelsize": 10, "axes.labelcolor": INK,
    "xtick.color": MUTED, "ytick.color": MUTED, "axes.edgecolor": "#b5b5b0", "figure.dpi": 110,
})


def _frame(title, subtitle, figsize=(9, 4.8)):
    """Başlık + kısa çıkarım (alt başlık, gerekirse iki satır) + sade eksenler."""
    import textwrap
    fig, ax = plt.subplots(figsize=figsize)
    lines = textwrap.wrap(subtitle, width=int(figsize[0] * 12.5))
    fig.text(0.075, 0.965, title, fontsize=13, fontweight="bold", color=INK, va="top")
    fig.text(0.075, 0.905, "\n".join(lines), fontsize=10, color=MUTED, va="top", linespacing=1.35)
    ax.grid(color=GRID, lw=0.8)
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)
    fig.subplots_adjust(left=0.075, right=0.86, top=0.82 - 0.045 * (len(lines) - 1), bottom=0.14)
    return fig, ax


def _end_label(ax, x, y, text, color):
    ax.annotate(text, (x, y), xytext=(8, 0), textcoords="offset points", va="center", fontsize=9,
                color=color, fontweight="bold")


def _mark_optimum(ax, x, y, text, dy=22):
    ax.plot(x, y, marker="o", ms=11, color=INK, mec="white", mew=2, zorder=6)
    ax.annotate(text, (x, y), xytext=(0, dy), textcoords="offset points", ha="center", fontsize=9,
                color=INK, fontweight="bold",
                bbox=dict(boxstyle="round,pad=0.35", fc="white", ec=GRID, lw=1))


# ---------------------------------------------------------------------------------------------
# 1) TAC – tüp sayısı
# ---------------------------------------------------------------------------------------------
def tac_vs_tubes(sweep, opt_label, Do, passes, T_water, selected=None):
    """sweep: esanjor_core.tac_vs_tubes_sweep çıktısı. Uygun bölge düz çizgi, sınır dışı kesikli."""
    ok = sweep[sweep.uygun]
    best = ok.loc[ok.TAC.idxmin()] if not ok.empty else None
    sub = ((f"En ucuz uygun nokta: {int(best.N)} tüp, L = {best.L / 1000:.1f} m, ${best.TAC / 1000:.1f} bin/yıl. "
            if best is not None else "")
           + f"Ø{Do} mm, {passes} geçiş, su çıkışı {T_water:.0f} °C sabit. Kesikli kısım kriterleri sağlamaz "
             "(solda alan yetersiz; sağda tüp hızı düşük; kabuk ΔP sınırı aşılan noktalar).")
    fig, ax = _frame("Yıllık toplam maliyet – tüp sayısı", sub)
    cols = {L: c for L, c in zip(sorted(sweep.L.unique()), (BLUE, GREEN, ORANGE, MUTED))}
    for L, d in sweep.groupby("L"):
        d = d.sort_values("N")
        col = cols[L]
        ax.plot(d.N, d.TAC / 1000, color=col, lw=1.6, ls=(0, (4, 3)), alpha=0.55)
        okd = d[d.uygun]
        if not okd.empty:
            ax.plot(okd.N, okd.TAC / 1000, color=col, lw=2.6, marker="o", ms=4.5, mec="white", mew=1)
        last = d.iloc[-1]
        _end_label(ax, last.N, last.TAC / 1000, f"L = {L / 1000:.1f} m", col)
    if selected is not None:                                   # (N, TAC $) seçilen tasarım
        ax.plot(selected[0], selected[1] / 1000, marker="D", ms=10, color=INK, mec="white", mew=1.5, zorder=7)
        _mark_optimum(ax, selected[0], selected[1] / 1000, f"{opt_label}: {int(selected[0])} tüp", dy=-36)
    elif best is not None:
        _mark_optimum(ax, best.N, best.TAC / 1000, f"{opt_label}: {int(best.N)} tüp", dy=-36)
    ax.set_xlabel("Tüp sayısı, N")
    ax.set_ylabel("Yıllık toplam maliyet [bin $/yıl]")
    return fig


# ---------------------------------------------------------------------------------------------
# 2) Sıcaklık – ısı yükü profili
# ---------------------------------------------------------------------------------------------
def temperature_profile(zones, T_dew, label="seçilen tasarım"):
    """En uygun tasarımda sıcak akış, soğutma suyu ve arayüz sıcaklığı (bölge sınırlarında)."""
    z = zones.copy()
    z["Q_b"] = z.dQ_kW.cumsum()
    z["Q_a"] = z.Q_b - z.dQ_kW
    Q = np.r_[z.Q_a.values, z.Q_b.values[-1]]
    Th = np.r_[z.T_a.values, z.T_b.values[-1]]
    tw = np.r_[z.t_su_a.values, z.t_su_b.values[-1]]
    dT = Th - tw
    i_pinch = int(np.argmin(dT))
    fig, ax = _frame(f"Sıcaklık – ısı yükü profili ({label})",
                     f"Yoğuşma {T_dew:.1f} °C'de başlar; en dar sıcaklık farkı {dT[i_pinch]:.1f} °C "
                     f"(Q = {Q[i_pinch]:.0f} kW). Arayüz, H₂ biriktikçe gazdan uzaklaşır.")
    ax.fill_between(Q, tw, Th, color=SHADE, lw=0, zorder=0)
    ax.plot(Q, Th, color=BLUE, lw=2.4)
    ax.plot(Q, tw, color=GREEN, lw=2.4)
    ok = z.T_i.notna()
    Qm = (0.5 * (z.Q_a + z.Q_b))[ok]
    ax.plot(Qm, z.T_i[ok], color=ORANGE, lw=2, ls=(0, (5, 3)))
    _end_label(ax, Q[-1], Th[-1], "Proses (kabuk)", BLUE)
    _end_label(ax, Q[-1], tw[-1], "Soğutma suyu (tüp)", GREEN)
    j = len(Qm) // 2
    ax.annotate("Arayüz sıcaklığı T_i", (Qm.iloc[j], z.T_i[ok].iloc[j]), xytext=(-10, -26),
                textcoords="offset points", ha="right", fontsize=9, color=ORANGE, fontweight="bold")
    ax.axhline(T_dew, color=MUTED, lw=1, ls=":")
    ax.text(Q[-1] * 0.55, T_dew + 1.2, f"çiy noktası {T_dew:.1f} °C", fontsize=8.5, color=MUTED)
    ax.annotate("", (Q[i_pinch], Th[i_pinch]), xytext=(Q[i_pinch], tw[i_pinch]),
                arrowprops=dict(arrowstyle="<->", color=INK, lw=1.2))
    ax.text(Q[i_pinch] - Q[-1] * 0.015, 0.5 * (Th[i_pinch] + tw[i_pinch]), f"ΔT min\n{dT[i_pinch]:.1f} °C",
            ha="right", va="center", fontsize=8.5, color=INK)
    ax.set_xlabel("Çekilen ısı, Q [kW]")
    ax.set_ylabel("Sıcaklık [°C]")
    ax.set_xlim(0, Q[-1] * 1.02)
    return fig


# ---------------------------------------------------------------------------------------------
# 3) Yıllık toplam maliyet – soğutma suyu çıkış sıcaklığı
# ---------------------------------------------------------------------------------------------
def tac_vs_water_temp(best, T_tested, selected=None):
    x, y = best.T_su_cikis.values, best.yillik_toplam_USD.values / 1000
    i = int(np.argmin(y))
    fig, ax = _frame("Yıllık toplam maliyet – soğutma suyu çıkış sıcaklığı",
                     f"En düşük maliyet su çıkışı {x[i]:.0f} °C'de: ${y[i]:.1f} bin/yıl "
                     f"({x[0]:.0f} °C'ye göre %{(1 - y[i] / y[0]) * 100:.0f} daha ucuz).")
    ax.plot(x, y, color=INK, lw=2.4, marker="o", ms=5, mec="white", mew=1)
    _mark_optimum(ax, x[i], y[i], f"En ucuz: {x[i]:.0f} °C · ${y[i]:.1f} bin")
    if selected is not None:                                   # (T_su_çıkış, TAC $) seçilen tasarım
        ax.plot(selected[0], selected[1] / 1000, marker="D", ms=10, color=INK, mec="white", mew=1.5, zorder=7)
        ax.annotate(f"Seçilen: {selected[0]:.0f} °C · ${selected[1] / 1000:.1f} bin", (selected[0], selected[1] / 1000),
                    xytext=(10, 14), textcoords="offset points", fontsize=8.8, fontweight="bold", color=INK,
                    bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="none"))
    _infeasible_band(ax, x, T_tested)
    ax.set_xlabel("Soğutma suyu çıkış sıcaklığı [°C]")
    ax.set_ylabel("Yıllık toplam maliyet [bin $/yıl]")
    ax.set_ylim(0, y.max() * 1.12)
    return fig


def _infeasible_band(ax, x_ok, x_tested):
    bad = [T for T in x_tested if T > max(x_ok)]
    if bad:
        ax.axvspan(min(bad) - 0.5, max(bad) + 0.5, color=SHADE, lw=0, zorder=0)
        ax.text(min(bad) - 0.3, 0.97, "uygun tasarım yok\n(F < 0,75)", transform=ax.get_xaxis_transform(),
                fontsize=8.5, color=MUTED, va="top")


# ---------------------------------------------------------------------------------------------
# 4) Yatırım – işletme dengesi
# ---------------------------------------------------------------------------------------------
def capex_opex_balance(best, T_tested):
    x = best.T_su_cikis.values
    cap = best.yillik_sermaye_USD.values / 1000
    ope = (best.elektrik_USD + best.su_USD).values / 1000
    i_opt = int(np.argmin(cap + ope))
    fig, ax = _frame("Yatırım – işletme dengesi",
                     "Su ısındıkça işletme maliyeti düşer, alan ve yatırım artar. Kesişim denge noktasıdır; "
                     f"en düşük toplam ise {x[i_opt]:.0f} °C'de (ikisi aynı olmak zorunda değil).")
    ax.plot(x, cap, color=BLUE, lw=2.4, marker="o", ms=5, mec="white", mew=1)
    ax.plot(x, ope, color=ORANGE, lw=2.4, marker="o", ms=5, mec="white", mew=1)
    _end_label(ax, x[-1], cap[-1], "Yatırım\n(yıllık sermaye payı)", BLUE)
    _end_label(ax, x[-1], ope[-1], "İşletme\n(elektrik + su)", ORANGE)
    d = ope - cap
    k = np.where(np.diff(np.sign(d)) != 0)[0]
    if len(k):
        k = k[0]
        f = d[k] / (d[k] - d[k + 1])
        xc, yc = x[k] + f * (x[k + 1] - x[k]), cap[k] + f * (cap[k + 1] - cap[k])
        _mark_optimum(ax, xc, yc, f"Denge: {xc:.1f} °C · ${yc:.1f} bin", dy=24)
    _infeasible_band(ax, x, T_tested)
    ax.set_xlabel("Soğutma suyu çıkış sıcaklığı [°C]")
    ax.set_ylabel("Yıllık maliyet [bin $/yıl]")
    ax.set_ylim(0, max(cap.max(), ope.max()) * 1.15)
    return fig


# ---------------------------------------------------------------------------------------------
# 5) En uygun modellerin karşılaştırması (kümülatif maliyet, başabaş)
# ---------------------------------------------------------------------------------------------
