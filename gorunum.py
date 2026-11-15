"""Uygulamanın görünüm ayarları: başlık, sekmeler, görseller, renkler.

Ayarlar gorunum.json dosyasında tutulur ve uygulamadaki "Görünümü düzenle"
panelinden değiştirilir; kod açmaya gerek yoktur.
"""
import ast
import json
import os
import re

import pandas as pd
import streamlit as st

KLASOR = os.path.dirname(os.path.abspath(__file__))
AYAR_DOSYASI = os.path.join(KLASOR, "gorunum.json")
GORSEL_KLASORU = os.path.join(KLASOR, "gorseller")
TEMA_DOSYASI = os.path.join(KLASOR, ".streamlit", "config.toml")

# (kimlik, emoji, varsayılan ad) — kimlikler koddaki sekmelere bağlıdır, değiştirilmez
SEKMELER = [
    ("tanitim", "ℹ️", "Tanıtım"),
    ("otomatik", "🧭", "Otomatik tasarım"),
    ("ozet", "📋", "Özet"),
    ("tema", "📄", "TEMA föyü"),
    ("egri", "🌡️", "Yoğuşma eğrisi"),
    ("adimlar", "🧮", "Hesap adımları"),
    ("bolgeler", "🔁", "Bölgeler & doğrulama"),
    ("mekanik", "🔩", "Mekanik"),
    ("maliyet", "💰", "Maliyet"),
    ("optimizasyon", "🏆", "Tasarım taraması"),
    ("dwsim", "🧪", "DWSIM karşılaştırma"),
]

VARSAYILAN = {
    "pencere_adi": "Eşanjör Tasarımı",
    "baslik": "Kabuk–Tüp Eşanjör Tasarımı — Kısmi Yoğuşturucu",
    "alt_yazi": "",
    "girdiler_basligi": "Girdiler",
    "emoji": True,
    "girdiler_gorseli_genislik": 100,     # kenar çubuğu genişliğinin yüzdesi
    "baslik_gorseli_genislik": 300,       # piksel
    "sekmeler": [{"kimlik": k, "ad": ad, "emoji": e, "goster": True} for k, e, ad in SEKMELER],
    "metinler": {},     # {"m03": {"metin": "...", "goster": true}} — yalnızca değiştirilenler
    "tema": {"base": "light", "primaryColor": "#1f5f8b", "backgroundColor": "#ffffff",
             "secondaryBackgroundColor": "#eef2f6", "textColor": "#1b1f24", "font": "sans serif"},
}

_EMOJI = re.compile("[\U0001F300-\U0001FAFF☀-➿️]+\\s*")


def yukle():
    ayar = json.loads(json.dumps(VARSAYILAN))
    if os.path.exists(AYAR_DOSYASI):
        try:
            with open(AYAR_DOSYASI, encoding="utf-8") as f:
                kayitli = json.load(f)
            ayar.update({k: v for k, v in kayitli.items() if k in ayar})
        except (OSError, ValueError):
            pass
    # kayıtlı listede olmayan (yeni eklenmiş) sekmeleri sona ekle, bilinmeyenleri at
    bilinen = {k for k, _, _ in SEKMELER}
    liste = [s for s in ayar["sekmeler"] if s.get("kimlik") in bilinen]
    var = {s["kimlik"] for s in liste}
    eksik = [{"kimlik": k, "ad": ad, "emoji": e, "goster": True} for k, e, ad in SEKMELER if k not in var]
    liste = [s for s in eksik if s["kimlik"] == "tanitim"] + liste + [s for s in eksik if s["kimlik"] != "tanitim"]
    ayar["sekmeler"] = liste
    global _AKTIF
    _AKTIF = ayar
    return ayar


def kaydet(ayar):
    with open(AYAR_DOSYASI, "w", encoding="utf-8") as f:
        json.dump(ayar, f, ensure_ascii=False, indent=2)


_AKTIF = None
TUR_ADI = {"markdown": "Metin", "info": "Bilgi kutusu", "caption": "Küçük not", "warning": "Uyarı",
           "success": "Başarı mesajı", "subheader": "Alt başlık"}


class _Bos(dict):
    """Biçimlendirmede bilinmeyen {alan} olduğu gibi kalsın (yazıyı düzenleyen yanlış ad yazarsa hata vermesin)."""
    def __missing__(self, k):
        return "{" + k + "}"


def yaz(fonksiyon, anahtar, varsayilan, *args, degerler=None, **kwargs):
    """Uygulamadaki yazıları gösterir; panelde değiştirilmiş veya gizlenmişse ona uyar.
    degerler: yazıdaki {ad} yer tutucularına konacak sayılar (ör. {"Q": "3 045"}). Panelde yazı değiştirilirken
    bu yer tutucular korunursa sayılar yine otomatik gelir."""
    ayar = _AKTIF or VARSAYILAN
    ozel = ayar.get("metinler", {}).get(anahtar, {})
    if not ozel.get("goster", True):
        return None
    metin = ozel.get("metin", varsayilan)
    if not metin.strip():
        return None
    if degerler:
        metin = metin.format_map(_Bos(degerler))
    return fonksiyon(temizle(metin, ayar), *args, **kwargs)


@st.cache_data(show_spinner=False)
def metin_listesi(_degisim):
    """app.py'deki gv.yaz(...) çağrılarını bulur: anahtar, tür, varsayılan yazı ve bulunduğu yer."""
    yol = os.path.join(KLASOR, "app.py")
    agac = ast.parse(open(yol, encoding="utf-8").read())
    adlar = {k: f"{ad}" for k, _, ad in SEKMELER}
    liste = []

    def gez(dugum, yer):
        if isinstance(dugum, ast.For) and isinstance(dugum.iter, ast.Call) and \
                getattr(dugum.iter.func, "attr", "") == "sekme" and len(dugum.iter.args) == 2:
            yer = adlar.get(getattr(dugum.iter.args[1], "value", ""), yer)
        if isinstance(dugum, ast.FunctionDef) and dugum.name == "tanitim_goster":
            yer = adlar["tanitim"]
        if isinstance(dugum, ast.With) and any("sidebar" in ast.unparse(i.context_expr) for i in dugum.items):
            yer = "Sol menü"
        if isinstance(dugum, ast.Call) and getattr(dugum.func, "attr", "") == "yaz" and len(dugum.args) >= 3 \
                and isinstance(dugum.args[1], ast.Constant):
            f = ast.unparse(dugum.args[0])
            liste.append(dict(anahtar=dugum.args[1].value, varsayilan=dugum.args[2].value,
                              tur=TUR_ADI.get(f.rsplit(".", 1)[-1], "Metin"),
                              yer="Sol menü" if "sidebar" in f else yer, satir=dugum.lineno))
        for alt in ast.iter_child_nodes(dugum):
            gez(alt, yer)

    gez(agac, "Genel (sekmelerin dışı)")
    return sorted(liste, key=lambda x: x["satir"])


def temizle(metin, ayar):
    """Emoji kapalıysa metindeki emojileri siler."""
    return metin if ayar["emoji"] else _EMOJI.sub("", metin).strip()


def gorsel(ad):
    """gorseller/<ad>.png|jpg varsa yolunu döndürür."""
    for uzanti in ("png", "jpg", "jpeg"):
        yol = os.path.join(GORSEL_KLASORU, f"{ad}.{uzanti}")
        if os.path.exists(yol):
            return yol
    return None


def _gorsel_kaydet(ad, dosya):
    os.makedirs(GORSEL_KLASORU, exist_ok=True)
    _gorsel_sil(ad)
    uzanti = dosya.name.rsplit(".", 1)[-1].lower().replace("jpeg", "jpg")
    with open(os.path.join(GORSEL_KLASORU, f"{ad}.{uzanti}"), "wb") as f:
        f.write(dosya.getbuffer())


def _gorsel_sil(ad):
    while (yol := gorsel(ad)):
        os.remove(yol)


# ------------------------------------------------------------------ sekmeler
def sekmeleri_kur(ayar):
    """Görünen sekmeleri sırasıyla açar; {kimlik: sekme} sözlüğü döndürür."""
    gorunen = [s for s in ayar["sekmeler"] if s["goster"]] or ayar["sekmeler"][:1]
    adlar = [temizle(f"{s['emoji']} {s['ad']}".strip(), ayar) for s in gorunen]
    return dict(zip([s["kimlik"] for s in gorunen], st.tabs(adlar)))


def sekme(sekmeler, kimlik):
    """`for _ in sekme(T, "ozet"):` — sekme gizliyse içerik hiç çalışmaz."""
    if kimlik in sekmeler:
        with sekmeler[kimlik]:
            yol = gorsel(f"sekme_{kimlik}")
            if yol:
                st.image(yol, width="content")
            yield


# ------------------------------------------------------------------ düzenleme paneli
def _tema_yaz(tema):
    os.makedirs(os.path.dirname(TEMA_DOSYASI), exist_ok=True)
    satirlar = ["# Görünümü düzenle panelinden yazıldı.", "[theme]"]
    satirlar += [f'{k} = "{v}"' for k, v in tema.items()]
    satirlar += ["", "[client]", "toolbarMode = \"minimal\"", "", "[server]", "headless = true", "", "[browser]", "gatherUsageStats = false", ""]
    with open(TEMA_DOSYASI, "w", encoding="utf-8") as f:
        f.write("\n".join(satirlar))


def _durum_temizle():
    for k in [k for k in st.session_state if str(k).startswith(("mt_", "sekme_tablosu"))]:
        del st.session_state[k]


def duzenleme_paneli(ayar):
    st.subheader("Görünümü düzenle")
    st.caption("Değişiklikleri yapıp en alttaki **Kaydet**'e bas. Ayarlar gorunum.json dosyasına yazılır; "
               "hesaplara dokunulmaz.")

    yazi, aciklama, sek, gor, renk = st.tabs(["Başlıklar", "Sayfa yazıları", "Sekmeler", "Görseller", "Renkler"])

    with yazi:
        ayar["baslik"] = st.text_input("Ana başlık", ayar["baslik"])
        ayar["alt_yazi"] = st.text_area("Başlığın altındaki yazı (boş bırakılabilir)", ayar["alt_yazi"], height=80)
        ayar["girdiler_basligi"] = st.text_input("Sol menü başlığı", ayar["girdiler_basligi"])
        ayar["pencere_adi"] = st.text_input("Pencere / tarayıcı sekmesi adı", ayar["pencere_adi"])
        ayar["emoji"] = st.toggle("Sekme ve başlıklarda emoji göster", ayar["emoji"])

    with aciklama:
        metinler = metin_listesi(os.path.getmtime(os.path.join(KLASOR, "app.py")))
        yerler_ = list(dict.fromkeys(m["yer"] for m in metinler))
        yer_ = st.selectbox("Hangi bölümün yazıları?", yerler_,
                            help="Uygulamadaki açıklama, bilgi kutusu ve alt başlıklar bölüm bölüm listelenir.")
        st.caption("Bir yazıyı açıp değiştirebilir ya da **Göster** kutusunu kaldırarak gizleyebilirsin. "
                   "Kalın yazı için kelimeyi \\*\\*iki yıldız\\*\\* arasına al. **{süslü parantez}** içindeki adlar "
                   "hesaptan gelen sayılardır; yerini değiştirebilir ya da silebilirsin, adını değiştirme.")
        ayar.setdefault("metinler", {})
        for m in [m for m in metinler if m["yer"] == yer_]:
            k = m["anahtar"]
            ozel = ayar["metinler"].get(k, {})
            mevcut = ozel.get("metin", m["varsayilan"])
            goster = ozel.get("goster", True)
            ozet = " ".join(mevcut.replace("*", "").split())[:70] or "(boş)"
            isaret = "" if goster else " — gizli"
            with st.expander(f"{m['tur']}: {ozet}{isaret}"):
                g_ = st.checkbox("Göster", goster, key=f"mt_g_{k}")
                y_ = st.text_area("Yazı", mevcut, key=f"mt_y_{k}", height=68 if len(mevcut) < 120 else 160,
                                  label_visibility="collapsed")
                if g_ and y_ == m["varsayilan"]:
                    ayar["metinler"].pop(k, None)
                else:
                    ayar["metinler"][k] = {"metin": y_, "goster": g_}

    with sek:
        st.caption("**Göster** kutusunu kaldırınca sekme gizlenir. **Sıra** numarasıyla yerini değiştirebilir, "
                   "**Ad** ve **Emoji** hücrelerine tıklayıp yazabilirsin.")
        df = pd.DataFrame([{"Sıra": i + 1, "Göster": s["goster"], "Emoji": s["emoji"], "Ad": s["ad"],
                            "kimlik": s["kimlik"]} for i, s in enumerate(ayar["sekmeler"])])
        df = st.data_editor(
            df, hide_index=True, width="stretch", key="sekme_tablosu", height=36 * (len(df) + 1) + 3,
            column_config={"Sıra": st.column_config.NumberColumn(min_value=1, step=1, width="small"),
                           "Göster": st.column_config.CheckboxColumn(width="small"),
                           "Emoji": st.column_config.TextColumn(width="small"),
                           "kimlik": None},
            disabled=["kimlik"])
        df = df.sort_values("Sıra", kind="stable")
        ayar["sekmeler"] = [{"kimlik": r.kimlik, "ad": str(r.Ad or "").strip() or r.kimlik,
                             "emoji": str(r.Emoji or "").strip(), "goster": bool(r.Göster)}
                            for r in df.itertuples()]
        if not any(s["goster"] for s in ayar["sekmeler"]):
            st.warning("En az bir sekme açık olmalı; hepsi kapalıysa ilk sekme gösterilir.")

    with gor:
        yerler = {"girdiler": "Sol menüde, Girdiler başlığının üstü", "baslik": "Ana başlığın üstü (logo)"}
        yerler.update({f"sekme_{s['kimlik']}": f"'{s['ad']}' sekmesinin en üstü" for s in ayar["sekmeler"]})
        yer = st.selectbox("Görsel nereye konsun?", list(yerler), format_func=yerler.get)
        mevcut = gorsel(yer)
        c1, c2 = st.columns([2, 1])
        with c1:
            yuklenen = st.file_uploader("Resim seç (png / jpg)", type=["png", "jpg", "jpeg"], key=f"yukle_{yer}")
            if yuklenen is not None and st.button("Bu resmi buraya koy", type="primary"):
                _gorsel_kaydet(yer, yuklenen)
                st.rerun()
        with c2:
            if mevcut:
                st.image(mevcut, caption="Şu anki görsel", width="stretch")
                if st.button("Görseli kaldır"):
                    _gorsel_sil(yer)
                    st.rerun()
            else:
                st.caption("Burada henüz görsel yok.")
        if yer == "girdiler":
            ayar["girdiler_gorseli_genislik"] = st.slider("Genişlik (%)", 30, 100, ayar["girdiler_gorseli_genislik"])
        elif yer == "baslik":
            ayar["baslik_gorseli_genislik"] = st.slider("Genişlik (piksel)", 80, 1200, ayar["baslik_gorseli_genislik"])

    with renk:
        t = ayar["tema"]
        t["base"] = st.radio("Temel tema", ["light", "dark"], index=["light", "dark"].index(t["base"]),
                             format_func={"light": "Açık", "dark": "Koyu"}.get, horizontal=True)
        c1, c2 = st.columns(2)
        t["primaryColor"] = c1.color_picker("Vurgu rengi (butonlar, seçili sekme)", t["primaryColor"])
        t["textColor"] = c2.color_picker("Yazı rengi", t["textColor"])
        t["backgroundColor"] = c1.color_picker("Arka plan", t["backgroundColor"])
        t["secondaryBackgroundColor"] = c2.color_picker("Sol menü ve kutular", t["secondaryBackgroundColor"])
        fontlar = ["sans serif", "serif", "monospace"]
        t["font"] = st.selectbox("Yazı tipi", fontlar, index=fontlar.index(t["font"]))
        st.info("Renkler uygulamayı **kapatıp yeniden açınca** uygulanır.")

    c1, c2 = st.columns(2)
    if c1.button("Kaydet", type="primary", width="stretch"):
        kaydet(ayar)
        _tema_yaz(ayar["tema"])
        _durum_temizle()
        st.toast("Görünüm kaydedildi.")
        st.rerun()
    if c2.button("Varsayılana dön", width="stretch"):
        if os.path.exists(AYAR_DOSYASI):
            os.remove(AYAR_DOSYASI)
        _tema_yaz(VARSAYILAN["tema"])
        _durum_temizle()
        st.rerun()
