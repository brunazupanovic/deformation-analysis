import numpy as np
from scipy.stats import f,t
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import json


# ==========================================
# 1. KLASE I METODE
# ==========================================
class mjerenja():
    def __init__(self, p_from,p_to,vis_raz,duljina):
        self.p_to = p_to
        self.p_from= p_from
        self.vis_raz = vis_raz
        self.duljina = duljina
class TestHomogenosti:

    def __init__(self, o1, o2, f1, f2, alpha=0.05):
        self.o1 = o1
        self.o2 = o2
        self.f1 = f1
        self.f2 = f2
        self.alpha = alpha

    def test(self):
        F_ratio = (self.o1**2) / (self.o2**2)
        F_critical = f.ppf(1 - self.alpha, self.f1, self.f2)

        if F_ratio <= F_critical:
            print(
                f"Prihvaća se H0: Mjerenja su homogene točnosti (F={F_ratio:.3f}"
                f" <= F_crit={F_critical:.3f})"
            )
            return True
        else:
            print(
                f"Prihvaća se HA: Mjerenja nisu homogene točnosti (F={F_ratio:.3f}"
                f" > F_crit={F_critical:.3f})"
            )
            return False

    def varijanca(self):
        s0_2 = (self.f1 * (self.o1**2) + self.f2 * (self.o2**2)) / (
            self.f1 + self.f2
        )
        return np.sqrt(s0_2)

class tocka():
    def __init__(self,name,visina):
        self.name = name
        self._visina = visina

    @property
    def visina(self):
        return self._visina


class izjednacenje:
    def __init__(self,name,prib_visine,unknown):
        self.name = name
        self.prib_visine = prib_visine
        self.unknown = unknown
        self.list_mjerenja=[]
        self.results = {
            "x_izj": None,
            "v": None,
            "s0": None,
            "Qxx": None,
            "df": None
        }

    def dodaj_mjerenja(self,sirovi_podaci):
            for od_tocke, do_tocke, razlika, duljina in sirovi_podaci:
                novo_m = mjerenja(od_tocke, do_tocke, razlika, duljina)
                self.list_mjerenja.append(novo_m)


    def matrica_a(self):
        idx={name: i for i, name in enumerate(self.unknown)}
        m=len(self.list_mjerenja)
        n=len(self.unknown)
        A = np.zeros((m, n), dtype=float)
        for r,obs in enumerate(self.list_mjerenja):
            A[r, idx[obs.p_from]] = -1.0  # od koga mjerimo ide -1
            A[r, idx[obs.p_to]] = +1.0  # prema kome mjrenom 1
        return A

    def matrica_P(self, shd_mm_km=1.50):
        sig_mm = []
        for obs in self.list_mjerenja:
            L_km = obs.duljina / 1000.00
            sigma_i_mm = shd_mm_km * np.sqrt(L_km)
            sig_mm.append(sigma_i_mm)

        sig_mm = np.array(sig_mm, dtype=float)
        P = np.diag(1.0 / (sig_mm ** 2))
        return P

    def izvrsi_izjednacenje(self):
        A = self.matrica_a()
        P = self.matrica_P()

        # Formiranje vektora l (mjerenja - približne razlike)
        l_list = []
        for obs in self.list_mjerenja:
            # obs je objekt mjerenja, pa koristimo njegove atribute
            # r = h_to - h_from
            r = self.prib_visine[obs.p_to] - self.prib_visine[obs.p_from]
            l_list.append(obs.vis_raz - r)

        l = np.array(l_list).reshape(-1, 1)

        # normalne jednadžbe
        N = A.T @ P @ A
        n_vec = A.T @ P @ l

        # Pseudo-inverz jer je mreža slobodna (datum nije fiksiran)
        N_plus = np.linalg.pinv(N)
        x_korekcije = N_plus @ n_vec
        v = A @ x_korekcije - l

        # Izjednačene visine
        # H_izj = H_približno + x_korekcija
        H0 = np.array([self.prib_visine[name] for name in self.unknown]).reshape(-1, 1)
        x_izj = H0 + x_korekcije


        # kontrole i stupnjevi slobode
        m = len(self.list_mjerenja)
        n = len(self.unknown)
        print("ovo je m i n",m,n)# broj točaka
        df = m - (n - 1)

        vPv = (v.T @ P @ v)
        lpv=(-l.T@P@v)
        APv=(A.T @ P @ v)


        s0 = np.sqrt(float(vPv) / df) if df > 0 else 0

        # rezultati
        self.results["x_izj"] = x_izj
        self.results["v"] = v
        self.results["s0"] = s0
        self.results["Qxx"] = N_plus
        self.results["df"] = df

        kontrol_1 = np.allclose(APv, 0, atol=1e-4)  # atol_apsolutna tolerancija
        kontrol_2 = np.round(vPv, 4) == np.round(lpv, 4)

        if kontrol_1 and kontrol_2:
            print("Izjednačenje uspješno provedeno, kontrole odgovaraju (tolerancija 4 decimale).")

        return self.results


class globalni_pomak:
    def __init__(self,izj_1,izj_2,qxx1,qxx2,f1,f2,s0,alpha=0.05):
        self.izj_1 = izj_1
        self.izj_2 = izj_2
        self.qxx1 = qxx1
        self.qxx2 = qxx2
        self.f1 = f1
        self.f2 = f2
        self.s0 = s0
        self.alpha = alpha
    def pomak(self):
        d=self.izj_1-self.izj_2
        qd=self.qxx1+self.qxx2
        qd_plus=np.linalg.pinv(qd)
        h=np.linalg.matrix_rank(qd)
        s0_sq = self.s0 ** 2
        theta_s=(d.T@qd_plus@d)[0,0]/h
        f_def=theta_s/s0_sq
        f_crit_def=f.ppf(1-self.alpha,h,self.f1+self.f2)
        if f_def < f_crit_def:
            print("Prihvaća se H0: Nema značajnih deformacija.")
        else:
            print("Odbacuje se H0: Postoje deformacije.")
        return d, f_def, f_crit_def

    def lokalizacija_hannover(self, imena_tocaka=None):
        #Lokalizacija nestabilnih točaka,

        d = (self.izj_2 - self.izj_1).flatten()
        Qdd = self.qxx1 + self.qxx2
        Pdd = np.linalg.pinv(Qdd)  # Matrica težina pomaka

        num_points = len(d)
        if imena_tocaka is None:
            imena_tocaka = [f"T{i + 1}" for i in range(num_points)]

        aktivne_tocke = list(range(num_points))
        nestabilne_tocke = []
        povijest_testiranja = []

        f_potpuno = self.f1 + self.f2

        while len(aktivne_tocke) > 1:
            sigma_sq_j = {}


            for idx in aktivne_tocke:
                # Za 1D visinsku mrežu h_j = 1 stupanj slobode po točki
                d_j = d[idx]
                p_jj = Pdd[idx, idx]
                sigma_j = (d_j ** 2) * p_jj
                sigma_sq_j[idx] = sigma_j

            #  Pronalaženje točke s maksimalnom vrijednošću
            najnestabilnija_idx = max(sigma_sq_j, key=sigma_sq_j.get)
            max_sigma_sq = sigma_sq_j[najnestabilnija_idx]

            # Formiranje F_test = sigma_test^2 / s0^2
            F_test = max_sigma_sq / (self.s0 ** 2)

            # Kritična vrijednost F(h_j, f1+f2)
            F_crit = f.ppf(1 - self.alpha, 1, f_potpuno)

            ime_najnestabilnije = imena_tocaka[najnestabilnija_idx]

            # 5. Odluka o stabilnosti / nastavku iteracije
            if F_test > F_crit:
                # Točka je nestabilna, izbacuje se i ide se u novu iteraciju
                nestabilne_tocke.append(ime_najnestabilnije)
                aktivne_tocke.remove(najnestabilnija_idx)

                povijest_testiranja.append(
                    {
                        "tocka": ime_najnestabilnije,
                        "F_test": F_test,
                        "F_crit": F_crit,
                        "status": "NESTABILNA (izbačena)",
                    }
                )
            else:
                # Sve preostale točke su stabilne, prekid iteracije
                for idx in aktivne_tocke:
                    povijest_testiranja.append(
                        {
                            "tocka": imena_tocaka[idx],
                            "F_test": sigma_sq_j[idx] / (self.s0 ** 2),
                            "F_crit": F_crit,
                            "status": "Stabilna",
                        }
                    )
                break

        return povijest_testiranja, nestabilne_tocke



    def crtaj_lokalizaciju(self, rezultati_lokalizacije, naslov=None):
        """Metoda za vizualizaciju Hannoverske lokalizacije pomaka."""
        tocke = [r["tocka"] for r in rezultati_lokalizacije]

        # Dohvaćanje statistike neovisno o nazivu ključa u rječniku
        t_stat = [
            r.get("F_test", r.get("T_stat", r.get("F_stat")))
            for r in rezultati_lokalizacije
        ]
        f_crit = rezultati_lokalizacije[0]["F_crit"]

        # Crvena za pomak, zelena za stabilne točke
        boje = ["#e74c3c" if t > f_crit else "#2ecc71" for t in t_stat]

        plt.figure(figsize=(10, 6))

        bars = plt.bar(
            tocke, t_stat, color=boje, edgecolor="black", linewidth=1.2, width=0.5
        )

        plt.axhline(
            y=f_crit,
            color="red",
            linestyle="--",
            linewidth=2,
            label=f"F-dist. granica ({f_crit:.2f})",
        )

        for bar in bars:
            yval = bar.get_height()
            plt.text(
                bar.get_x() + bar.get_width() / 2.0,
                yval + (max(t_stat) * 0.02),
                f"{yval:.2f}",
                ha="center",
                va="bottom",
                fontsize=11,
            )

        if naslov is None:
            naslov = "Lokalizacija pomaka - Hannover metoda (F-test)"

        plt.title(naslov, fontsize=14, fontweight="bold", pad=15)
        plt.xlabel("Točke", fontsize=12, labelpad=10)
        plt.ylabel("F-statistika", fontsize=12, labelpad=10)
        plt.grid(axis="y", linestyle=":", alpha=0.6)

        legend_elements = [
            Patch(facecolor="#e74c3c", edgecolor="black", label="Otkriven pomak"),
            Patch(facecolor="#2ecc71", edgecolor="black", label="Referentna baza"),
            plt.Line2D(
                [0],
                [0],
                color="red",
                lw=2,
                ls="--",
                label=f"F-dist. granica ({f_crit:.2f})",
            ),
        ]
        plt.legend(handles=legend_elements, loc="upper left", frameon=True)

        plt.ylim(0, max(max(t_stat), f_crit) * 1.15)
        plt.tight_layout()
        plt.show()
# ==========================================
# 2. GLAVNI PROGRAM (Izvršavanje)
# ==========================================
# Ulazni podaci iz epone 2
unknown_order = ["A2", "A3", "O1", "O2", "O3", "O4"]
pocetne_visine = {
    "A2": 100.00000,
    "A3": 99.95660,
    "O1": 100.32220,
    "O2": 100.30360,
    "O3": 100.26520,
    "O4": 100.33180,
}
popis_mjerenja_2 = [
    ("A2", "A3", -0.04330, 23.0),
    ("A2", "O1", 0.32280, 40.0),
    ("A2", "O2", 0.30340, 50.0),
    ("A3", "O1", 0.36610, 35.0),
    ("A3", "O2", 0.34620, 40.0),
    ("O1", "O2", -0.01930, 22.0),
    ("O2", "O3", -0.03970, 8.0),
    ("O3", "O4", 0.06590, 22.0),
    ("O4", "O1", -0.00750, 8.0)
]

#ulazni podaic epohe 1
popis_mjerenja = [
    ("A2", "A3", -0.04340, 23.0),
    ("A2", "O1", 0.32220, 40.0),
    ("A2", "O2", 0.30330, 50.0),
    ("A3", "O1", 0.36600, 35.0),
    ("A3", "O2", 0.34650, 40.0),
    ("O1", "O2", -0.01860, 22.0),
    ("O2", "O3", -0.03840, 8.0),
    ("O3", "O4", 0.06660, 22.0),
    ("O4", "O1", -0.00960, 8.0)
]
koordinate_2d = {
    "A2": (457700.000, 5074620.000),  # Procijenjena pozicija blizu A3
    "A3": (457719.901, 5074623.290),
    "O1": (457716.500, 5074657.500),
    "O2": (457736.300, 5074656.000),
    "O3": (457736.900, 5074664.400),
    "O4": (457718.000, 5074665.700),
}



epoha1=izjednacenje("izjednacenje epohe 1",pocetne_visine,unknown_order)
epoha1.dodaj_mjerenja(popis_mjerenja)
res_oop=epoha1.izvrsi_izjednacenje()
print(res_oop)

epoha2=izjednacenje("izjednacenje epohe 2",pocetne_visine,unknown_order)
epoha2.dodaj_mjerenja(popis_mjerenja_2)
res_oop2=epoha2.izvrsi_izjednacenje()
print(res_oop2)
# stvaramo objekt testa homogenosti s dobivenim rezultatima
homo_test = TestHomogenosti(
    o1=res_oop["s0"],
    o2=res_oop2["s0"],
    f1=res_oop["df"],
    f2=res_oop2["df"],
    alpha=0.05,
)

# pokrece se test
je_homogeno = homo_test.test()
s0_zajednicki = homo_test.varijanca()
print(f"Zajednički s0 iz obje epohe: {s0_zajednicki*1000:.4f} mm")

# 1. Pozivamo globalni test s podacima koje smo već izračunali
test_deformacija = globalni_pomak(
    izj_1=res_oop["x_izj"],  # Vektor visina 1. epohe
    izj_2=res_oop2["x_izj"],  # Vektor visina 2. epohe
    qxx1=res_oop["Qxx"],  # Matrica Qxx 1. epohe
    qxx2=res_oop2["Qxx"],  # Matrica Qxx 2. epohe
    f1=res_oop["df"],  # df1 = 4
    f2=res_oop2["df"],  # df2 = 4
    s0=s0_zajednicki,  # Zajednički s0 iz F-testa
    alpha=0.05,
)

# 2. Pokrećemo izračun
d, F_def, F_crit = test_deformacija.pomak()

# 3. Ispisujemo rezultate
print(f"Izračunata F vrijednost : {F_def:.3f}")
print(f"Kritična F vrijednost   : {F_crit:.3f}")

rezultati_hannover, izbacene = test_deformacija.lokalizacija_hannover(
    imena_tocaka=["A2", "A3", "O1", "O2", "O3", "O4"]
)

print("\n--- HANNOVERSKA LOKALIZACIJA (PANDA MODEL) ---")
for r in rezultati_hannover:
    print(
        f"Točka {r['tocka']}: F_test = {r['F_test']:.3f} | F_crit = {r['F_crit']:.3f} -> {r['status']}"
    )

print("\nKonačno identificirane nestabilne točke:", izbacene)

#grafički prikaz
test_deformacija.crtaj_lokalizaciju(rezultati_hannover)

#IZOZ ZA QGIS
def izvoz_geojson_za_qgis(
    filename, unknown_order, koordinate_2d, res_oop, res_oop2, rezultati_hannover
):
    stat_dict = {r["tocka"]: r for r in rezultati_hannover}

    features = []
    for i, name in enumerate(unknown_order):
        east, north = koordinate_2d.get(name, (0.0, 0.0))
        h1 = res_oop["x_izj"][i, 0]
        h2 = res_oop2["x_izj"][i, 0]

        # Pomak u milimetrima i apsolutna vrijednost za veličinu simbola/kružnice
        pomak_mm = (h2 - h1) * 1000.0
        abs_pomak_mm = abs(pomak_mm)

        # Određivanje smjera vertikalnog pomaka
        if abs_pomak_mm < 0.01:
            smjer = "Bez pomaka"
        elif pomak_mm > 0:
            smjer = "Dizanje (+)"
        else:
            smjer = "Slijeganje (-)"

        stat_podaci = stat_dict.get(name, {})
        f_test = stat_podaci.get("F_test", 0.0)
        status = stat_podaci.get("status", "Nepoznato")

        # Izrada GeoJSON objekta za svaku točku
        feature = {
            "type": "Feature",
            "geometry": {"type": "Point", "coordinates": [east, north]},
            "properties": {
                "Tocka": name,
                "H_Epoha1_m": round(float(h1), 5),
                "H_Epoha2_m": round(float(h2), 5),
                "Pomak_mm": round(float(pomak_mm), 2),
                "Abs_Pomak_mm": round(float(abs_pomak_mm), 2),
                "Smjer": smjer,
                "F_test": round(float(f_test), 2),
                "Status": status,
            },
        }
        features.append(feature)

    geojson_data = {
        "type": "FeatureCollection",
        "crs": {
            "type": "name",
            "properties": {"name": "urn:ogc:def:crs:EPSG::3765"},
        },
        "features": features,
    }

    with open(filename, "w", encoding="utf-8") as f:
        json.dump(geojson_data, f, ensure_ascii=False, indent=4)

    print(f"GeoJSON datoteka '{filename}' uspješno je kreirana za QGIS!")

izvoz_geojson_za_qgis(
    "deformacije_mreze.geojson",
    unknown_order,
    koordinate_2d,
    res_oop,
    res_oop2,
    rezultati_hannover,
)
def izvoz_txt_izvjestaj(
    filename, unknown_order, res_oop, res_oop2, rezultati_hannover
):
    with open(filename, "w", encoding="utf-8") as f:
        f.write("=" * 65 + "\n")
        f.write("   IZVJEŠTAJ ANALIZE DEFORMACIJA NIVELMANSKE MREŽE (1D)\n")
        f.write("=" * 65 + "\n\n")

        # Hannover metoda - testiranje pojedinačnih točaka
        f.write("--- LOKALIZACIJA DEFORMACIJA (HANNOVER METODA) ---\n")
        f.write(
            f"{'Točka':<8} | {'H_Ep1 [m]':<10} | {'H_Ep2 [m]':<10} | {'d [mm]':<8} | {'F-test':<8} | {'Status'}\n"
        )
        f.write("-" * 65 + "\n")

        stat_dict = {r["tocka"]: r for r in rezultati_hannover}

        for i, name in enumerate(unknown_order):
            h1 = res_oop["x_izj"][i, 0]
            h2 = res_oop2["x_izj"][i, 0]
            pomak_mm = (h2 - h1) * 1000.0

            stat = stat_dict.get(name, {})
            f_test = stat.get("F_test", 0.0)
            status = stat.get("status", "Nije testirano")

            f.write(
                f"{name:<8} | {h1:<10.4f} | {h2:<10.4f} | {pomak_mm:<8.2f} | {f_test:<8.2f} | {status}\n"
            )

        f.write("=" * 65 + "\n")

    print(f"TXT izvještaj je uspješno spremljen pod: '{filename}'")


izvoz_txt_izvjestaj(
    "izvjestaj_deformacije.txt",
    unknown_order,
    res_oop,
    res_oop2,
    rezultati_hannover,
)


