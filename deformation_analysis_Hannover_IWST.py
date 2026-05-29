import numpy as np
from scipy.stats import f, t
import matplotlib.pyplot as plt


# POMOCNE FUNKCIJE

def vektorl(mjerenja, prib_visine):
    """
    l_i =dh-H0,racuna vektor mjerenja l
    """
    l = []
    for (p_from, p_to, dh, _L) in mjerenja:
        r = prib_visine[p_to] - prib_visine[p_from]  # približna razlika
        l.append(dh - r)
    return np.array(l, dtype=float).reshape(-1, 1)


def colvec(x):
    return np.array(x, dtype=float).reshape(-1, 1)


def matrica_A(mjerenja, unknown_order):
    idx = {name: i for i, name in enumerate(
        unknown_order)}  # idx funkcija uzima A2 i vraca poziciju u ovom slucjau 0 i takoms e stvara matrica
    m = len(mjerenja)
    n = len(unknown_order)
    A = np.zeros((m, n), dtype=float)
    for r, (p_from, p_to, _dh, _L) in enumerate(mjerenja):
        A[r, idx[p_from]] = -1.0  # od koga mjrenimo ide -1
        A[r, idx[p_to]] = +1.0  # prema kome mjrenom 1
    return A
    print("ovo je marica a,", A)


def matrica_P(mjerenja, shd_mm_km=1.50):
    sig_mm = []
    for (_from, _to, _dh, L_m) in mjerenja:
        L_km = L_m / 1000.0
        sigma_i_mm = shd_mm_km * np.sqrt(L_km)  # u mm
        sig_mm.append(sigma_i_mm)

    sig_mm = np.array(sig_mm, dtype=float)
    P = np.diag(1.0 / (sig_mm ** 2))
    return P, sig_mm


unknown_order = ["A2", "A3", "O1", "O2", "O3", "O4"]


def izjednačenje(mjerenja, prib_visine, unknown_order):
    A = matrica_A(mjerenja, unknown_order);
    l = vektorl(mjerenja, prib_visine);
    P, sig_mm = matrica_P(mjerenja, shd_mm_km=1.5)  # ????? mislimd a je u pandi 1.5  tocnost pa se preuzeo isti broj
    # normlane jednazbe

    N = A.T @ P @ A
    n_vec = A.T @ P @ l

    N_plus = np.linalg.pinv(N)  # N je singularna
    x = N_plus @ n_vec  # korekcije visina
    v = A @ x - l

    # izjednacenje visina
    x0 = colvec([pocetne_visine[k] for k in unknown_order])
    x_izj = x0 + x

    m = A.shape[0]
    n = A.shape[1]
    df = m - (n - 1)

    s0 = float(np.sqrt((v.T @ P @ v) / df))

    Qxx = N_plus  # kofaktorska matrica

    # Kontrole
    check1 = A.T @ P @ v
    check2a = float(v.T @ P @ v)
    check2b = float(-l.T @ P @ v)

    return {
        "unknown_order": unknown_order,
        "N": N,
        "A": A,
        "P": P,
        "sigma_mm": sig_mm,
        "l": l,
        "x_izj": x_izj,
        "v": v,
        "s0": s0,
        "Qxx": Qxx,
        "df": df,
        "check_ATPv": check1,
        "check_vTPv": check2a,
        "check_minus_lTPv": check2b
    }


# -------------------------9
# EPOHA 1 - PODACI
# -------------------------


mjerenja = [
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

mjerenja_11 = colvec([row[2] for row in mjerenja])

pocetne_visine = {
    "A2": 100.00000,
    "A3": 99.95670,
    "O1": 100.32280,
    "O2": 100.30350,
    "O3": 100.26380,
    "O4": 100.33030
}

# -------------------------
# izjednacenje 1 epohe
# -------------------------
res1 = izjednačenje(mjerenja, pocetne_visine, unknown_order)

print("unknown_order =", res1["unknown_order"])
print("s0 =", res1["s0"], "df =", res1["df"])
print("Kontrola A^T v (treba biti ~0):", res1["check_ATPv"].flatten())
print("Kontrola v^T v i -l^T v:", res1["check_vTPv"], res1["check_minus_lTPv"])
print("x_izj =", res1["x_izj"].flatten())

# -------------------------
# EPOHA 2 - PODACI
# -------------------------

mjerenja_2 = [
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

mjerenja_22 = colvec([row[2] for row in mjerenja_2])

pocetne_visine_2 = {
    "A2": 100.00000,
    "A3": 99.95670,
    "O1": 100.32280,
    "O2": 100.30350,
    "O3": 100.26380,
    "O4": 100.33030
}

# -------------------------
# izjednacenje 1 epohe
# -------------------------
res2 = izjednačenje(mjerenja_2, pocetne_visine_2, unknown_order)

print("unknown_order =", res2["unknown_order"])
print("s0 =", res2["s0"], "df =", res2["df"])
print("Kontrola A^T v (treba biti ~0):", res2["check_ATPv"].flatten())
print("Kontrola v^T v i -l^T v:", res2["check_vTPv"], res2["check_minus_lTPv"])
print("x_izj =", res2["x_izj"].flatten())

# x0 = početne visine nepoznanica (bez A2) kao STUPČASTI vektor (5x1)
# redoslijed ispisa A: [A3, O1, O2, O3, O4]
unknown_order = ["A2", "A3", "O1", "O2", "O3", "O4"]

# -------------------------
# TEST HOMOGENOSTI VARIJANCI (F-test)
# -------------------------
so = res1["s0"]
so_2 = res2["s0"]
alpha = 0.05
df1 = res1["df"]
df2 = res2["df"]
F_ratio = (so ** 2) / (so_2 ** 2)
F_critical = f.ppf(1 - alpha, df1, df2)
o = (df1 * so ** 2 + df2 * so_2 ** 2) / (df1 + df2)

if F_ratio <= F_critical:
    print("Varijance su homogene, varijaca iznosi o =", o)
else:
    print("Varijance nisu homogene  =", o)

# IZEJDNACENJE MJRENJA
izj_mjerenja = mjerenja_11 + res1["v"]
izj_mjerenja_2 = mjerenja_22 + res2["v"]

razlika_L = izj_mjerenja - izj_mjerenja_2
print("Razlika izjednačenih mjerenja (9x1) shape:", razlika_L, razlika_L.shape)

# Pomaci točaka (OVO je d za deformacije / IWST)
x_izj_1 = res1["x_izj"]
x_izj_2 = res2["x_izj"]
d = x_izj_1- x_izj_2
print("d (x_izj_1 - x_izj_2) shape:", d, d.shape)
print("d =", d.flatten())

# Qd
Qd = res1["Qxx"] + res2["Qxx"]
print(Qd, "ovo je qd")

Qd_plus = np.linalg.pinv(Qd)
h = np.linalg.matrix_rank(Qd)

theta_s = (d.T @ Qd_plus @ d)[0, 0] / h
F_def = theta_s / o

Fcrit_def = f.ppf(1 - alpha, h, df1 + df2)

print("F_def =", F_def, "Fcrit_def =", Fcrit_def)

if F_def < Fcrit_def:
    print("Prihvaća se H0: Nema značajnih deformacija.")
else:
    print("Odbacuje se H0: Postoje deformacije.")



# =========================================================================
# LOKALIZACIJA -
# =========================================================================

# 1. Podjela točaka prema PANDA reportu
reference_points = ["A2", "A3", "O1", "O2"]
object_points = ["O3", "O4"]

# Izračun kritične vrijednosti (alfa=0.05, f1=1, f2=df1+df2)

F_crit_panda = f.ppf(0.95, 1, df1 + df2)

print("\n" + "=" * 90)
print(f"{' LOKALIZACIJA POMAKA ':^90}")
print("=" * 90)

# 2. Testiranje OBJEKTNIH točaka (onih koje pratimo)
print(f"{'Br.':<5} {'Točka':<10} {'dz [mm]':>10} {'T-stat':>12} {'F-crit':>12} {'Rezultat':<15}")
print("-" * 90)

for i, name in enumerate(unknown_order):
    if name in object_points:
        dz_mm = d[i, 0] * 1000.0
        # T-statistika: kvadrat pomaka / (varijanca * kofaktor razlike)
        T_val = (d[i, 0] ** 2) / (o * Qd[i, i])

        status = "!!! POMAK !!!" if T_val > F_crit_panda else "Stabilna"
        print(f"{i + 1:<5} {name:<10} {dz_mm:>10.2f} {T_val:>12.2f} {F_crit_panda:>12.2f} {status:<15}")

# 3. Prikaz REFERENTNIH točaka (onih koje čine datum)
print("-" * 90)
print(f"{'Tip':<15} {'Točka':<10} {'dz [mm]':>10} {'T-stat':>12} {'Status':<15}")
print("-" * 90)
for i, name in enumerate(unknown_order):
    if name in reference_points:
        dz_mm = d[i, 0] * 1000.0
        T_val = (d[i, 0] ** 2) / (o * Qd[i, i])

        print(f"{'Referenca':<15} {name:<10} {dz_mm:>10.2f} {T_val:>12.2f} {'Referetna točka'}")

print("=" * 90)


# pregledi st dostupanaj pomaka i 95% intervala
def pregled_intervali_ispis(res1, res2, o, df_total, conf=0.95):
    names = res1["unknown_order"]

    # visine (izjednačene)
    H1 = res1["x_izj"].flatten()
    H2 = res2["x_izj"].flatten()

    # pomak
    d = (H1 - H2).reshape(-1, 1)

    # dijagonale kofaktorske matrice
    q1 = np.diag(res1["Qxx"])
    q2 = np.diag(res2["Qxx"])

    # standardna odstupanja visina (epoha 1 i epoha 2)
    sH1 = np.sqrt(np.clip(o * q1, 0, np.inf)).reshape(-1, 1)
    sH2 = np.sqrt(np.clip(o * q2, 0, np.inf)).reshape(-1, 1)

    # standardno odstupanje pomaka
    sd = np.sqrt(np.clip(o * (q1 + q2), 0, np.inf)).reshape(-1, 1)

    # t-faktor za konf. interval
    alpha = 1 - conf
    k = t.ppf(1 - alpha / 2, df_total)

    half = k * sd

    print(f"\n--- PREGLED (1D 'elipse' = intervali) ---")
    for i, nm in enumerate(names):
        print(
            f"{nm:>2} | σH1={sH1[i, 0]:.6f} m | σH2={sH2[i, 0]:.6f} m | "
            f"d={d[i, 0]:+.6f} m | {int(conf * 100)}% ±{half[i, 0]:.6f} m"
        )


pregled_intervali_ispis(res1, res2, o, df_total=df1 + df2, conf=0.95)

# 1. Definiramo indekse stabilnih točaka prema PANDA izvještaju
# To su A2, A3, O1, O2
stabilni_idx = [0, 1, 2, 3]
n_tocka = len(unknown_order)
m_stab = len(stabilni_idx)

# 2. Kreiramo  S-matricu
S_han = np.eye(n_tocka)
for i in range(n_tocka):
    for j in stabilni_idx:
        S_han[i, j] -= 1.0 / m_stab

# 3. Transformiramo pomake (dz u Pandi)
d_han = S_han @ d
Qd_han = S_han @ Qd @ S_han.T

print(        "\n--- REZULTATI HANNOVER LOKALIZACIJE  ---")
print("\n" + "=" * 90)
print(f"{' LOKALIZACIJA POMAKA ':^90}")
print("=" * 90)

# 2. Testiranje OBJEKTNIH točaka (onih koje pratimo)
print(f"{'Br.':<5} {'Točka':<10} {'dz [mm]':>10} {'T-stat':>12} {'F-crit':>12} {'Rezultat':<15}")
print("-" * 90)
for i, name in enumerate(unknown_order):
    dz_mm = d[i, 0] * 1000.0
    Ti_han = (d_han[i, 0]**2) / (o * Qd_han[i, i])
    status = "!!! POMAK !!!" if Ti_han > F_crit_panda else "Stabilna"
    print(f"{i + 1:<5} {name:<10} {dz_mm:>10.2f} {Ti_han:>12.2f} {F_crit_panda:>12.2f} {status:<15}")


W = np.eye(len(unknown_order))
# grafovii
# =========================================================================
# GRAFIČKI PRIKAZ LOKALIZACIJE
# =========================================================================

def plot_panda_localisation(unknown_order, d, o, Qd, F_crit, object_points):
    names = unknown_order
    T_values = []
    colors = []

    # Izračun T-vrijednosti za sve točke radi grafa
    for i, name in enumerate(names):
        T = (d_han[i, 0]**2) / (o * Qd_han[i, i])
        T_values.append(T)


        # Ako je u object_points i T > F_crit -> Crvena (Pomak)
        # Ako je u object_points i T <= F_crit -> Plava (Stabilna objektna)
        # Ako je u reference_points -> Zelena (Referentna točka)
        if name in object_points:
            if T > F_crit:
                colors.append('#e74c3c')  # Crvena
            else:
                colors.append('#3498db')  # Plava
        else:
            colors.append('#2ecc71')  # Zelena

    plt.figure(figsize=(10, 6))
    bars = plt.bar(names, T_values, color=colors, edgecolor='black', alpha=0.8)

    # Crta kritične vrijednosti
    plt.axhline(y=F_crit, color='red', linestyle='--', linewidth=2, label=f'Kritična granica (F-crit = {F_crit:.2f})')

    # Estetika grafa
    plt.title('Lokalizacija pomaka - Hannover metoda (T-test)', fontsize=14, fontweight='bold')
    plt.ylabel('T-statistika', fontsize=12)
    plt.xlabel('Točke', fontsize=12)
    plt.grid(axis='y', linestyle=':', alpha=0.7)

    # Dodavanje vrijednosti iznad stupaca
    for bar in bars:
        yval = bar.get_height()
        plt.text(bar.get_x() + bar.get_width() / 2, yval + 1, f'{yval:.2f}', ha='center', va='bottom', fontsize=10)

    # Legenda
    from matplotlib.lines import Line2D
    custom_lines = [Line2D([0], [0], color='#e74c3c', lw=4),
                    Line2D([0], [0], color='#2ecc71', lw=4),
                    Line2D([0], [0], color='red', linestyle='--')]
    plt.legend(custom_lines, ['Otkriven pomak', 'Referentna baza', 'F-dist. granica (5.33)'])

    plt.tight_layout()
    plt.show()


# Poziv funkcije
plot_panda_localisation(unknown_order, d, o, Qd, F_crit_panda, object_points)


# IWST METODA - ITERATIVNI ISPIS SVIH KORAKA


print(f"\nPOKRETANJE IWST METODE (Iterative Weighted Similarity Transformation)")

d_iwst = d.copy()
n_points = len(unknown_order)  #
W = np.eye(n_points)  # Početne težine su 1
G = np.ones((n_points, 1))

max_iter = 10
tolerancija = 0.0001

for i in range(max_iter):
    d_old = d_iwst.copy()

    # 1. Izračun matrice S-transformacije (R)
    # R = I - G * inv(G^T * W * G) * G^T * W
    inv_part = np.linalg.inv(G.T @ W @ G)
    S = np.eye(n_points) - (G @ inv_part @ G.T @ W)

    # 2. Transformacija pomaka
    d_iwst = S @ d_iwst

    # ISPIS TRENUTNE ITERACIJE
    print(f"\nITERACIJA {i + 1}:")
    print("-" * 20)
    for idx, name in enumerate(unknown_order):
        print(f"Točka {name}: pomak = {d_iwst[idx, 0]:.6f} m, težina = {W[idx, idx]:.2f}")

    # 3. Ažuriranje težina za SLJEDEĆU iteraciju (w = 1 / |d|)

    weights = 1.0 / (np.abs(d_iwst.flatten()) + 1e-6)
    W = np.diag(weights)
    print(f"\nNove težine nakon ieteracije su")

    for idx, name in enumerate(unknown_order):
        print(f"Točka {name}: pomak={d_iwst[idx,0]:0.6f} m, težine= {W[idx, idx]:.4f}")


    # 4. Provjera konvergencije
    razlika = np.max(np.abs(d_iwst - d_old))
    if razlika < tolerancija:
        print(f"\nIWST je konvergirao u {i + 1}. iteraciji jer je promjena manja od {tolerancija} m.")
        break

# FINALNI TEST NAKON IWST-a
# Transformacija kofaktorske matrice Qd
Qd_iwst = S @ Qd @ S.T
Qd_plus_iwst = np.linalg.pinv(Qd_iwst)
h_iwst = np.linalg.matrix_rank(Qd_iwst)

theta_s_iwst = (d_iwst.T @ Qd_plus_iwst @ d_iwst) / h_iwst
F_def_iwst = float((theta_s_iwst / o).ravel()[0])

print(f"FINALNI REZULTAT NAKON IWST TRANSFORMACIJE:")
print(f"Novi F_def = {F_def_iwst:.4f}")
print(f"Kritični F_crit = {Fcrit_def:.4f}")

if F_def_iwst < Fcrit_def:
    print("ZAKLJUČAK: Mreža je stabilna (nema značajnih deformacija).")
else:
    print("ZAKLJUČAK: Otkrivene su značajne deformacije na objektu.")

# poslije IWST-a, zadnja iteracija
w_after = np.diag(W).astype(float)
print(w_after)

x = np.arange(len(unknown_order))
width = 0.35

plt.figure()
plt.bar(x + width / 2, w_after, width, label="Nakon IWST", alpha=0.6)
plt.xticks(x, unknown_order)
plt.xlabel("Točka")
plt.ylabel("Težina")
plt.title("Promjena težina točaka nakon IWST metode")
plt.grid(True, axis="y", which="both")
plt.legend()
plt.show()













