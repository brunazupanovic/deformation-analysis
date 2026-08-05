# Analiza deformacija nivelmanske mreže — Hannover metoda

Ovaj repozitorij sadrži računalne programe pisane u programskom jeziku Python za izjednačenje 1D nivelmanske mreže po metodi najmanjih kvadrata te detekciju i lokalizaciju nestabilnih točaka u sklopu deformacijske analize (PANDA/DEFANA model).

## Vizualizacija rezultata

![Prikaz QGIS](Prikaz_QGIS.png)

## Sadržaj repozitorija

Repozitorij obuhvaća slijedeće komponente:

1. **`deformacijska_analiza(OOP_20.07).py` (Glavna skripta):**
   - Objektno-orijentirana arhitektura koda s klasama za unos mjerenja, izjednačenje metodom najmanjih kvadrata i statističko testiranje.
   - Primjena Hannover metode za identifikaciju i eliminaciju nestabilnih točaka iz referentne osnove.
   - Automatski izvoz izračunatih pomaka u `.geojson` format za vizualizaciju u GIS okruženju (QGIS) te `.txt` izvještaj s numeričkim pokazateljima.
  
2. **`deformation_analysis_Hannover_IWST.py` (Prethodna verzija koda):**
   - Verzija koda koja provodi izjednačenje 1D mreže metodom najmanjih kvadrata i statističko testiranje
   - Primjena Hannover i IWST metode za identifikaciju nestabilnih točaka referetne osnove
   - Hannover metoda izbacuje nestabilne točke dok IWST metoda nestabilnim točkama smanjuje težinu i time ne utječe na geometriju mreže.

4. **Popratne datoteke:**
   - `panda_report.pdf` — Izvorni izvještaji izjednačenja i deformacijske analize iz programa PANDA/DEFANA.
   - `deformacije_mreze.geojson` — Generirani vektorski sloj s vektorima pomaka za QGIS.
   - `izvjestaj_deformacije.txt` — Tekstualni izvoz izračunatih parametara i statističkih testova.

## Metodologija i statistički testovi

- **Izjednačenje slobodne mreže:** Rješavanje rang-defekta mreže primjenom Moore-Penroseovog pseudo-inverza.
- **Statistička analiza:**
  - Test homogenosti varijanci iz dviju epoha mjerenja .
  - Globalni test deformacija mreže.
  - Hannoverska lokalizacija pomaka pojedinih točaka uz zadanu razinu značajnosti $\alpha = 0.05$.

## Rezultati analize

Analizom mjernih podataka iz dviju epoha identificirane su sljedeće nestabilne točke u mreži:

* **A3** — uvrštena u skupinu pomaknutih točaka (izdizanje)
* **O1** — uvrštena u skupinu pomaknutih točaka (izdizanje)
* **O3** — uvrštena u skupinu pomaknutih točaka (slijeganje)
* **O4** — uvrštena u skupinu pomaknutih točaka (slijeganje)

```text
=================================================================
   IZVJEŠTAJ ANALIZE DEFORMACIJA NIVELMANSKE MREŽE (1D)
=================================================================

--- LOKALIZACIJA DEFORMACIJA (HANNOVER METODA) ---
Točka    | H_Ep1 [m]  | H_Ep2 [m]  | d [mm]   | F-test   | Status
-----------------------------------------------------------------
A2       | 100.0001   | 100.0005   | 0.38     | 2.06     | Stabilna
A3       | 99.9566    | 99.9572    | 0.60     | 5.81     | NESTABILNA (izbačena)
O1       | 100.3223   | 100.3231   | 0.84     | 26.02    | NESTABILNA (izbačena)
O2       | 100.3035   | 100.3039   | 0.37     | 4.81     | Stabilna
O3       | 100.2651   | 100.2643   | -0.85    | 20.35    | NESTABILNA (izbačena)
O4       | 100.3318   | 100.3305   | -1.34    | 50.03    | NESTABILNA (izbačena)
=================================================================

```

Dobiveni rezultati podudaraju se s rezultatima softvera PANDA/DEFANA u izjednačenju iznosa pomaka i globalnom pomaku, dok do neslaganja dolazi kod lokalizacije točaka, budući da PANDA koristi t-test, a u ovom je radu primijenjen samo F-test
