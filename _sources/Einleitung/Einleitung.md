---
jupytext:
  cell_metadata_filter: -all
  formats: md:myst
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.16.7
kernelspec:
  display_name: Python 3 (ipykernel)
  language: python
  name: python3
---

# Einleitung
## Thema des Projekts
Die induktive Erwärmung ist ein faszinierendes Phänomen an der Schnittstelle von Elektromagnetismus und Thermodynamik, das in zahlreichen industriellen Anwendungen eine zentrale Rolle spielt. Von der Materialbearbeitung, wie dem Härten von Metallen oder dem Schmelzen, bis hin zum alltäglichen Induktionskochfeld – die Fähigkeit, Materialien berührungslos, effizient und gezielt zu erwärmen, bietet erhebliche Vorteile. Das Verständnis und die präzise Modellierung der zugrundeliegenden physikalischen Prozesse sind entscheidend für die Optimierung und Entwicklung solcher Systeme.

Das vorliegende Projekt widmet sich der numerischen Untersuchung der induktiven Erwärmung. Ziel ist es, ein gekoppeltes Problem zu modellieren, bei dem zunächst das elektromagnetische Feld mittels der Maxwell-Gleichungen beschrieben und die daraus resultierenden Wirbelströme in leitfähigen Materialien berechnet werden. Diese induzierten Ströme führen aufgrund des ohmschen Widerstands zur Dissipation von Energie in Form von Wärme. Anschliessend wird diese durch die Wirbelströme erzeugte Wärmequelle ermittelt und die resultierende stationäre Temperaturverteilung im Material mithilfe der Wärmeleitungsgleichung simuliert.

Im Rahmen dieser Arbeit werden die relevanten Maxwell-Gleichungen für das Wirbelstromproblem sowie die Wärmeleitungsgleichung betrachtet. Ein besonderer Fokus liegt auf der Herleitung der schwachen Formulierungen dieser Gleichungen, welche die Grundlage für die Anwendung der Methode der finiten Elemente (FEM) bilden. Die Untersuchung erfolgt exemplarisch an rotationssymmetrischen Geometrien, was eine Vereinfachung der Problemstellung ermöglicht, ohne die wesentlichen physikalischen Aspekte zu vernachlässigen. Diese Arbeit dokumentiert den Prozess der Modellierung, die numerische Implementierung und die Analyse der Ergebnisse dieses gekoppelten elektromagnetisch-thermischen Problems.

## Parameter und Konstanten
Die folgenden Parameter und Konstanten werden im gesamten Projekt verwendet.
Sie sind definiert in `global_parameters.py`. Die Werte basieren auf den Quellen {cite:p}`Stingelin2025HANAProjekt` und {cite:p}`GrossWassertechnikLeitfaehigkeit`. 

Die magnetische Permeabilität $\mu$ ist gegeben durch $\mu = \mu_r \cdot \mu_0$. Die relative Permeabilität ist ein Mass für die Feldverstärkung im Material (vgl. {numref}`mat_constants`). Vom Gebiet abhängige Konstanten können mit Hilfe von `CoefficientFunction()` definiert werden (vgl. {numref}`global_constants_include`).

```{table} Materialparameter
:name: mat_constants
| Material                     | Luft                               | Windungen                          | Kern / Topf                               | Wasser |
|------------------------------|:------------------------------------:|:------------------------------------:|:------------------------------------:|:------------------------------------:|
| el. Leitfähigkeit $\sigma$   | $0 \text{ S/m}$                  | $56 \cdot 10^6 \text{ S/m}$      | $56 \cdot 10^6 \text{ S/m}$      | $0.01 \text{ S/m}$ (Leitungswasser typ.: $0.005\, \ldots\, 0.05 \text{ S/m}$) |
| relative mag. Permeabilität $\mu_r$ | $1$                              | $1$                              | $1$                    | $1$ |
| Wärmeleitfähigkeit $\lambda$ | $0.0262 \text{ W/(m K)}$         | $400 \text{ W/(m K)}$            | $400 \text{ W/(m K)}$            | $0.6 \text{ W/(m K)}$            |
```

```{table} Weitere Parameter und Konstanten
:name: other_constants
| Parameter                             | Wert             |
|------------------------------|:-------------:|
| Frequenz $f$               | $50 \text{ Hz}$ |
| magnetische Feldkonstante $\mu_0$ | $4\pi \cdot 10^{-7}$ |
```

```{literalinclude} ../ImportExport/global_parameters.py
:language: python
:name: global_constants_include
:caption: Inhalt von global_parameters.py
```
