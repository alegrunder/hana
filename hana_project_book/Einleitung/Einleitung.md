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

## Parameter und Konstanten
Die folgenden Parameter und Konstanten werden im gesamten Projekt verwendet.
Sie sind definiert in `global_parameters.py`.

Die magnetische Permeabilität $\mu$ ist gegeben durch $\mu = \mu_r \cdot \mu_0$. Die relative Permeabilität ist ein Mass für die Feldverstärkung im Material (vgl. {numref}`mat_constants`). Vom Gebiet abhängige Konstanten können mit Hilfe von `CoefficientFunction()` definiert werden (vgl. {numref}`global_constants_include`).

```{table} Materialparameter
:name: mat_constants
| Material                     | Luft                               | Windungen                          | Kern                               |
|------------------------------|:------------------------------------:|:------------------------------------:|:------------------------------------:|
| el. Leitfähigkeit $\sigma$   | $0 \text{ S/m}$                  | $56 \cdot 10^6 \text{ S/m}$      | $56 \cdot 10^6 \text{ S/m}$      |
| relative mag. Permeabilität $\mu_r$ | $1$                              | $1$                              | $1$                              |
| Wärmeleitfähigkeit $\lambda$ | $0.0262 \text{ W/(m K)}$         | $400 \text{ W/(m K)}$            | $400 \text{ W/(m K)}$            |
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
:caption: Inhalt von project_constants.py
```
