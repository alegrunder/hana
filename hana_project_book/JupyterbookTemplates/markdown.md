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

# Ein paar spannende Dinge
$$\begin{split}
-u''(x) & = f(x)\quad\forall\ x\in (0,1)\\
u(0) & = u(1) = 0
\end{split}$$(eq:eindimrwp)

Das zur Gleichung {eq}`eq:eindimrwp` analoge Randwertproblem ist in dem Fall für die Poisson Gleichung gegeben durch

```{code-cell} ipython3
from netgen.geom2d import unit_square
from ngsolve import *
from ngsolve.webgui import Draw
```

```{code-cell} ipython3
mesh = Mesh(unit_square.GenerateMesh(maxh=0.1))
Draw(mesh);
```

```{table} Weitere Parameter und Konstanten
:name: other_constants_copy
| Parameter                             | Wert             |
|------------------------------|-------------|
| Frequenz $f$               | $50 \text{ Hz}$ |
| magnetische Feldkonstante $\mu_0$ | $4\pi 10^{-7}$ |
```

Siehe {numref}`mat_constants`.

(ref:IntroPoissonSchwacheGleichung)=
## Schwache Gleichung
Hier ein Verweis auf dieses Kapitel: {ref}`ref:IntroPoissonSchwacheGleichung`.



# Markdown Files

Whether you write your book's content in Jupyter Notebooks (`.ipynb`) or
in regular markdown files (`.md`), you'll write in the same flavor of markdown
called **MyST Markdown**.
This is a simple file to help you get started and show off some syntax.

## What is MyST?

MyST stands for "Markedly Structured Text". It
is a slight variation on a flavor of markdown called "CommonMark" markdown,
with small syntax extensions to allow you to write **roles** and **directives**
in the Sphinx ecosystem.

For more about MyST, see [the MyST Markdown Overview](https://jupyterbook.org/content/myst.html).

## Sample Roles and Directives

Roles and directives are two of the most powerful tools in Jupyter Book. They
are like functions, but written in a markup language. They both
serve a similar purpose, but **roles are written in one line**, whereas
**directives span many lines**. They both accept different kinds of inputs,
and what they do with those inputs depends on the specific role or directive
that is being called.

Here is a "note" directive:

```{note}
Here is a note
```

It will be rendered in a special box when you build your book.

Here is an inline directive to refer to a document: {doc}`markdown-notebooks`.

+++

## Citations

You can also cite references that are stored in a `bibtex` file. For example,
the following syntax: `` {cite}`holdgraf_evidence_2014` `` will render like
this: *removed*

Moreover, you can insert a bibliography into your page with this syntax:
The `{bibliography}` directive must be used for all the `{cite}` roles to
render properly.
For example, if the references for your book are stored in `references.bib`,
then the bibliography is inserted with:


## Learn more

This is just a simple starter to get you started.
You can learn a lot more at [jupyterbook.org](https://jupyterbook.org).
