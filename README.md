# HANA Projekt
HANA Numerik SEP Projekt: Induktive Erwärmung

[>> Link zur Webpage (gh-pages)](https://alegrunder.github.io/hana/intro.html)

# Jupyter-Book Documentation
[>> Link Webseite Jupyter-Book](https://jupyterbook.org/en/stable/start/overview.html)

## Jupyter-Book installieren
1. In hana venv wechseln, etwa so:
    - Mac: `$ source myvenv/hana/bin/activate`
    - Windows: `start_venv.bat` script in venv_control
2. `$ pip install -U jupyter-book`

## Book als HTML builden
1. In hana venv wechseln, etwa so: `$ source myvenv/hana/bin/activate` oder mit `start_venv.bat` script in venv_control
2. `$ cd xy/Github/hana` in HANA-Github Folder
3. `$ jupyter-book build hana_project_book/`
4. ./hana_project_book/_build/html/index.html öffnen

## HTML auf gh-pages publishen
1. Sicherstellen dass in /hana/ folder
2. Book als HTML builden (siehe oben)
3. `$ ghp-import -n -p -f hana_project_book/_build/html`
    - `-n`: Creates a .nojekyll file (important for GitHub Pages to serve sites correctly, especially if they contain files/folders starting with underscores).
    - `-p`: Pushes the gh-pages branch to GitHub.
    - `-f`: Forces the push (overwrites the existing gh-pages branch).
4. Username wird abgefragt: `alegrunder`
5. Password: Token verwenden
6. Ein bisschen Geduld haben, Github macht Checks (Status unter "Branches" im Repo sichtbar)

## Quickly add YAML metadata for MyST Notebooks
If you have a markdown file and you’d like to quickly add YAML metadata to it, so that Jupyter Book will treat it as a MyST Markdown Notebook, run the following command:

`$ jupyter-book myst init hana_project_book/markdownfile.md`

> Keine Ahnung, ob das notwendig ist

## PDF Export
### Install
1. `$ pip install playwright`
2. `$ playwright install --with-deps chromium`

### Build
1. In ./hana_project_book/_static/myfile.css irgendwelche Dinge mit `@media print` ergänzen
2. `$ jupyter-book build mybookname/ --builder pdfhtml`
3. ./hana_project_book/_build/pdf/book.pdf öffnen
