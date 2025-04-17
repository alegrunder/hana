# HANA Projekt
HANA Numerik SEP Projekt: Induktive Erwärmung

# Jupyter-Book 
[>> Link Webseite Jupyter-Book](https://jupyterbook.org/en/stable/start/overview.html)

## Installieren
1. In hana venv wechseln, etwa so: `$ source myvenv/hana/bin/activate` oder mit `start_venv.bat` script in venv_control
2. `$ pip install -U jupyter-book`

## Book als HTML builden
1. In hana venv wechseln, etwa so: `$ source myvenv/hana/bin/activate` oder mit `start_venv.bat` script in venv_control
2. `$ cd xy/Github/hana` in HANA-Github Folder
3. `$ jupyter-book build hana_project_book/`
4. ./hana_project_book/_build/html/index.html öffnen


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
