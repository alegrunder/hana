cd ..

call venv\Scripts\activate

cd venv_control

call pip freeze > requirements.txt

call deactivate