cd ./docs 
sphinx-build -b html source build/html
sphinx-build -b html source build/html
cd ..
python3 -m build
python3 -m twine upload --repository pypi dist/*