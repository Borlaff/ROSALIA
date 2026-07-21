mkdir dist
rm dist/*
python -m build
python -m twine upload --repository pypi dist/* --verbose
