#!/usr/bin/env bash
pipenv run python setup.py sdist bdist_wheel

repo=testpypi
if [[ "$1" -eq "prod" ]]; then
	repo=pypi
fi
[ -f dist ] && rm -rf dist
[ -f build ] && rm -rf build
pipenv run twine upload --repository $repo dist/*
