#!/usr/bin/env bash
[ -d dist ] && rm -rf dist
[ -d build ] && rm -rf build
[ -d pamogk-egg.info ] && rm -rf pamogk-egg.info

pipenv run python setup.py sdist bdist_wheel

repo=testpypi
if [[ "$1" -eq "prod" ]]; then
	repo=pypi
fi

pipenv run twine upload --repository $repo dist/*
