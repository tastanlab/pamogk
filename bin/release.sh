#!/usr/bin/env bash
pipenv run python setup.py sdist bdist_wheel

repo=testpypi
if [[ "$1" -eq "prod" ]]; then
	repo=pypi
fi
pipenv run twine upload --repository $repo dist/*
