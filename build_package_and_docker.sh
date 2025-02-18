#!/bin/bash

## Get version from Python project spec

build_package=true
build_docker=true
repository="pypi"
repourl="https://pypi.org/simple"

while test $# -gt 0; do
	case $1 in 
		--nopackage) 
			build_package=false
			shift
			;;
		--nodocker)
			build_docker=false
			shift
			;;
		--test)
			repository="testpypi"
			repourl="https://test.pypi.org/simple"
			shift
			;;
		*)
			echo "$0 Usage: "
			echo " "
			echo "--nopackage		Skip Python package build"
			echo "--nodocker		Skip Docker image build"
			echo "--test			Upload to testpypi"
			exit 0
			;;
	esac
done		

version="$(grep version pyproject.toml | cut -d' ' -f3)"
version="$(echo $version | tr -d '"')"

## Build the python package

if $build_package; then
	python -m build
	python3 -m twine upload -r $repository dist/*${version}*
	
fi

## Build Docker image

if $build_docker; then
	## CACHEBUST build arg doesn't do anything except force the pip install of sc_ssGSEA to
	## execute every time this image is built while using the cache for the rest of the layers.
	## See https://stackoverflow.com/questions/35134713/disable-cache-for-specific-run-commands

	## Dual build Mac (linux/arm64) Linux (linux/amd64)
	## For full debug log, add --progress=plain
	docker build -f Dockerfile --platform linux/amd64,linux/arm64 -t atwenzel/sc_ssgsea:v$version  \
		--build-arg CACHEBUST=$(date +%s) .
fi