#!/bin/bash
docker build -t "$CI_REGISTRY_IMAGE":base -f Dockerfile.base .
docker push "$CI_REGISTRY_IMAGE":base
