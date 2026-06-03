#!/usr/bin/env bash

# Below is copied from 2x2 sim set up scripts (but using SPINE container bc of python dependencies). We may want to unify these at some point, but for now it's easier to just copy and paste and edit as needed.
export ARCUBE_RUNTIME=SHIFTER
export ARCUBE_CONTAINER=ghcr.io/deeplearnphysics/spine:0.13.1

source ../../2x2_sim/util/reload_in_container.inc.sh
source ../../2x2_sim/util/init.inc.sh