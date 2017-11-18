#!/bin/bash

gfortran -shared -fPIC -g -o EDAllomUnitWrap.o EDAllomUnitWrap.F90

gfortran EDAllomUnitWrap.o -shared -fPIC -DALLOMUNITTEST -g -o FatesAllometryMod.o FatesAllometryMod.F90

