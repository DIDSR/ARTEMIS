#!/bin/bash

icc -c -fast -axSSSE3 -align -falign-functions -finline-functions Source/transEHP.c -o Source/transEHP.o

ifort -fast -axSSSE3 -align -falign-functions -finline-functions Source/penEasy_EMfield_EDEtally.f Source/tallyEnergyDepositionEvents.f Source/penelope.f Source/pengeom.f Source/penvared.f Source/penfield_uniform_z.f Source/transEHP.o -o ARTEMIS_v1.0.x

