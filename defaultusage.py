#! /usr/bin/python3

# Default usage of the "agent natomizer".

from anatomizer import AgentAnatomy


protein = AgentAnatomy('P00533')

features = protein.getfeatures()

print(features)

