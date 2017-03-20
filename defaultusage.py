#! /usr/bin/python3

# Default usage of the "agent anatomizer".

from anatomizer import AgentAnatomy


protein = AgentAnatomy('P19174')

features = protein.getfeatures()

print(features)

