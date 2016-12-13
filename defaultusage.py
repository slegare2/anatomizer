#! /usr/bin/python3

# Default usage of the "agent anatomizer".

from anatomizer import AgentAnatomy


protein = AgentAnatomy('Q9NWT1')

features = protein.getfeatures()

print(features)

