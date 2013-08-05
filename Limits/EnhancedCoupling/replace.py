import os, sys

file = open("test", "r")
lines = file.readlines()
file.close()

file = open("input2.txt", "w")
for line in lines:
    newline = line.replace("\t", " ")
    file.write(newline)

file.close()
