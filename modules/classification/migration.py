from rdkit import Chem
from index import Init


strInChI = "InChI=1S/H4N2.2NO3.Ni/c1-2;2*2-1(3)4;/h1-2H2;;;/q;2*-1;+2"
print(strInChI, "\n")
mol = Chem.inchi.MolFromInchi(strInChI)

# Code START
param = Init(mol)
param["iPartNum"] = strInChI.count("/")

# seperate each layer
temp = strInChI.split("/")[0]
param["strTFormula"] = strInChI.split("/")[1]

for i in range(2, param["iPartNum"] + 1):
    temp = strInChI.split("/")[i]
    c1 = temp[0]

    # ---main layer
    # atom layer
    c2 = "c"
    if c1 == c2:
        param["strTAtom"] += temp
    # hydrogen
    c2 = "h"
    if c1 == c2:
        param["strTHydrogen"] += temp

    # ---charge layer
    c2 = "p"
    if c1 == c2:
        param["strProton"] += temp
    c2 = "q"
    if c1 == c2:
        param["strTCharge"] += temp

    # --stereotypes
    c2 = "b"
    if c1 == c2:
        param["strBond"] += temp
    c2 = "m"
    if c1 == c2:
        param["strTetra"] += temp
    c2 = "t"
    if c1 == c2:
        param["strTetra"] += temp
    c2 = "s"
    if c1 == c2:
        param["strStereo"] += temp

    if param["strTAtom"] == "":
        param["strTAtom"] = "c1"

param["iNumSubStruc"] = param["strTFormula"].count(".") + 1

iNumSub1 = 0
iNumSub2 = 0
iNumSub3 = 0
iNumSub4 = 0
iMul = 0

# iNumSub1 = param["strTFormula"].count(".") + 1
# if param["strTAtom"] != "":
#     iNumSub2 = param["strTAtom"].count(";")
# if param["strTHydrogen"] != "":
#     iNumSub2 = param["strTHydrogen"].count(";")
# if param["strTCharge"] != "":
#     iNumSub2 = param["strTCharge"].count(";")


str1 = ""
str2 = ""
cTemp = ""
k = 0
for i in range(0, param["iNumSubStruc"]):
    str1 = ""
    k = 0
    param["strSubFormula"][i] = param["strTFormula"].split(".")[i]
    while True:
        cTemp = param["strSubFormula"][i][k]

        if ord(cTemp) > 47 and ord(cTemp) < 58:
            str1 = str1 + cTemp
            k = k + 1
        else:
            break

    if str1 == "":
        iMul = 0
    else:
        iMul = int(str1) - 1
    str2 = str1 + "*"

    param["iMulSub"][i] = iMul + 1
    param["strSubConnect"][i] = param["strTAtom"].split(";")[i + iNumSub2]
    param["strSubHydrogen"][i] = param["strTHydrogen"].split(";")[i + iNumSub3]
    param["strSubCharge"][i] = param["strTCharge"].split(";")[i + iNumSub4]

    if param["strSubConnect"][i] == "" or param["strSubConnect"][i] == "c":
        iNumSub2 = iNumSub2 + iMul
    if param["strSubHydrogen"][i] == "" or param["strSubHydrogen"][i] == "h":
        iNumSub3 = iNumSub3 + iMul
    if param["strSubCharge"][i] == "" or param["strSubCharge"][i] == "q":
        iNumSub4 = iNumSub4 + iMul

print("iMulSub", ":", param["iMulSub"])
print("strSubFormula", ":", param["strSubFormula"])
print("strSubConnect", ":", param["strSubConnect"])
print("strSubHydrogen", ":", param["strSubHydrogen"])
print("strSubCharge", ":", param["strSubCharge"], "\n")

for i in range(1, param["iNumSubStruc"]):
    if param["strSubFormula"][i] != "":
        if param["strSubFormula"][i].count("*") > 0:
            str1 = param["strSubFormula"][i].split("*")[1]
            param["strSubFormula"][i] = str1
    if param["strSubConnect"][i] != "":
        if param["strSubConnect"][i].count("*") > 0:
            str1 = param["strSubConnect"][i].split("*")[1]
            param["strSubConnect"][i] = str1
    if param["strSubHydrogen"][i] != "":
        if param["strSubHydrogen"][i].count("*") > 0:
            str1 = param["strSubHydrogen"][i].split("*")[1]
            param["strSubHydrogen"][i] = str1
    if param["strSubCharge"][i] != "":
        if param["strSubCharge"][i].count("*") > 0:
            str1 = param["strSubCharge"][i].split("*")[1]
            param["strSubCharge"][i] = str1

print("strSubFormula", ":", param["strSubFormula"])
print("strSubConnect", ":", param["strSubConnect"])
print("strSubHydrogen", ":", param["strSubHydrogen"])
print("strSubCharge", ":", param["strSubCharge"], "\n")
