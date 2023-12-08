import sys, os

sys.path.append((os.path.dirname(__file__)))
from modules import NIST_UNIFAC

if len(sys.argv) < 2:
    result = NIST_UNIFAC.cal_activity_coefficient(
        "c1ccccc1",
        "CCCO",
        0.3,
        0.7,
        303,
    )
    print("a12 : ", result[0])
    print("a21 : ", result[1])
else:
    try:
        SMILES1 = str(sys.argv[1])
        SMILES2 = str(sys.argv[2])
        x1 = float(sys.argv[3])
        x2 = float(sys.argv[4])
        T = float(sys.argv[5])

        result = NIST_UNIFAC.cal_activity_coefficient(
            SMILES1,
            SMILES2,
            x1,
            x2,
            T,
        )
        print(f"Result:{result[0]},{result[1]}", end="")
    except Exception as ex:
        print(f"Error:{ex}", end="")
