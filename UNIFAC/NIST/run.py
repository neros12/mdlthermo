import sys
import json

from modules import cal_activity_coefficient


if __name__ == "__main__" and len(sys.argv) > 2:
    try:
        SMILES1 = str(sys.argv[1])
        SMILES2 = str(sys.argv[2])
        x1 = float(sys.argv[3])
        x2 = float(sys.argv[4])
        T = float(sys.argv[5])
        try:
            activity_coefficient = cal_activity_coefficient(
                SMILES1,
                SMILES2,
                x1,
                x2,
                T,
            )
            print(
                json.dumps(
                    {
                        "success": True,
                        "SMILES1": SMILES1,
                        "SMILES2": SMILES2,
                        "x1": x1,
                        "x2": x2,
                        "gamma1": activity_coefficient[0],
                        "gamma2": activity_coefficient[1],
                        "temperature": T,
                    }
                )
            )
        except Exception as error_message:
            print(
                json.dumps(
                    {
                        "success": False,
                        "SMILES1": SMILES1,
                        "SMILES2": SMILES2,
                        "message": error_message,
                    }
                )
            )
    except:
        print(json.dumps({"success": False, "message": "wrong input values"}))
