import sys
import json

from _fragment import get_feature_matrces


nfm, efm = get_feature_matrces(sys.argv[1])
print(json.dumps({"nfm": nfm.tolist(), "efm": efm.tolist()}, indent=2), end="")
