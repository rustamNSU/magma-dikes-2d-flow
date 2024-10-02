python
import sys
import gdb
sys.path.insert(0, '/home/rustam/repositories/magma-dikes-2d-flow/.vscode')
from eigengdb import register_eigen_printers
register_eigen_printers (None)

from nlohmann_json import json_lookup_function
gdb.pretty_printers.append(json_lookup_function)
end

# source /home/rustam/repositories/magma-dikes-2d-flow/.vscode/nlohmann_json.py