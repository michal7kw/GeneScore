import warnings
import pandas as pd

# Filter out the specific FutureWarning from pandas
warnings.filterwarnings('ignore', category=FutureWarning,
                       message='The default of observed=False is deprecated')

# Import your original script
import sys
import os

# Get the directory of the current script
current_dir = os.path.dirname(os.path.abspath(__file__))

# Add the directory to Python path
sys.path.append(current_dir)

# Now import and run your original script
import importlib.util
spec = importlib.util.spec_from_file_location(
    "original_script",
    os.path.join(current_dir, "2_herring_atac_processing.py")
)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
