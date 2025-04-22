import os
import sys
from dotenv import load_dotenv

print("Before loading .env:")
print(f"PROJECT_FUNCTIONS_PATH: {os.getenv('PROJECT_FUNCTIONS_PATH')}")

load_dotenv()

print("\nAfter loading .env:")
print(f"PROJECT_FUNCTIONS_PATH: {os.getenv('PROJECT_FUNCTIONS_PATH')}")

if os.getenv('PROJECT_FUNCTIONS_PATH'):
    sys.path.insert(0, os.getenv('PROJECT_FUNCTIONS_PATH'))
    print(f"\nAdded to sys.path: {os.getenv('PROJECT_FUNCTIONS_PATH')}")
    print(f"Current sys.path: {sys.path}")
    
    # Check if the directory exists
    if os.path.exists(os.getenv('PROJECT_FUNCTIONS_PATH')):
        print(f"\nDirectory exists: {os.getenv('PROJECT_FUNCTIONS_PATH')}")
        
        # List files in the directory
        print("\nFiles in directory:")
        for file in os.listdir(os.getenv('PROJECT_FUNCTIONS_PATH')):
            print(f"  - {file}")
            
        # Check if grn_helpers.py exists
        grn_helpers_path = os.path.join(os.getenv('PROJECT_FUNCTIONS_PATH'), 'grn_helpers.py')
        if os.path.exists(grn_helpers_path):
            print(f"\ngrn_helpers.py exists: {grn_helpers_path}")
            
            # Try to import
            try:
                sys.path.insert(0, os.getenv('PROJECT_FUNCTIONS_PATH'))
                print("\nTrying to import grn_helpers...")
                import grn_helpers
                print("Successfully imported grn_helpers")
            except Exception as e:
                print(f"Error importing grn_helpers: {e}")
        else:
            print(f"\ngrn_helpers.py does NOT exist at: {grn_helpers_path}")
    else:
        print(f"\nDirectory does NOT exist: {os.getenv('PROJECT_FUNCTIONS_PATH')}")
else:
    print("\nPROJECT_FUNCTIONS_PATH is not set after loading .env file") 