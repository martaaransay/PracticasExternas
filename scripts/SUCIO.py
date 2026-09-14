import os
actual_dir = os.path.dirname(__file__)
HATTCI_DIR = os.path.join(actual_dir, "T3_integronfiltering", "HattCI")

os.environ["PATH"] = f"{HATTCI_DIR}{os.pathsep}{os.environ['PATH']}"
REPO_ROOT = os.path.join(actual_dir, os.pardir, os.pardir)

print(REPO_ROOT)
print(HATTCI_DIR)