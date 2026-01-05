import subprocess
from pathlib import Path
import platform

def main():
    if platform.system() != "Darwin":
        print("This command is only needed on macOS.")
        return

    pkg_root = Path(__file__).resolve().parent
    bin_dir = pkg_root / "wrapper" / "bin"

    if not bin_dir.exists():
        print(f"Binary directory not found: {bin_dir}")
        return

    print(f"Removing quarantine attribute from: {bin_dir}")
    subprocess.run(
        ["xattr", "-dr", "com.apple.quarantine", str(bin_dir)],
        check=False
    )

    print("Done. You should now be able to run the executables.")

if __name__ == "__main__":
    main()
