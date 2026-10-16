import sys
from pathlib import Path

def process_file(filepath):
    try:
        with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
            lines = f.readlines()

        new_lines = []
        for line in lines:
            # If the line starts with 6 spaces, replace them with 4 spaces
            if line.startswith('      '):
                new_lines.append('    ' + line[6:])
            # If the line starts with '!', remove up to two spaces right after it
            elif line.startswith('!'):
                rest = line[1:]
                removed = 0
                while removed < 2 and rest.startswith(' '):
                    rest = rest[1:]
                    removed += 1
                new_lines.append('!' + rest)
            else:
                new_lines.append(line)

        with open(filepath, 'w', encoding='utf-8') as f:
            f.writelines(new_lines)
        print(f"Updated: {filepath}")
    except Exception as e:
        print(f"Error processing {filepath}: {e}", file=sys.stderr)

def main():
    current_dir = Path('.')
    # Recursively find all .f90 files (case-insensitive search can also be handled if needed)
    f90_files = list(current_dir.rglob('*.f90')) + list(current_dir.rglob('*.F90'))

    if not f90_files:
        print("No .f90 files found in the current or subdirectories.")
        return

    for filepath in f90_files:
        process_file(filepath)

    print(f"\nDone! Processed {len(f90_files)} file(s).")

if __name__ == '__main__':
    main()