import numpy as np
import matplotlib.pyplot as plt

def parse_xpm(xpm_file):
    """Parses a GROMACS XPM matrix into a 2D NumPy array."""
    color_map = {}
    data = []
    in_matrix = False
    
    with open(xpm_file, 'r') as f:
        for line in f:
            if 'static char *' in line:
                in_matrix = True
                continue
            if not in_matrix:
                continue
            
            line = line.strip()
            # Look for lines containing matrix strings
            if line.startswith('"') and not line.startswith('/*'):
                row_str = line.split('"')[1]
                
                # 1. Extract color mappings
                if 'c #' in line or 'c None' in line:
                    char = row_str[0]
                    # Map background (None/White) to 0, and bonds (Red/Black) to 1
                    color_map[char] = 0 if ('None' in line or '#FFFFFF' in line) else 1
                    
                # 2. Skip the dimensions line (e.g., "10000 16 2 1")
                elif len(row_str.split()) >= 3 and row_str.replace(" ", "").isdigit():
                    continue 
                    
                # 3. Extract matrix rows
                else:
                    data.append([color_map.get(c, 0) for c in row_str])
                    
    # Enforce uniform row length to prevent inhomogeneous shape errors
    if not data:
        raise ValueError("No data parsed from the XPM file.")
    max_len = max(len(row) for row in data)
    cleaned_data = [row for row in data if len(row) == max_len]
    
    return np.array(cleaned_data)

# Parse matrix
matrix = parse_xpm('hbmap.xpm')

# Calculate occupancy percentage for each hydrogen bond
occupancy = (np.sum(matrix, axis=1) / matrix.shape[1]) * 100

# Plot Existence Map
plt.figure(figsize=(10, 6))
plt.imshow(matrix, aspect='auto', cmap='binary', origin='lower', interpolation='none')
plt.xlabel('Frame / Time Step')
plt.ylabel('Hydrogen Bond Index')
plt.title('Hydrogen Bond Existence Map')
plt.colorbar(label='0 = Absent, 1 = Present')
plt.tight_layout()
plt.savefig('hbmap_plot.png', dpi=300)
print("Saved plot as hbmap_plot.png")

# Print the most persistent hydrogen bonds (up to 5)
top_bonds = np.argsort(occupancy)[::-1][:min(5, len(occupancy))]
print("\nTop Persistent Hydrogen Bonds:")
for idx in top_bonds:
    print(f"HBond Index {idx}: {occupancy[idx]:.2f}% occupancy")