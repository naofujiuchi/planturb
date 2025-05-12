#%%
from trimesh.creation import cylinder, icosphere, box
from trimesh.util import concatenate
import numpy as np

# blockMesh is 1.8m x 1.8m x 2m

# Parameters in real 
small_leaf_length = 0.06  # Semi-major axis length (small leaf length) [m]
small_leaf_width = 0.03  # Semi-minor axis length (small leaf width) [m]
a = small_leaf_length / 2  # Semi-major axis length (small leaf length) [m]
b = small_leaf_width / 2  # Semi-minor axis length (small leaf width) [m]
small_leaf_area = a * b * np.pi # leaf area of the small leaf [m2]

# Small leaf (ellipse)
# Create a circular cylinder
cyl = trimesh.creation.cylinder(radius=small_leaf_width, height=0.02, sections=64)
# Scale in x and y to make it elliptical (e.g., a=2, b=1)
scale_matrix = np.diag([2.0, 1.0, 1.0, 1.0])  # a=2, b=1
cyl.apply_transform(scale_matrix)

# Leaf with 9 small leaves
# Small leaf locations: (0.05, -0.09), (0.05, 0.09), (0.10, -0.09), (0.10, 0.09), (0.15, -0.09), (0.15, 0.09), (0.20, -0.09), (0.20, 0.09), (0.25, 0)
leaf = None
small_leaf_locations = [
    [0.06, -0.08, 0], # 1
    [0.06, 0.08, 0], # 2
    [0.12, -0.12, 0], # 3
    [0.12, 0.12, 0], # 4
    [0.18, -0.10, 0], # 5
    [0.18, 0.10, 0], # 6
    [0.24, -0.08, 0], # 7
    [0.24, 0.08, 0], # 8
    [0.3, 0, 0], # 9
]
# Create a rotation matrix (e.g., 45 degrees around Z axis)
angle = np.deg2rad(90)
rotation_matrix = trimesh.transformations.rotation_matrix(
    angle, [0, 0, 1], point=[0, 0, 0]
)
for location in small_leaf_locations:
    cyl_copy = cyl.copy()
    if location != small_leaf_locations[-1]:
        cyl_copy.apply_transform(rotation_matrix) # Apply the rotation
    cyl_copy.apply_translation(location)
    leaf = concatenate([leaf, cyl_copy])

# leaf.show()
# leaf.export("tomato_leaf.stl")

# Tomato plant
# Phyllotaxy: 137.5 degrees
# Leaf locations: z = 0.06, 0.12, 0.18, 0.24, 0.3, 0.36, 0.42, 0.48, 0.54, 0.60, 0.66, 0.72
z_locations = [0.06, 0.12, 0.18, 0.24, 0.3, 0.36, 0.42, 0.48, 0.54, 0.60, 0.66, 0.72]  # in meters
phyllotaxy_angle_deg = 137.5
# List to hold all transformed leaves
leaves = []
for i, z in enumerate(z_locations):
    # Compute rotation angle for this leaf
    angle_deg = i * phyllotaxy_angle_deg
    angle_rad = np.deg2rad(angle_deg)    
    # Create a copy of the leaf
    leaf_copy = leaf.copy()    
    # Rotate around Z axis (stem axis)
    rotation_matrix = trimesh.transformations.rotation_matrix(
        angle_rad, [0, 0, 1], point=[0, 0, 0]
    )
    leaf_copy.apply_transform(rotation_matrix)
    # Translate to the correct z location
    translation_matrix = trimesh.transformations.translation_matrix([0, 0, z])
    leaf_copy.apply_transform(translation_matrix)
    # Add to the list
    leaves.append(leaf_copy)
# Combine all leaves into one mesh
plant = trimesh.util.concatenate(leaves)

# Rotate 90 degrees around y axis
plant.apply_transform(trimesh.transformations.rotation_matrix(np.deg2rad(90), [0, 1, 0], point=[0, 0, 0]))
plant.apply_translation([0.2, 0.9, 1.0])

# # Optionally, show the plant
# plant.show()
# Save as STL
plant.export("tomato_plant.stl")
