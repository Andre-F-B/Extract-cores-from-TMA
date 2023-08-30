# Extracts cores from Tissue microarrays (TMA files) with their respective annotations

# from __future__ import annotations
import os
import openslide
from tqdm.notebook import tqdm
from shapely.geometry import shape, Polygon, Point
import numpy as np
import cv2
from pathlib import Path
from paquo.projects import QuPathProject
from tqdm.autonotebook import tqdm
from distutils.util import strtobool


# First of all, your image must be part of a QuPath project.
# Then you must define the cores using the TMA dearrayer tool
# (what matters most is the center of each core, the radius makes no difference since it is defined here in the code)
# Export the TMA measurements as a txt file.
# The annotations must also be saved on the image inside your project.


# ------------------------------------------ Parameters to edit -------------------------------------------------------
annotated = True # if you want to export the cores with their annotations
ignore_unannotated = True # if you would like to skip/ignore cores with no annotations
# TMA measurements as txt
txt_filename = "path/to/txt_file.txt"
# Path to the project
PROJECT_FILE = Path(f"project.qpproj")
PROJECT_PATH = Path(os.path.join(r"Path/to/project", f"{PROJECT_FILE}"))
# Image index on the project (in the list "QuPathProject(PROJECT_PATH, mode='r').images")
image_index = 2
output_dir = Path('path/to/output/dir')
tma_file = Path('path/to/tma_file.mrxs')
level = 0
radius_pixels = 5500 # in pixels
tmaspot_size = 2*radius_pixels # in pixels
# ---------------------------------------------------------------------------------------------------------------------


slide = openslide.OpenSlide(tma_file)
mpp_ratio = float(slide.properties['openslide.mpp-y'])
radius_microns = radius_pixels*float(mpp_ratio)
print(f'Core radius in microns: {radius_microns}')

# Now we define the dataset with the centroids of each core, read the annotations, and convert them to GeoJSON.

# Getting the TMA information
dataset = np.loadtxt(txt_filename, dtype=str, skiprows=1)

# Getting the annotations
qp = QuPathProject(PROJECT_PATH, mode='r')
print(f"Opened project ‘{qp.name}’ ")
print(f"Project has {len(qp.images)} image(s).")
print(qp.images)
image = qp.images[image_index]
tma_filename = image.image_name[:-5]
annotations = image.hierarchy.annotations # annotations drawn by pathologist
annotations_geojson = image.hierarchy.to_geojson()

# We define bounding boxes for each core based on their centroid and on radius_microns

core_centroids = []

for row in tqdm(dataset):
    fname, object_id, label, missing, x, y = row
    core_centroids.append((float(x), float(y)))

# print(len(core_centroids))

core_points = [Point(center) for center in core_centroids]
# print(len(core_points))

core_boxes = [Point(center).buffer(radius_microns) for center in core_centroids]
# print(len(core_boxes))

# The core coordinates are in micrometers and the annotation coordinates are in pixels, so we convert the coordinates of the cores to pixels

# Iterate over each core box and divide its coordinates by the mpp ratio
for i, core_box in enumerate(core_boxes):
    # Get the current coordinates of the core box
    coordinates = list(core_box.exterior.coords)

    # Divide each coordinate by mpp_ratio
    updated_coordinates = [(x / mpp_ratio, y / mpp_ratio) for x, y in coordinates]

    # Update the core box with the new coordinates
    core_box = Polygon(updated_coordinates)
    core_boxes[i] = core_box


# We convert the annotations into shapely objects and detect which annotations fall within each core

shapely_annotations = [shape(annotation['geometry']) for annotation in annotations_geojson]

# recursive function to convert each geometry inside possible MultiPolygons into a Polygon
def breakdown_multipolygons(shapely_annotations):

    count_multi = 0
    count_something_else = 0
    for index, annotation in enumerate(shapely_annotations):
        if annotation.geom_type != 'Polygon':
            if annotation.geom_type == 'MultiPolygon':
                count_multi += 1
                del shapely_annotations[index]
                for geometry in annotation:
                    subshape = shape(geometry)
                    shapely_annotations.append(subshape)
            else:
                count_something_else += 1

    if count_multi == 0 and count_something_else == 0:
        # everything is a polygon
        return None
    
    if count_something_else != 0:
        print(f'There are {count_something_else} unidentified shapes here')
        return None

    # checking if there is another weird shape
    count = 0
    for index, annotation in enumerate(shapely_annotations):
        if annotation.geom_type != 'Polygon':
            count += 1
            print(f'Possible problematic shape found: "{annotation.geom_type}"')
    if count == 0:
        print('All shapes are Polygons')
    else:
        breakdown_multipolygons(shapely_annotations)


breakdown_multipolygons(shapely_annotations)

annotations_in_cores = [[] for _ in core_points]  # Empty lists to store annotations for each core

for annotation in shapely_annotations:
    for i, core_box in enumerate(core_boxes):
        if core_box.contains(annotation):
            annotations_in_cores[i].append(annotation)
            break  # Skip checking other cores if annotation is already assigned to one


# Exports each individual core in a TMA file to a TIFF image


print(f"Extracting cores from {tma_file} into {output_dir}")


if not os.path.isdir(f"{output_dir}"):
    os.mkdir(f"{output_dir}")


# Print the slide's downsample info
level_dims = slide.level_dimensions[level]
print(f'Slide level dimensions: {slide.level_dimensions}')
level_downsample = slide.level_downsamples[level]
print(f'Downsample at level {level} is: {level_downsample}')
print(f'WSI dimensions at level {level} are: {level_dims}.')

bounds_x = float(slide.properties['openslide.bounds-x']) if ("openslide.bounds-x") in slide.properties.keys() else 0
bounds_y = float(slide.properties['openslide.bounds-y']) if ("openslide.bounds-y") in slide.properties.keys() else 0

ratio_x = 1.0/float(slide.properties['openslide.mpp-x'])
ratio_y = 1.0/float(slide.properties['openslide.mpp-y'])
print(slide.properties['openslide.mpp-x'])
print(f'Ratios: {(ratio_x, ratio_y)}')
print(f'Bounds: {(bounds_x, bounds_y)}')

dataset = np.loadtxt(txt_filename, dtype=str, skiprows=1)
print(f"Number of rows in txt file ：{len(dataset)}")  # note that you aren't guaranteed to get exactly this many spots out, simply because some may be set to false or are missing


i = 0
for row in tqdm(dataset):
    fname, object_id, label, missing, x, y=row
    
    # skip unnanotated cores
    if ignore_unannotated and (not annotations_in_cores[i]):
        i += 1
        continue

    if(not strtobool(missing)):
        x = (float(x)*ratio_x) + bounds_x
        y = (float(y)*ratio_y) + bounds_y
        print(f"Extracting spot {label} at location", (x, y))
        scaled_spot_size = (int(tmaspot_size/level_downsample), int(tmaspot_size/level_downsample))
        top_left_x, top_left_y = int(x - tmaspot_size*0.5), int(y-tmaspot_size*0.5)
        tmaspot = np.asarray(slide.read_region((top_left_x, top_left_y), level, scaled_spot_size))[:, :, 0:3]
        tmaspot = cv2.cvtColor(tmaspot,cv2.COLOR_RGB2BGR)

        if annotated:
            # draws each annotation on the extracted patch
            translated_annotations = []
            for annotation in annotations_in_cores[i]:
                annotation_array = np.expand_dims(np.array(annotation.exterior.coords, np.int32), axis=0)
                translated_annotation = annotation_array - [top_left_x, top_left_y] + [bounds_x, bounds_y]
                translated_annotation = translated_annotation.astype(int)
                translated_annotations.append(translated_annotation)
            img_mod = cv2.polylines(tmaspot, translated_annotations, True, (0, 0, 255),3)

        
        cv2.imwrite(f"{output_dir}/{tma_filename}_({label}).tiff", tmaspot)
    else:
        print(f'The spot {label} is missing, skipping!')

    i += 1

print('Extracted all the spots!')



