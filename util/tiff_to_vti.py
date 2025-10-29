import rasterio
from vtk.util import numpy_support
import vtk
import numpy as np

in_file  = "../data/Klimadatastyrelsen/DanmarksDybdeModel/2024/ddm_50m.dybde.tiff"
out_file = "dem.vti"
print("Reading →", in_file)
with rasterio.open(in_file) as src:
    data = src.read(1).astype(np.float32)
    height, width = data.shape
    transform = src.transform

print("Converting to VTI...")
image = vtk.vtkImageData()
image.SetDimensions(width, height, 1)
image.AllocateScalars(vtk.VTK_FLOAT, 1)

# numpy → vtk
print("Filling VTI data...")
vtk_array = numpy_support.numpy_to_vtk(
    data.ravel(order="C"),
    deep=True,
    array_type=vtk.VTK_FLOAT
)
image.GetPointData().SetScalars(vtk_array)

# georeference
print("Setting georeference...")
image.SetOrigin(transform.c, transform.f, 0)
image.SetSpacing(transform.a, transform.e, 1)

print("Writing →", out_file)

writer = vtk.vtkXMLImageDataWriter()
writer.SetFileName(out_file)

writer.SetDataModeToBinary()
writer.EncodeAppendedDataOn()
#writer.SetDataModeToAppended()
#writer.EncodeAppendedDataOff()

writer.SetInputData(image)
writer.Write()

print("OK →", out_file)
