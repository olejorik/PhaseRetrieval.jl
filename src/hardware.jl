# dictionaries with typical harware for ease of use with experimental data

camerasdict = Dict()
lensesdict = Dict()


# Cameras
camerasdict["UI1490"] = PhaseRetrieval.CameraChip(pixelsize = 1.67E-3, imagesize = (3840, 2748), bitdepth = 8)
camerasdict["UI1540"] = PhaseRetrieval.CameraChip(pixelsize = 5.2E-3, imagesize = (1280, 1024), bitdepth = 8)
camerasdict["UI2210"] = PhaseRetrieval.CameraChip(pixelsize = 9.9E-3, imagesize = (640, 480), bitdepth = 8)

# Lenses
lensesdict["F300A25"] = PhaseRetrieval.ImagingLens(300, 25)
lensesdict["F700A75"] = PhaseRetrieval.ImagingLens(700, 75)