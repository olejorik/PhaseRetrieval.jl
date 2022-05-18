# dictionaries with typical harware for ease of use with experimental data

camerasdict = Dict()
lensesdict = Dict()


# Cameras
camerasdict["UI1490"] = PhaseRetrieval.CameraChip(pixelsize = 1.67E-3, imagesize = (3840, 2748), bitdepth = 8, channelbitdepth = 8)
camerasdict["UI1540"] = PhaseRetrieval.CameraChip(pixelsize = 5.2E-3, imagesize = (1280, 1024), bitdepth = 8, channelbitdepth = 8)
camerasdict["UI2210"] = PhaseRetrieval.CameraChip(pixelsize = 9.9E-3, imagesize = (640, 480), bitdepth = 8, channelbitdepth = 8)
camerasdict["UI3260"] = PhaseRetrieval.CameraChip(pixelsize = 5.86E-3, imagesize = (1936, 1216), bitdepth = 12, channelbitdepth = 16) #Bit depth of the images is 16

# Lenses
lensesdict["F300A25"] = PhaseRetrieval.ImagingLens(300, 25)
lensesdict["F700A75"] = PhaseRetrieval.ImagingLens(700, 75)