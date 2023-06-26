# dictionaries with typical harware for ease of use with experimental data

const m = 1
const mm = 1e-3
const um = 1e-6
const μm = 1e-6
const nm = 1e-9

"""
    `camerasdict` is a dictionary with often used cameras.

Use `keys(camerasdict)` to get the list of implemented cameras.
"""
camerasdict = Dict()

"""
    `lensesdict` is a dictionary with often used cameras.

Use `keys(lensesdict)` to get the list of implemented cameras.
"""
lensesdict = Dict()


# Cameras
camerasdict["UI1490"] = PhaseRetrieval.CameraChip(pixelsize = 1.67um, imagesize = (3840, 2748), bitdepth = 8, channelbitdepth = 8)
camerasdict["UI1540"] = PhaseRetrieval.CameraChip(pixelsize = 5.2um, imagesize = (1280, 1024), bitdepth = 8, channelbitdepth = 8)
camerasdict["UI2210"] = PhaseRetrieval.CameraChip(pixelsize = 9.9um, imagesize = (640, 480), bitdepth = 8, channelbitdepth = 8)
camerasdict["UI3260"] = PhaseRetrieval.CameraChip(pixelsize = 5.86um, imagesize = (1936, 1216), bitdepth = 12, channelbitdepth = 16) #Bit depth of the images is 16
camerasdict["UI3860"] = PhaseRetrieval.CameraChip(pixelsize = 2.9um, imagesize = (1936, 1096), bitdepth = 12, channelbitdepth = 16) #Bit depth of the images is 16
camerasdict["MC203MG"] = PhaseRetrieval.CameraChip(pixelsize = 2.74um, imagesize = (4504, 4504), bitdepth = 12, channelbitdepth = 16) #Bit depth of the images is 16

# Lenses
lensesdict["F300A25"] = PhaseRetrieval.ImagingLens(300mm, 25mm)
lensesdict["F700A75"] = PhaseRetrieval.ImagingLens(700mm, 75mm)



export camerasdict, lensesdict, m, mm, um, μm, nm