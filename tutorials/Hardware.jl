# ```@meta
# CurrentModule = PhaseRetrieval
# DocTestSetup = quote
#     using PhaseRetrieval
# end
# ```

# # Hardware types

# ## Camera chip
# [`CameraChip`](@ref) represents a device that outputs sampled field intensity, quantized 
# and packed in 8- or 16-bit integer per color channel (currently, only monochrome 
# cameras are implemented).

# This code sets up a camera with 5.2μm pixel, 1280 × 1024 frame, bitdepth and channeldepth
#  of 8 bits:

using PhaseRetrieval
using PhasePlots
using CairoMakie
using ImageCore
cam = CameraChip(;
    pixelsize=5.2um, imagesize=(1280, 1024), bitdepth=8, channelbitdepth=8
)


# Cameras can be also created using [`camerasdict`](@ref) dictionary.
# Let's consider examples of 8-bit and 12-bit cameras:
cam8bit = camerasdict["UI1240"]

#  and
cam12bit = camerasdict["UI3860"]

#  and compare the results of using thes two cameras to measure this 
#  field created by a circular aperture:
ap,_ = PhaseRetrieval.aperture(-1:.025:1,-1:.025:1,0.9)

showarray(ap)

# The field itself is represented by complex array, which has this amplitude and phase (we
# can see that due to the aperture symmentry, the field in the focal plane is real, actually,
# and the phase take only values 0 and π):
imf = PhaseRetrieval.toimageplane(ap)
showarray(abs.(imf)) |> display
showphasetight(angle.(imf)); current_figure()

# And now we can quantize this field by our 8-bit and 12-bit example cameras
cam8bit(imf) |> collect 
# 
cam12bit(imf) |> collect

# !!! note
# We have used `collect` here, as `CameraChip` returns a mapped view of the field.

# We can see that for our 12-bit camera, the measurements are represented by the last 
# 12 bits of a 16-bit number
# (`ImageCore.N4f12`(@ref)).

# We can hardly see the difference in the camera outputs. Compate 8bit:
cam8bit(imf) .|> Gray
# with 12 bit:
cam12bit(imf) .|> Gray

# The difference is  visible in logarithmic scale. compare 8bit:
cam8bit(imf) .|> float |> PhaseRetrieval.logrescale .|> Gray
#  and 12 bit
cam12bit(imf) .|> float |> PhaseRetrieval.logrescale .|> Gray


# !!! note
# Here we have used already sampled optical filed `imf`, so there is no difference in 
#  the sampling rate performed by the cameras. This effect will be demonstrated in the next example.
