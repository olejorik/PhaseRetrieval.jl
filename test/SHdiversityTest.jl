using PhaseRetrieval
using Plots, Images
default(aspect_ratio =:equal)

# repeating my simulations in MMA
dim =520 # total grid size for padding 
s=5.2
p=150
# cell=p/s
cell = 64
d=1.2
r=cell*8/dim # aperture of 512 pixels, relative size

ph1 = SHdiversity(dim, cell, (12,-8.5), .1); # offset and rotation are chosen arbitrary
ph1 = 2*pi*ph1./(4* dim)


heatmap(phwrap.(ph1), c=:cyclic_mygbm_30_95_c78_n256, size=(dim,dim)) |> display

ap1 = [ i^2 + j^2 < r^2 for i in range(-1,1,length= dim), j in range(-1,1, length=dim)]
display(heatmap(ap1))
hp1ref = psf(ap1, ph1)
heatmap(hp1ref, c=:grays, size = (540,540), cb=:none)  |> display

# using Images instead of heatmap
Gray.(PhaseRetrieval.rescale(phwrap.(ph1))) |> display
Gray.(PhaseRetrieval.rescale(hp1ref)) |> display

# for 4 times more lenslets (32 x 32) we need upsample it 4 times
dim=1080
q = 4
dimup = dim*q
phup = SHdiversity(dimup, cell, (12,-8.5), .1); # offset and rotation are chosen arbitrary
phup = 2*pi*phup./(2* dimup)
Gray.(PhaseRetrieval.rescale(phwrap.(phup))) |> display
apup = [ i^2 + j^2 < r^2 for i in range(-1,1,length= dimup), j in range(-1,1, length=dimup)]
hpupref = imresize(psf(apup, phup*2), ratio=1/q)
Gray.(PhaseRetrieval.rescale(hpupref)) |> display
hp2=restrict(restrict(psf(apup, -phup*2)))
Gray.(PhaseRetrieval.rescale(hp2)) |> display

# in MMA I had dimgrid 5850, q = 5, ap= 4805, number of pixel per s/a is 145
dim=1170
q = 5
dimup = dim*q
r= 961/1170
cellup =145
phup = SHdiversity(dimup, cellup, (12,-8.5), .1); # offset and rotation are chosen arbitrary
phup = 2*pi*phup./(2* dimup)
Gray.(PhaseRetrieval.rescale(phwrap.(phup))) |> display
apup = [ i^2 + j^2 < r^2 for i in range(-1,1,length= dimup), j in range(-1,1, length=dimup)]
hpupref = imresize(psf(apup, phup), ratio=1/q)
Gray.(PhaseRetrieval.rescale(hpupref)) |> display
save("p150f10_ref.png",Gray.(PhaseRetrieval.rescale(hpupref))) 