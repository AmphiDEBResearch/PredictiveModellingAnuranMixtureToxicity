# ======================================== #
# helper script to convert figures to tiff
# ======================================== #

using Pkg; Pkg.activate(".")
using FileIO, ImageIO

Pkg.add("ImageIO")

root = "manuscript/v01d_first_submission"

files = readdir(root)
pngs = filter(f->occursin(".png", f), files)

for png in pngs
    name = split(png, ".")[1]
    tiffname = join([name, ".tiff"])
    save(joinpath(root, tiffname), load(joinpath(root, png)))
end

