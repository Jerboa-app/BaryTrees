include("BaryTrees/BaryTrees.jl")
using Main.BaryTrees, Plots

# construct an empty tree from a Barycentric coordinate system defined
# by the triangle with vertices (0,sqrt(3)/2), (-0.5,0), (0.5,0)
T = BaryTree(
    BaryCentric(Point((0,sqrt(3)/2)),Point((-0.5,0)),Point((0.5,0.0)))
)

# random points in a square: (-1,-1), (-1,1), (1,1), (1,-1)
p = (rand(2^9,2) .- 0.5).*2;

pin = zeros(1,2)
pout = zeros(1,2)

for i in 1:size(p,1)
    b = Cartesian2BaryCentric(Point((p[i,:]...)),T.boundary.b)

    # test if in trig
    if b[1] >= 0 && b[2] >= 0 && b[1]+b[2] < 1
        pin = cat(pin,p[i:i,:],dims=1)
    else
        pout = cat(pout,p[i:i,:],dims=1)
    end
end

# insert only points inside trig
for i in 1:size(pin,1)
    insert!(Point((pin[i,:]...)),T)
end

scatter(pin[:,1],pin[:,2],label="",dpi=300)
Draw!(T)
plot!()

savefig("barytree.png")
