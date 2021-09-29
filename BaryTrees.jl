module BaryTrees

    using Plots, LinearAlgebra

    import Base.∈, Base.∉, Plots.Shape

    export BaryTree, Point, insert!, ∈, ∉, Draw, Draw!, Query

    const Point = Tuple{Float64,Float64}

    Point(x::Float64,y::Float64) = Point((x,y))

    struct BaryCentric
        """
            Represents the 3 cartesian (2D) points of a Barycentric
            coordinate system
        """
        r₁::Point
        r₂::Point
        r₃::Point
    end

    function BaryCentric2Cartesian(λ₁,λ₂,λ₃,p::BaryCentric)
        return Point(
            λ₁*p.r₁[1]+λ₂*p.r₂[1]+λ₃*p.r₃[1],
            λ₁*p.r₁[2]+λ₂*p.r₂[2]+λ₃*p.r₃[2],
        )
    end

    function Cartesian2BaryCentric(p::Point,b::BaryCentric)
        x,y = p
        x1,x2,x3 = b.r₁[1],b.r₂[1],b.r₃[1]
        y1,y2,y3 = b.r₁[2],b.r₂[2],b.r₃[2]

        l1 = ( (y2-y3)*(x-x3)+(x3-x2)*(y-y3) ) / ( (y2-y3)*(x1-x3) + (x3-x2)*(y1-y3) )
        l2 = ( (y3-y1)*(x-x3) + (x1-x3)*(y-y3) ) / ( (y2-y3)*(x1-x3) + (x3-x2)*(y1-y3) )
        return (
            l1,l2,1.0-l1-l2,
        )
    end

    function Medians(b::BaryCentric)::Vector{Point}
        """
            Returns the three mid points of each triangular edge of
            the triangle defining b in cartesian units
        """
        return [
            BaryCentric2Cartesian(1/2,0,1/2,b)
            BaryCentric2Cartesian(0,1/2,1/2,b)
            BaryCentric2Cartesian(1/2,1/2,0,b)
        ]
    end

    struct AxisAlignedTriangle
        """
            wrapper around BaryCentric
        """
        b::BaryCentric
    end

    Medians(t::AxisAlignedTriangle)::Vector{Point} = Medians(t.b)

    mutable struct BaryTree
        root::Bool
        boundary::AxisAlignedTriangle
        point::Union{Nothing,Point}
        c::Union{Nothing,BaryTree}
        t::Union{Nothing,BaryTree}
        l::Union{Nothing,BaryTree}
        r::Union{Nothing,BaryTree}
    end

    function ∈(p::Point,t::AxisAlignedTriangle)::Bool
        """
            Test if a cartesian point is inside a triangle
                using BaryCentric Coordinates of the triangle
        """
        b = Cartesian2BaryCentric(p,t.b)
        if (b[1] >= 0 && b[2] >= 0 && b[1]+b[2] < 1)
            return true
        else
            return false
        end
    end

    function points(a::AxisAlignedTriangle)::Vector{Point}
        """
            Short hand for the corners of an AxisAlignedTriangle
        """
        return  Point.([
                a.b.r₁,
                a.b.r₂,
                a.b.r₃
            ])
    end

    function intersects(a::AxisAlignedTriangle,b::AxisAlignedTriangle)::Bool
        """
            Tests if any of the corners of triangle a are in triangle b
        """
        for p in points(a)
            if p ∈ b
                return true
            end
        end
        return false
    end

    ∉(p::Point,t::AxisAlignedTriangle)::Bool = !∈(p,t)

    function subdivide(t::AxisAlignedTriangle)::Vector{AxisAlignedTriangle}
        """
            Subdivides a triangle (as in the Sierpinski construction) and returns
            them each initialised with their own coordinate systems
        """
        m1,m2,m3 = Medians(t)
        return AxisAlignedTriangle.([
            BaryCentric(t.b.r₁,m3,m1),
            BaryCentric(m3,m2,m1),
            BaryCentric(m3,t.b.r₂,m2),
            BaryCentric(m1,m2,t.b.r₃)
        ])

    end

    ∈(p::Point,q::BaryTree)::Bool = ∈(p,q.boundary)
    ∉(p::Point,q::BaryTree)::Bool = ∉(p,q.boundary)
    isempty(q::BaryTree)::Bool = q.point == nothing

    BaryTree(a::AxisAlignedTriangle)::BaryTree = BaryTree(false,a,nothing,nothing,nothing,nothing,nothing)
    BaryTree(r::Bool,a::AxisAlignedTriangle)::BaryTree = BaryTree(r,a,nothing,nothing,nothing,nothing,nothing)

    function subdivide!(q::BaryTree)::Nothing
        """
            Subdivide this tree node as in the Sierpinski construction
            and store it in t,c,l,r child nodes
        """
        trigs = subdivide(q.boundary)
        q.t = BaryTree(trigs[1])
        q.c = BaryTree(trigs[2])
        q.l = BaryTree(trigs[3])
        q.r = BaryTree(trigs[4])
        nothing
    end

    function insert!(p::Point,q::BaryTree)::Bool
        """
            Trial an insert of point p (cartesian) in q's
            triangular boundary
        """
        if (p ∉ q)
            return false
        end

        if isempty(q) && q.t == nothing
            subdivide!(q)
            if p ∈ q.t
                q.t.point = p
            end

            if p ∈ q.c
                q.c.point = p
            end
            if p ∈ q.l
                q.l.point = p
            end
            if p ∈ q.r
                q.r.point = p
            end
            return true
        end

        if (q.t == nothing)
            subdivide!(q)
        end

        if insert!(p,q.t)
            return true
        end
        if insert!(p,q.c)
            return true
        end
        if insert!(p,q.l)
            return true
        end
        if insert!(p,q.r)
            return true
        end

        return false
    end

    Shape(a::AxisAlignedTriangle) = Shape(
        Shape([a.b.r₁[1],a.b.r₂[1],a.b.r₃[1]],[a.b.r₁[2],a.b.r₂[2],a.b.r₃[2]])
    )

    function Draw!(a::AxisAlignedTriangle)
        plot!(Shape(a),label="",fillalpha=0.0)
    end

    function Draw(Q::BaryTree,points=true)
        p = plot(aspect_ratio=:equal,label="")
        Draw!(Q.boundary)
        for q in [Q.t,Q.c,Q.l,Q.r]
            q != nothing ? Draw!(q,points) : nothing
        end
        return p
    end

    function Draw!(Q::BaryTree,points=true)
        Draw!(Q.boundary)
        for q in [Q.t,Q.c,Q.l,Q.r]
            q != nothing ? Draw!(q) : nothing
        end
    end

    function query(a::AxisAlignedTriangle, Q::BaryTree)
        if ~intersects(a,Q.boundary)
            return Vector{Point}([])
        end

        if Q.point == nothing
            return Vector{Point}([])
        end

        result = Vector{Point}()

        if (Q.point ∈ a)
            push!(result,Q.point)
        end

        for q in [Q.t,Q.c,Q.l,Q.r]
            if q != nothing
                for p in query(a,q)
                    push!(result,p)
                end
            end
        end

        return result
    end

    function size(Q::BaryTree)
        s = 1
        for q in [Q.t,Q.c,Q.l,Q.r]
            if q != nothing
                s += size(q)
            end
        end
        return s
    end
end
