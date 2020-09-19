using Javis

function ground(args...) 
    background("white") # canvas background
    sethue("black") # pen color
end

function object(;p=O, color="black")
    sethue(color)
    s = ngon(p, 25, 4)
    poly(s, :fill, close=true)
    return p
end

myvideo = Video(500, 500)

actions = [Action((args...)-> object(;color="green"), 
           Translation(Point(-250, 0), Point(250, 0))),
]

javis(
    myvideo,  
    [
        BackgroundAction(1:500, ground), actions...
    ],
    pathname="particle.gif"
)
