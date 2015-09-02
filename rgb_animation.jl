# this import seemps to take a long time

for name in ["Images","ArgParse"]
    try
        require(name)
    catch
        println("Could not find package $(name), installing...")
        Pkg.add(name)
    end
end

using Images
using ArgParse

function main()
    println("starting...")
    parsed_args = parse_commandline()
    println("running with following arguments : ")
    for pa in parsed_args
        println("  $(pa[1])  =  $(pa[2])")
    end

    const OUT_DIR = parsed_args["output_dir"]
    const SIDE = parsed_args["side"]
    const FRAMES = parsed_args["frames"]
    const NUM_SEEDS = parsed_args["num_seeds"]
    const MAX_SAMPLES = parsed_args["max_samples"]
    const SORT_TYPE = parsed_args["sort_type"]
    const NEIGHBOUR_TYPE = parsed_args["neighbour_type"]

    if !ispath(OUT_DIR)
        mkdir(OUT_DIR)
    end

    const colors = get_color_list(SIDE, FRAMES, SORT_TYPE)
    
    println("setting up voxel cube")
    const a = [Voxel(convert(Int32,x),convert(Int32,y),convert(Int32,z)) for x in 1:SIDE, y in 1:SIDE, z in 1:FRAMES]

    start_time = time()

    #  we keep track of a list of available points (we want to 'randomly' sample, so a list-like 
    #  structure is needed here).  Instead of deleting available points when we use them, which
    #  be expensive, we keep track of those that we have deleted, and periodically refresh the
    #  available set list without the deleted points.  Finally, we keep track of those that we 
    #  have ever added to make sure we 
    available_points = Voxel[]
    deleted_points = Set{Voxel}()
    ever_seen = Set{Voxel}() 
    

    score_sum :: Int64 = 0

    for point in get_random_points(a, NUM_SEEDS)
        push!(available_points, point)
        push!(ever_seen, point)
    end

    iter :: Int32 = 0

    const big_val :: Int32 = 1073741824

    println("starting voxel placing...")

    for c in colors
        iter += 1
        
        # TODO : get rid of index, score stuff
        best_pt ,index, score = find_best_free_point(available_points, deleted_points, c, big_val, MAX_SAMPLES)
        score_sum += score

        set_color(best_pt, c)

        push!(deleted_points, best_pt)

        #if deleted point makes up more than 1% of the available points, we reset
        #the arrays
        if length(deleted_points) / length(available_points) > 0.01
            available_points = reset_array(available_points, deleted_points)
            deleted_points = Set{Voxel}()
        end

        add_colors_to_neighbours(NEIGHBOUR_TYPE, a, best_pt, c, available_points, ever_seen)

        if iter %10000 == 0
            elapsed_time = (time() - start_time) * 1.0 / 60
            estimated_time = (SIDE * SIDE * FRAMES - iter) * elapsed_time / iter

            println(iter,  ",  available points : " ,
                length(available_points), #- length(deleted_points),
                ", elapsed time : ", @sprintf("%0.3f", elapsed_time), " minutes,",
                " estimated remaining time : ", @sprintf("%0.3f", estimated_time), " minutes")
        end
    end
    
    make_images(a, OUT_DIR)

end

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "output_dir"
            help = "directory to save images when done (this will be created if it doesn't exist)"
            required = true
        "--side", "-s"
            help = "side length of square frames (must be power of 2)"
            arg_type = Int
            default = 128
        "--frames", "-f"
            help = "number of frames (must be power of 2)"
            arg_type = Int
            default = 128
        "--num_seeds"
            help = "number of random seed points"
            arg_type = Int
            default = 10
        "--max_samples"
            help = "maximum number of samples to take per iteration"
            arg_type = Int
            default = 1024
        "--sort_type"
            help = "how do we sort the colors prior to iterating, currently supported : rgb, random, hsv"
            arg_type = ASCIIString
            default = "random"
        "--neighbour_type"
            help = "which points do we count as neighbouring, the 6 closest (cross), or the 26 closest (cube)"
            arg_type = ASCIIString
            default = "cube"
    end

    return parse_args(s)
end

type Color_
    r :: Uint8
    g :: Uint8
    b :: Uint8
end

Color_() = Color_(-1,-1,-1)

# we store the color of the voxel, and the sum of the norms, and sum of 
# red, green, blue components.  This means calculating the distance to a color can
# be done in one pass, rather than iterating through all the neighbours
type Voxel
    x :: Int16
    y :: Int16
    z :: Int16
    color :: Color_
    num_neighbours :: Int16
    color_square_sum :: Int32
    red_sum :: Int32
    green_sum :: Int32
    blue_sum :: Int32
    is_coloured :: Bool
end   
    
Voxel(x :: Int32, y :: Int32, z :: Int32) = Voxel(x,y,z,Color_(), 0 ,0 ,0, 0 ,0, false)

function get_color_list(SIDE, FRAMES, SORT_TYPE)

    STEP = max(256 / SIDE, 1)
    all_colors = [c for c in [ Color_((x * STEP) % 256,(y * STEP)  % 256, (z* STEP)  % 256) 
                                            for x in 0 : SIDE - 1, y in 0 : SIDE - 1, z in 0: FRAMES - 1]]

    println("sorting list")
    if SORT_TYPE == "rgb"
        const colors = all_colors
    elseif SORT_TYPE == "random"
        const colors = shuffle(all_colors)
    elseif SORT_TYPE == "hue"
        const colors = sort(all_colors, by = get_hue)
    elseif SORT_TYPE == "brightness"
        const colors = sort(all_colors, by = get_brightness)
    else
        error("bad sort type argument : $(SORT_TYPE)")
    end
    colors
end

function get_random_points(a :: Array{Voxel}, n :: Int64)
    out = Set{Voxel}()
    X,Y,Z = size(a)
    for i in 1:n
        x = rand(1:X)
        y = rand(1:Y)
        z = rand(1:Z)
        push!(out, a[x, y, z])
    end
    out
end

# by default we don't check if the point is in the deleted pointed set, as this ends up being a bottleneck,
# instead we check to see if the best point is in the deleted set at the end, and if so, we rerun, being more
# careful that we don't take deleted points

function find_best_free_point(available_points, deleted_points, c,  big_val, MAX_SAMPLES)
    best_pt ,index, score = find_best_point(available_points, deleted_points, c,  big_val, MAX_SAMPLES, false)

    if best_pt.x == 0
        println("num available points : ", length(available_points))
        find_best_point(available_points, c, big_val, true)
        error("zero point?")
    end

    if in(best_pt, deleted_points)
        best_pt ,index, score = find_best_point(available_points, deleted_points, c,  big_val, MAX_SAMPLES, true)
    end
    best_pt ,index, score
end

function find_best_point(available_pts :: Array{Voxel,1}, deleted_points :: Set{Voxel}, 
        c :: Color_, big_val :: Int32, MAX_SAMPLES, check_deleted :: Bool)
    best_pt :: Voxel = Voxel(convert(Int32,0),convert(Int32,0),convert(Int32,0))
    best_score :: Int32 = big_val
    best_index = 0
    n = 1

    #  I found that randomly sampling the available points was too slow (may be that 
    #  I was using the wrong function though).  Instead, I choose an 'arithmetic progression'
    #  of the indices at random.  The list is getting shuffled semi regularly, so hopefully
    #  there isn't too much dependence between available points here

    step_size = max(itrunc(length(available_pts) / MAX_SAMPLES) , 1)
    index = rand(1 : step_size)

    while index <= length(available_pts) 
        vxl = available_pts[index]
        index += step_size
        if !check_deleted || !in(vxl, deleted_points)

            score = get_color_distance(vxl, c, true)
            if score <= best_score 
                if score < best_score
                    n = 1
                    best_score = score
                    best_pt = vxl
                    best_index = index
                
                #reservoir sampling, to make the order of available points less important
                else
                    if rand(1:n) == 1
                        n += 1
                        best_score = score
                        best_pt = vxl
                        best_index = index
                    end
                end
            end
        end
    end
    best_pt, best_index, best_score
end

function get_color_distance(vxl :: Voxel, c :: Color_, avg :: Bool)
    if vxl.num_neighbours == 0
        return -1073741824
    end
    out = vxl.color_square_sum - c.r * vxl.red_sum - c.g * vxl.green_sum - c.b * vxl.blue_sum

    # integer division for small speedup
    #div(out, vxl.num_neighbours)

    # uncomment next line if you want average assuming all the unfilled neighbours are black (gives odd effect)
    out
end

function set_color(vxl :: Voxel, c:: Color_)
    if vxl.is_coloured
        error("trying to colour voxel that has been coloured already")
    end
    vxl.color = c
    vxl.is_coloured = true
end

function reset_array(array, to_delete)
    out = Voxel[]
    for vxl in array
        if !in(vxl, to_delete)
            push!(out, vxl)
        end
    end
    shuffle(out)
end

function add_colors_to_neighbours(NEIGHBOUR_TYPE, a, best_pt, c, available_points, ever_seen)
    X,Y,Z = size(a)
    x_range = max(1, best_pt.x - 1): min(X, best_pt.x + 1)
    y_range = max(1, best_pt.y - 1): min(Y, best_pt.y + 1)
    z_range =  best_pt.z - 1 : best_pt.z + 1
    if NEIGHBOUR_TYPE == "cube"
        for x = x_range, y = y_range, z = z_range
            n = a[x, y, mod(z, Z)]
            color_neighbour_and_update_sets(n, c, available_points, ever_seen) 
        end

    elseif NEIGHBOUR_TYPE == "cross"
        #probably a nicer way of doing this...
        for x = x_range
            n = a[x, best_pt.y, mod(best_pt.z , Z)]
            color_neighbour_and_update_sets(n, c, available_points, ever_seen) 
        end

        for y = y_range
            n = a[best_pt.x, y, mod(best_pt.z , Z)]
            color_neighbour_and_update_sets(n, c, available_points, ever_seen) 
        end
        
        for z = z_range
            n = a[best_pt.x, best_pt.y, mod(z , Z)]
            color_neighbour_and_update_sets(n, c, available_points, ever_seen) 
        end
    else
        error("bad neighbour type argument : $(NEIGHBOUR_TYPE)")
    end
end

function color_neighbour_and_update_sets(n, c, available_points, ever_seen)
    if !n.is_coloured
        add_neighbour_color(n, c)
    end
    if !in(n, ever_seen)
        push!(available_points, n)
    end
    push!(ever_seen, n)
end

function add_neighbour_color(vxl :: Voxel, c:: Color_)
    vxl.num_neighbours += 1
    vxl.color_square_sum += c.r * c.r + c.g * c.g + c.b * c.b
    vxl.red_sum += 2 * c.r
    vxl.green_sum += 2 * c.g
    vxl.blue_sum += 2 * c.b
end

function mod(n, r)
    if n %r == 0
        return r
    else
        return (n + r) % r
    end
end

function make_images(a :: Array{Voxel,3}, out_dir)
    X,Y,Z = size(a)
    image_start_time = time()
    println("making images")
    for z in 1:Z
        vxl_array = a[1:X, 1:Y, z]
        im_array :: Array{Uint8,3} = [get_colour(vxl_array[x,y],i) for x = 1 : X, y = 1: Y, i = 1 : 3]
        Images.imwrite(Images.colorim(im_array),"$(out_dir)/test_$(dec(z, 4)).png")
    end

    println("time to make images : ", 
        @sprintf("%0.3f", (time() - image_start_time) * 1.0 / 60), " minutes.")
    println("finished!")
end

function get_colour(vxl :: Voxel, i )
    c = vxl.color
    if i == 1
        c.r
    elseif i == 2
        c.g
    elseif i == 3
        c.b
    else 
        error("bad color index")
    end
end

##  TODO : check this
function get_hue(c)
    r,g,b = c.r, c.g, c.b
    R = r / 256
    G = g / 256
    B = b / 256

    m = min(R,G,B)
    M = max(R,G,B)
    C = M - m
    
    if C == 0
        h = 0
    elseif M == R
        h = (G - B) / C
    elseif M == G
        h = 2 + (B - R) / C
    elseif M == B
        h = 4 + (R - G) / C
    end

    if h < 0
       h += 6
    end

    h * 60
end

function get_brightness(c)
    r,g,b = c.r, c.g, c.b
    r + g + b
end

main()

# for profiling code...
#Profile.init(10^7, 0.01)

#@profile main()

#Profile.print()
