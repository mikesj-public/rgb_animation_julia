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
            default = "cross"
        "--smoothing_grid_side"
            help = "how many cells to split search field into (must divide 256)"
            arg_type = Int
            default = 8
    end

    return parse_args(s)
end


function main_run()
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
    const SMOOTH_GRID_SIDE = parsed_args["smoothing_grid_side"]

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
    
    #main delete array
    mda = Set{Voxel}()
    ever_seen = Set{Voxel}()

    println("generating lookup grid")

    grid = [Set{Voxel}() for x in 1:SMOOTH_GRID_SIDE, y = 1:SMOOTH_GRID_SIDE, z =1:SMOOTH_GRID_SIDE]
    
    score_sum :: Int64 = 0

    for point in get_random_points(a, NUM_SEEDS)
        push!(mda, point)
        add_to_grid(grid, point, SMOOTH_GRID_SIDE)
        push!(ever_seen, point)
    end

    iter :: Int32 = 0

    const big_val :: Int32 = 1073741824

    println("starting voxel placing...")

    for c in colors
        iter += 1
        
        # TODO : get rid of index, score stuff
        #best_pt ,index, score, fail, big_fail = find_best_free_point(mda, grid,  c, big_val, MAX_SAMPLES)
        best_pt = find_best_free_point(mda, grid,  c, big_val, MAX_SAMPLES, SMOOTH_GRID_SIDE)

        set_color(best_pt, c)

        delete!(mda, best_pt)

        delete_from_grid(grid, best_pt, SMOOTH_GRID_SIDE)

        add_colors_to_neighbours(NEIGHBOUR_TYPE, a, best_pt, c, mda, ever_seen, grid, SMOOTH_GRID_SIDE)

        if iter %10000 == 0
            elapsed_time = (time() - start_time) * 1.0 / 60
            estimated_time = (SIDE * SIDE * FRAMES - iter) * elapsed_time / iter

            println(iter,  ",  available points : " ,
                length(mda),
                ", elapsed time : ", @sprintf("%0.3f", elapsed_time), " minutes,",
                " estimated remaining time : ", @sprintf("%0.3f", estimated_time), " minutes")
                #", failures : $(failures), big failures : $(big_failures), resets $(resets)")
        end
        if iter > 5000000
            break
        end
    end
    
    #make_images(a, OUT_DIR)

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

# type DeleteArray
#     available_points :: Array{Voxel,1}
#     available_size :: Int
#     deleted_points :: Set{Voxel} 
# end

# DeleteArray() = DeleteArray(Voxel[], 0, Set{Voxel}())

# function delete_array_delete(delete_array ::DeleteArray, pt :: Voxel)
#     reset = 0
#     push!(delete_array.deleted_points, pt)
#     if length(delete_array.deleted_points) / delete_array.available_size > 0.01 && delete_array.available_size > 2  && length(delete_array.deleted_points) > 2
#         reset_delete_array(delete_array)
#         reset = 1
#     end
#     reset
# end

# function reset_delete_array(da :: DeleteArray)
#     to_delete :: Set{Voxel} = da.deleted_points
#     pts = da.available_points
#     index = 1
#     for i = 1:da.available_size
#         next = pts[i]
#         if !in(next, to_delete)
#             pts[index] = next
#             index += 1
#         end
#     end
#     da.available_size = index - 1

#     da.deleted_points = Set{Voxel}()
    
#     if length(da.deleted_points) > 0
#         error("deleted size error")
#     end
#     if da.available_size > length(da.available_points)
#         error("available size error")
#     end
# end

# function delete_array_add(da :: DeleteArray, pt ::Voxel)
#     if in(pt, da.deleted_points)
#         error("adding deleted point")
#     end
#     if length(da.available_points) > da.available_size
#         da.available_points[da.available_size + 1] = pt
#     else
#         push!(da.available_points, pt)
#     end
#     da.available_size += 1
#     if da.available_size > length(da.available_points)
#         error("size error")
#     end

# end

# function delete_array_length(da :: DeleteArray)
#     da.available_size - length(da.deleted_points)
# end

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

function find_best_free_point(da :: Set{Voxel}, grid :: Array{Set{Voxel},3}, c,  big_val, MAX_SAMPLES, SMOOTH_GRID_SIDE)
    
    fail = 0
    big_fail = 0

    best_vxl_array = get_best_vxl_set(grid, c, SMOOTH_GRID_SIDE)
    if length(best_vxl_array) > 8
        array = best_vxl_array
        samples = 128
    else
        array = da
        samples = MAX_SAMPLES
    end

    best_pt = find_best_point(array, c,  big_val, samples)

    if best_pt.x == 0
        best_pt  = find_best_point(da, c,  big_val, MAX_SAMPLES)
        big_fail = 1
    end

    if best_pt.x == 0
        #samples = length(array.available_points)
        println("num available points : ", delete_array_length(best_vxl_array),
            ", ", delete_array_length(da) )

        println("check one " , length(best_vxl_array.available_points) - length(best_vxl_array.deleted_points))
        println("check two " , length(da.available_points) - length(da.deleted_points))
        find_best_point(da, c,  big_val, MAX_SAMPLES)

        error("zero point?")
    end

    best_pt
end

function get_best_vxl_set(grid ::Array{Set{Voxel},3}, c :: Color_, SMOOTH_GRID_SIDE)
    side_length = 256 / SMOOTH_GRID_SIDE
    r = div(c.r, side_length) + 1
    g = div(c.g, side_length) + 1
    b = div(c.b, side_length) + 1
    grid[r,g,b]
end

function get_best_vxl_set(grid ::Array{Set{Voxel},3}, r,g,b, SMOOTH_GRID_SIDE)
    side_length = 256 / SMOOTH_GRID_SIDE
    r2 = div(r, side_length) + 1
    g2 = div(g, side_length) + 1
    b2 = div(b, side_length) + 1
    grid[r2,g2,b2]
end

const empty_vxl = Voxel(convert(Int32,0),convert(Int32,0),convert(Int32,0))

function find_best_point(array :: Set{Voxel}, 
        c :: Color_, big_val :: Int32, MAX_SAMPLES)

    best_pt :: Voxel = empty_vxl
    
    best_score :: Int32 = big_val
    best_index = 0

    #  I found that randomly sampling the available points was too slow (may be that 
    #  I was using the wrong function though).  Instead, I iterate through an arithmetic progression
    #  of the indices at random.  The list is getting shuffled semi regularly, so hopefully
    #  there isn't too much dependence between available points here  
    counter = 0

    for vxl in array
        score = get_color_distance(vxl, c, true)

        if score <= best_score 
            if score < best_score
                best_score = score
                best_pt = vxl
            end
        end

        counter += 1
        if counter == MAX_SAMPLES
            break
        end
    end
    best_pt
end

function get_color_distance(vxl :: Voxel, c :: Color_, avg :: Bool)
    neibs = vxl.num_neighbours
    if vxl.num_neighbours == 0
        return -1073741824
    end
    out = vxl.color_square_sum + c.r * (c.r * neibs - vxl.red_sum)  + c.g * (c.g * neibs - vxl.green_sum) + (c.b * neibs - c.b * vxl.blue_sum)

    # integer division for small speedup
    #div(out, neibs )

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

function delete_from_grid(grid, pt, SMOOTH_GRID_SIDE)
    nn = pt.num_neighbours
    if nn == 0
        r,g,b = 128,128,128
    else
        r = div(pt.red_sum , 2 * nn)
        g = div(pt.green_sum , 2 * nn)
        b = div(pt.blue_sum , 2 * nn)
    end
    set :: Set{Voxel} = get_best_vxl_set(grid,r,g,b, SMOOTH_GRID_SIDE)
    delete!(set, pt)
end

function add_to_grid(grid, pt, SMOOTH_GRID_SIDE)
    nn = pt.num_neighbours
    if nn == 0
        r,g,b = 128,128,128
    else
        r = div(pt.red_sum , 2 * nn)
        g = div(pt.green_sum , 2 * nn)
        b = div(pt.blue_sum , 2 * nn)
    end
    set :: Set{Voxel} = get_best_vxl_set(grid,r,g,b, SMOOTH_GRID_SIDE)
    push!(set, pt)
end

function add_colors_to_neighbours(NEIGHBOUR_TYPE, a :: Array{Voxel,3}, best_pt :: Voxel, c :: Color_, da :: Set{Voxel}, ever_seen, grid, SMOOTH_GRID_SIDE)
    X,Y,Z = size(a)
    x_range = max(1, best_pt.x - 1): min(X, best_pt.x + 1)
    y_range = max(1, best_pt.y - 1): min(Y, best_pt.y + 1)
    z_range =  best_pt.z - 1 : best_pt.z + 1


    if NEIGHBOUR_TYPE == "cube"
        for x = x_range, y = y_range, z = z_range
            pt = a[x, y, mod(z, Z)]
            color_neighbour_and_update_sets(pt, c, da, ever_seen, grid, SMOOTH_GRID_SIDE) 
        end

    elseif NEIGHBOUR_TYPE == "cross"
        #probably a nicer way of doing this...
        for x = x_range
            pt = a[x, best_pt.y, mod(best_pt.z , Z)]
            color_neighbour_and_update_sets(pt, c, da, ever_seen, grid, SMOOTH_GRID_SIDE) 
        end

        for y = y_range
            pt = a[best_pt.x, y, mod(best_pt.z , Z)]
            color_neighbour_and_update_sets(pt, c, da, ever_seen, grid, SMOOTH_GRID_SIDE) 
        end
        
        for z = z_range
            pt = a[best_pt.x, best_pt.y, mod(z , Z)]
            color_neighbour_and_update_sets(pt, c, da, ever_seen, grid, SMOOTH_GRID_SIDE) 
        end
    else
        error("bad neighbour type argument : $(NEIGHBOUR_TYPE)")
    end

end

function color_neighbour_and_update_sets(pt :: Voxel, c :: Color_, da :: Set{Voxel}, ever_seen, grid, SMOOTH_GRID_SIDE)

    if !pt.is_coloured
        nn = pt.num_neighbours
        if nn == 0
            r1,g1,b1 = 128,128,128
        else
            r1 = div(pt.red_sum , 2 * nn)
            g1 = div(pt.green_sum , 2 * nn)
            b1 = div(pt.blue_sum , 2 * nn)
        end
        set1 = get_best_vxl_set(grid,r1,g1,b1, SMOOTH_GRID_SIDE)
        add_neighbour_color(pt, c)
        r2 = div(pt.red_sum , 2 * pt.num_neighbours)
        g2 = div(pt.green_sum , 2 * pt.num_neighbours)
        b2 = div(pt.blue_sum , 2 * pt.num_neighbours)
        set2 = get_best_vxl_set(grid,r2,g2,b2, SMOOTH_GRID_SIDE)
        if set1 != set2
            delete!(set1, pt)
            push!(set2, pt)
        end
    end

    if !in(pt, ever_seen)
        push!(da, pt)
    end
    push!(ever_seen, pt)
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

#main_run()

# for profiling code...
Profile.init(10^7, 0.1)


@profile main_run()

Profile.print()