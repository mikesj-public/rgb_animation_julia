import StatsBase
using Images

type Color_
    r :: Uint8
    g :: Uint8
    b :: Uint8
end

Color_() = Color_(-1,-1,-1)

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

function set_color(vxl :: Voxel, c:: Color_)
    if vxl.is_coloured
        error("trying to colour voxel that has been coloured already")
    end
    vxl.color = c
    vxl.is_coloured = true
end

function add_neighbour_color(vxl :: Voxel, c:: Color_)
    vxl.num_neighbours += 1
    vxl.color_square_sum += c.r * c.r + c.g * c.g + c.b * c.b
    vxl.red_sum += 2 * c.r
    vxl.green_sum += 2 * c.g
    vxl.blue_sum += 2 * c.b
end

const empty_pen = convert(Int32, 100)

function get_color_distance(vxl :: Voxel, c :: Color_, avg :: Bool)
    if vxl.num_neighbours == 0
        return 0
    end
    out = vxl.color_square_sum - c.r * vxl.red_sum - c.g * vxl.green_sum - c.b * vxl.blue_sum #- vxl.num_neighbours * empty_pen
    out
end

function get_random_points(a :: Array{Voxel}, n :: Int64)
    out = Set{Voxel}()
    X,Y,Z = size(a)
    for i in 1:n
        x = rand(1:X)
        y = rand(1:Y)
        z = rand(1:Z)
        push!(out, a[x, y, z])
        #println(size(out))
    end
    out
end

const MAX_SAMPLES = 5000

function find_best_point(available_pts :: Array{Voxel,1}, c :: Color_, big_val :: Int32)
    best_pt :: Voxel = Voxel(convert(Int32,0),convert(Int32,0),convert(Int32,0))
    best_score :: Int32 = big_val
    best_index = 0
    n = 1

    step_size = max(itrunc(length(available_pts) / MAX_SAMPLES) , 1)
    index = rand(1 : step_size)

    while index <= length(available_pts) 
        vxl = available_pts[index]
        score = get_color_distance(vxl, c, true)
        if score <= best_score 
            if score == best_score
                n = 1
                best_score = score
                best_pt = vxl
                best_index = index
            
            else
                if rand(1:n) == 1
                    n += 1
                    best_score = score
                    best_pt = vxl
                    best_index = index
                end
            end
        end
        index += step_size
    end
    best_pt, best_index, best_score
end

function find_best_point(available_pts :: Set{Voxel}, c :: Color_, big_val :: Int32)
    best_pt :: Voxel = Voxel(convert(Int32,0),convert(Int32,0),convert(Int32,0))
    best_score :: Int32 = big_val
    index = 1
    count = 1
    n = 1

    for vxl in available_pts
        score :: Int32 = get_color_distance(vxl, c, true)
        if score <= best_score 
            if score == best_score
                n = 1
                best_score = score
                best_pt = vxl
                index = count
            
            else
                if rand(1:n) == 1
                    n += 1
                    best_score = score
                    best_pt = vxl
                    index = count
                end
            end
        end
        count += 1
        
    end
    best_pt, index, best_score
end

function mod(n :: Int64, r :: Int64)
    if n %r == 0
        return r
    else
        return (n + r) % r
    end
end


function update_available(available_points, deleted_points)
    new_avail = Set{Voxel}()
    for vxl in available_points
        if !in(vxl, deleted_points)
            push!(new_avail, vxl)
        end
    end
    println(length(available_points) - length(deleted_points), "_", length(new_avail))
    new_avail
end

function get_colour(vxl :: Voxel)
    c = vxl.color
    [c.r, c.g, c.b]
end

function make_images(a :: Array{Voxel,3})
    X,Y,Z = size(a)
    for z in 1:Z
        im_array :: Array{Uint8,3} = [get_colour(a[x,y,z])[i] for x = 1 : X, y = 1: Y, i = 1 : 3]
        Images.imwrite(Images.colorim(im_array),"images/test_$(dec(z, 4)).png")
    end
end

function insert_into(array :: Array{Int64, 1}, el :: Int64 ,  index :: Int64)
    push!(array, array[index])
    array[index] = el
end

function insert_into(array, el)
    index = rand(1: length(array) + 1)
    head = array[1 : index - 1]
    tail = array[index : length(array)]
    vcat(head, el, tail)
end

function run()

    # const SIDE = 128
    # const STEP = 256 / SIDE

    # const all_colors = shuffle([c for c in [ Color_(x * STEP % 256, y * STEP  % 256, z* STEP  % 256) for x in 0 : SIDE - 1, 
    #                                                         y in 0 : SIDE - 1, z in 0:SIDE - 1]])

    const SIDE = 512
    const FRAMES = 256
    const STEP = 1

    const all_colors = shuffle([c for c in [ Color_((x * STEP) % 256, (y * STEP)  % 256, (z* STEP)  % 256) for x in 0 : SIDE - 1, 
                                                            y in 0 : SIDE - 1, z in 0: FRAMES - 1]])

    a = [Voxel(convert(Int32,x),convert(Int32,y),convert(Int32,z)) for x in 1:SIDE, y in 1:SIDE, z in 1:FRAMES]

    start_time = time()

    #available_points = Set{Voxel}()
    available_points = Voxel[]
    ever_seen = Set{Voxel}() 
    #deleted_points = Set{Voxel}()

    score_sum :: Int64 = 0

    for point in get_random_points(a, 2)
        push!(available_points, point)
        push!(ever_seen, point)
    end

    iter :: Int32 = 0
    const X = size(a)[1]
    const Y = size(a)[2]
    const Z = size(a)[3]

    const big_val :: Int32 = 1073741824

    for c in all_colors
        iter += 1
        best_pt ,index, score = find_best_point(available_points, c, big_val)
        score_sum += score

        set_color(best_pt, c)

        if best_pt.x == 0
            println("num available points : ", length(available_points))
            error("zero point?")

        end
        
        x_range = max(1, best_pt.x - 1): min(X, best_pt.x + 1)
        y_range = max(1, best_pt.y - 1): min(Y, best_pt.y + 1)
        z_range =  best_pt.z - 1 : best_pt.z + 1

        splice!(available_points, index)
        for x = x_range, y = y_range, z = z_range
            n = a[x,y,mod(z , Z)]
            if !n.is_coloured
                add_neighbour_color(n, c)
            end
            if !in(n, ever_seen)
                # available_points = insert_into(available_points, n)  
                push!(available_points, n)
            end
            push!(ever_seen, n)
        end
        
        
        if iter %10000 == 0
            elapsed_time = (time() - start_time) * 1.0 / 60
            estimated_time = (SIDE * SIDE * FRAMES - iter) * elapsed_time / iter

            println(iter, "  " , "avg  score : ", score_sum * 1. / iter, "  " ,
                length(available_points), #- length(deleted_points),
                " elapsed time : ", @sprintf("%0.3f", elapsed_time), " minutes",
                " estimated remaining time : ", @sprintf("%0.3f", estimated_time), " minutes")
        end

    end

    make_images(a)
end

run()
#@profile run()

#Profile.print()
