include("autopar.jl")
#input the Original Image `Img_O`, Image with noise  `Img_N`, sets of parameters `Par`
#use Iterative method to get Estimated Image `Img_E`
#output the Estimated Image `Img_E` in the form of Array{Float32,3} H x W x 3 (3 = RGB Channels)

include("blockMatching.jl") #SLOWWWWWWW
include("patEstimation.jl") #SLOWWWW because of ↓↓↓↓↓↓
include("WNNM.jl")          #CORE and its SLOWWWWWWWW

function denoising(Img_O::Array{Float32,3}, Img_N::Array{Float32,3}, Par::PAR)
    #use Iterative method to get Estimated Image
    Img_E = Img_N

    (H,W,D) = size(Img_N)

    #the total number of patches
    TotalPatNum = (H - Par.patsize + 1) * (W - Par.patsize + 1)

    (Neighbor_arr, Num_arr, Self_arr) = neighborIndex(H, W, Par)
    NL_mat              =   zeros(Int32,   Par.patnum[], length(Num_arr))

    #PSNR reflect the quality of `Img_E`
    PSNR_arr            =   zeros(Float64, 1, Par.Iter)

    for iter in 1 : Par.Iter

        #Iterative method
        Img_E =	Img_E + Par.delta * (Img_N - Img_E)

        #image to patch `CurPat` and estimate local noise variance `Sigma_arr`
        (CurPat, Sigma_arr)	= im2patch(Img_E, Img_N, Par)

        if (mod(iter - 1, Par.Innerloop) == 0)

            #Lower noise level, less None-local patches
            Par.patnum[] = Par.patnum[] - 10

            #Caculate Non-local similar patches for each
            NL_mat  =  blockMatching(CurPat, Par, Neighbor_arr, Num_arr, Self_arr)

            #First Iteration use the input noise parameter
            if(iter == 1)
                Sigma_arr = ones(Float32, size(Sigma_arr)) .* Par.nSig
            end
        end
        #Estimate all the patches
        (EPat, WPat) = patEstimation(NL_mat, Self_arr, Sigma_arr, permutedims(CurPat, (2,3,1)), Par)

        #get Estimated image from patches
        Img_E = patch2Im(EPat, WPat, Par.patsize, H, W, D)

        #record PSNR in each Iteration
        PSNR_arr[iter]  = psnr(Img_O, Img_E)
        println("in Iteration:", iter, " PSNR = ", PSNR_arr[iter])

        #saturation situation
        if (iter > 1) && (PSNR_arr[iter] - PSNR_arr[iter-1] < 0.01)
            break
        end

    end
    Img_E
end

function neighborIndex(H, W, Par::PAR)
    SW          =   Par.SearchWin
    s           =   Par.step
    patsr       =   H - Par.patsize + 1
    patsc       =   W - Par.patsize + 1

    #Generate the GridIndex of patches
    r_Grid	    =   Vector(1 : s : patsr)::Array{Int64,1}
    r_Grid	    =   [r_Grid; (last(r_Grid) + 1):patsr]
    c_Grid	    =   Vector(1 : s : patsc)::Array{Int64,1}
    c_Grid	    =   [c_Grid; (last(c_Grid) + 1):patsc]

    #index of patches ( col-first↓ )
    IdxV        =   Vector(1 : patsr * patsc):: Array{Int64,1}
    Idx         =   reshape(IdxV, patsr, patsc):: Array{Int64,2}

    r_Gridn     =   length(r_Grid)
    c_Gridm     =   length(c_Grid)

    Neighbor_arr    =   zeros(Int32, (2 * SW + 1)^2, r_Gridn * c_Gridm)
    Num_arr         =   zeros(Int32,  r_Gridn * c_Gridm)
    SelfIndex_arr   =   zeros(Int32,  r_Gridn * c_Gridm)

    for i in 1 : r_Gridn
        for j in 1 : c_Gridm
            #offsrt r and c of pats
            OffsetR     =   r_Grid[i] :: Int64
            OffsetC     =   c_Grid[j] :: Int64

            #Offset1 is the index of pats(the left-up conor of grid)
            Offsetpats   	=  (OffsetC - 1) * patsr + OffsetR :: Int64

            #Offset2 is the index of grid
            Offsetgrid   	=  (j - 1) * r_Gridn + i :: Int64

            #make sure the windows is in the valid region
            top         =   max(OffsetR - SW, 1)
            button      =   min(OffsetR + SW, patsr)
            left        =   max(OffsetC - SW, 1)
            right       =   min(OffsetC + SW, patsc)

            NL_IdxV     =   Idx[top:button, left:right] :: Array{Int64,2}
            NL_Idx     =   NL_IdxV[:]                   :: Array{Int64,1}

            Num_arr[Offsetgrid]  =  length(NL_Idx)

            Neighbor_arr[1:Num_arr[Offsetgrid], Offsetgrid]  =  NL_Idx[1:Num_arr[Offsetgrid]]
            SelfIndex_arr[Offsetgrid] = Offsetpats
        end
    end
    (Neighbor_arr, Num_arr, SelfIndex_arr)
end

function im2patch(Img_E::Array{Float32,3}, Img_N::Array{Float32,3}, Par::PAR)

    (H, W, D) = size(Img_E)
    N           =   H - Par.patsize + 1
    M           =   W - Par.patsize + 1
    TotalPatNum =   N * M

    Y           =   zeros(Float32, D, Par.patsize^2, TotalPatNum)       #Current Patches
    N_Y         =   zeros(Float32, D, Par.patsize^2, TotalPatNum)       #Noisy   Patches

    k = 0
    for i in 1 : Par.patsize
        for j in 1 : Par.patsize
            k += 1

            n = H - Par.patsize + i
            m = W - Par.patsize + j

            E_patch      =  Img_E[i:n, j:m, :]::Array{Float32,3}
            N_patch      =  Img_N[i:n, j:m, :]::Array{Float32,3}

            Y[:,k,:]     .=  reshape(E_patch, (N*M, D))'
            N_Y[:,k,:]   .=  reshape(N_patch, (N*M, D))'

        end
    end

    #Estimated Local Noise Level
    SigmaArr = Par.lamada .* sqrt.(abs.(dropdims(mean((N_Y .- Y).^2, dims=2), dims=2)' .- Par.nSig^2 ))

    (Y, SigmaArr)
end

function patch2Im(ImPat::Array{Float32,3}, WPat::Array{Int32,3}, PatSize, H, W, D)
    N       =   H - PatSize + 1
    M       =   W - PatSize + 1
    TempOffsetR  =   Vector(1:N)
    TempOffsetC  =   Vector(1:M)

    Img_E  	=  zeros(Float32,H, W, D)
    Img_W 	=  zeros(Float32,H, W, D)

    k = 0
    for i in 1 : PatSize
        for j in 1 : PatSize
            k += 1

            Img_E[TempOffsetR .+ (i-1), TempOffsetC .+ (j-1), :]  += reshape( ImPat[k,:,:], (N, M, D))
            Img_W[TempOffsetR .+ (i-1), TempOffsetC .+ (j-1), :]  += reshape( WPat[k,:,:],  (N, M, D))

        end
    end
    Img_E   =  Float32.(Img_E ./ (Img_W .+ eps()))
end
