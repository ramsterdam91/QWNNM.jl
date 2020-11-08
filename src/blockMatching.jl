function blockMatching(X::Array{Float32,3}, Par::PAR, Neighbor_arr::Array{Int32,2}, Num_arr::Array{Int32,1}, SelfIndex_arr::Array{Int32,1})::Array{Int32,2}

    L = length(Num_arr)

    Init_Index = zeros(Int32, Par.patnum[], L)

    for  i in 1 : L
        Patch = view(X, :, :, SelfIndex_arr[i])
        Nj = Num_arr[i]

        Dist = zeros(Float32,Nj)
        for j in 1:Nj
            @inbounds Neighbors = view(X,:,:,Neighbor_arr[j, i])
            Dist[j] = distance(Neighbors, Patch)
        end

        index = Int32.(sortperm(Dist, alg = QuickSort))

        Init_Index[:, i] .= Neighbor_arr[view(index,1:Par.patnum[]), i]
    end

    Init_Index
end

function distance(A::AbstractMatrix{T}, B::AbstractMatrix{T}) where T
    @assert size(A) == size(B) "size of matrices does not match, got $(size(A)) and $(size(B))"
    out = zero(T)
    for k=1:size(B, 2)
        for l=1:size(B, 1)
            @inbounds x = A[l,k] - B[l,k]
            out += abs2(x)
        end
    end
    return out
end
