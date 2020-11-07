function blockMatching(X::Array{Float32,2}, Par::PAR, Neighbor_arr::Array{Int32,2}, Num_arr::Array{Int32,2}, SelfIndex_arr::Array{Int32,2})::Array{Int32,2}

    L = length(Num_arr)

    Init_Index = zeros(Int32, Par.patnum, L)

    for  i in 1 : L
        Patch = X[:, SelfIndex_arr[i], :]

        Neighbors = X[:, Neighbor_arr[1:Num_arr[i], i], :]
        
        Dist = zeros(Float32,size(Neighbors)[2])
        for j in 1:size(Neighbors)[2]
            Dist[j] = sum(float32.(Neighbors[:,j,:] - Patch).^2)
        end
      
        index = Int32.(sortperm(Dist, alg = QuickSort))

        Init_Index[:, i] = Neighbor_arr[index[1:Par.patnum], i] 
    end

    Init_Index
end