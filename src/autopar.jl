#Generate Parameters of WNNM-algorithm
#input the nSig  the noise Level

mutable struct PAR{T<:Int64, F<:Float32}
    nSig        ::T                         #level of Gaussian noise
    SearchWin   ::T                         #the search range for Neighbour Patches
    delta       ::F                         #parameter of iterative method
    c           ::F              
    Innerloop   ::T

    patsize     ::T                         #the size of Patches is: patsize * patsize
    patnum      ::T                   
    Iter        ::T                         #parameter of iterative method
    lamada      ::F
    step        ::T
end

function autopar(nSig)
  
    SearchWin =   30  
    delta     =   0.1
    c         =   3*√(2)
    Innerloop =   2

    if          nSig <= 20
        patsize       =   6             
        patnum        =   70            
        Iter          =   8            
        lamada        =   0.54 
    elseif      nSig <= 40
        patsize       =   7             
        patnum        =   90           
        Iter          =   12            
        lamada        =   0.56  
    elseif      nSig <= 60
        patsize       =   8             
        patnum        =   120          
        Iter          =   14            
        lamada        =   0.58    
    else
        patsize       =   9             
        patnum        =   140           
        Iter          =   14          
        lamada        =   0.58   
    end
    step              =   patsize ÷ 2

    PAR(nSig, SearchWin, delta, c, Innerloop, patsize, patnum , Iter, lamada, step)
end