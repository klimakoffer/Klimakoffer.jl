
struct Discretization{LType, UType, PType}
  L::LType
  U::UType
  p::PType
  NT::Int64
  mesh::Mesh
  model::Model
end


function Discretization(mesh, model, NT)
    A = ComputeMatrix(mesh,NT,model)
    LUdec = lu(A)
    L = sparse(LUdec.L)
    U = sparse(LUdec.U)

    Discretization(L, U, LUdec.p, NT, mesh, model)
end
