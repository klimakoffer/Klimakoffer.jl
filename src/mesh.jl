

struct Mesh
  coordinates::Array{Float64, 3}
end


function Base.show(io::IO, mesh::Mesh)
  print(io, "Mesh() with ", size(mesh.coordinates, 2), "×", size(mesh.coordinates, 3), " points")
end
