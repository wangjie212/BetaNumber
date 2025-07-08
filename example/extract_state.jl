using JuMP
using COSMO
using LinearAlgebra

s = Matrix{Complex{Int8}}[[1 0; 0 1], [0 1; 1 0], [0 -im; im 0], [1 0; 0 -1]]
# a = [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [3, 0, 0, 0], [0, 3, 0, 0], [3, 3, 3, 0], [2, 3, 2, 1], [2, 2, 1, 1], [2, 1, 3, 3]]
a = [[1, 0, 0], [0, 1, 0], [3, 0, 0], [0, 3, 0], [2, 3, 1], [2, 1, 3], [3, 2, 3], [1, 2, 1], [2, 2, 0]]
q = 3
c = sqrt.(moment[1,2:10])
lower = zeros(2^9)
for j = 1:2^9
# println(j)
l = zeros(Int, 9)
p = j
for k = 1:9
    if p >= 2^k
        l[k] = 1
        p -= 2^k
    end
end
model = Model(optimizer_with_attributes(COSMO.Optimizer, "eps_abs" => 1e-8, "eps_rel" => 1e-8))
set_optimizer_attribute(model, MOI.Silent(), true)
pos = @variable(model, [1:2^q, 1:2^q] in HermitianPSDCone())
@constraint(model, tr(pos) == 1)
@objective(model, Min, sum((real(tr(kron(s[a[i][1]+1], s[a[i][2]+1], s[a[i][3]+1])*pos)) - (-1)^l[i]*c[i])^2 for i = 1:9))
optimize!(model)
objv = objective_value(model)
# println("optimum = $objv")
rho1 = value.(pos)
# println(sum([real(tr(kron(s[a[i][1]+1], s[a[i][2]+1], s[a[i][3]+1])*rho1))^2 for i = 1:9].*coe))
lower[j] = sum([real(tr(kron(s[a[i][1]+1], s[a[i][2]+1], s[a[i][3]+1])*rho1))^2 for i = 1:9].*coe)
end

io = open("D:/project/betanumber/rho27.txt", "w")
write(io, "$rho1")
close(io)

v = [0.5104979739, -0.5104979739im, 0.2113764731, -0.2113764731im, 0.1765626796, -0.1765626796im, 0.1224534324, -0.1224534324im, -0.0383245452im,
0.0383245452, 0.2114509722im, -0.2114509722, -0.2956187661im, 0.2956187661, -0.1223594828im, 0.1223594828]

rho = v*v'

b = [real(tr(kron(s[a[i][1]+1], s[a[i][2]+1], s[a[i][3]+1], s[a[i][4]+1])*rho)) for i = 1:9]
