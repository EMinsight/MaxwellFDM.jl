using MaxwellFDM
using BenchmarkTools
using StaticArrays
using JLD2

xprim = [
-400.0000
-390.0000
-380.0000
-370.0000
-360.0000
-350.0000
-340.0000
-330.0000
-320.0000
-310.0000
-300.0000
-290.0000
-279.4727
-268.9455
-258.4182
-247.8910
-237.3637
-226.8365
-216.3092
-205.7820
-195.2547
-184.7275
-174.2002
-163.6730
-153.1457
-142.6185
-132.0912
-121.5640
-111.5640
-105.9406
-102.7783
-101.0000
-100.0000
 -99.0000
 -98.0000
 -97.0000
 -96.0000
 -95.0000
 -94.0000
 -93.0000
 -92.0000
 -91.0000
 -90.0000
 -89.0000
 -88.0000
 -87.0000
 -86.0000
 -85.0000
 -84.0000
 -83.0000
 -82.0000
 -81.0000
 -80.0000
 -79.0000
 -78.0000
 -77.0000
 -76.0000
 -75.0000
 -74.0000
 -73.0000
 -72.0000
 -71.0000
 -70.0000
 -69.0000
 -68.0000
 -67.0000
 -66.0000
 -65.0000
 -64.0000
 -63.0000
 -62.0000
 -61.0000
 -60.0000
 -59.0000
 -58.0000
 -57.0000
 -56.0000
 -55.0000
 -54.0000
 -53.0000
 -52.0000
 -51.0000
 -50.0000
 -49.0000
 -48.0000
 -47.0000
 -46.0000
 -45.0000
 -44.0000
 -43.0000
 -42.0000
 -41.0000
 -40.0000
 -39.0000
 -38.0000
 -37.0000
 -36.0000
 -35.0000
 -34.0000
 -33.0000
 -32.0000
 -31.0000
 -30.0000
 -29.0000
 -28.0000
 -27.0000
 -26.0000
 -25.0000
 -24.0000
 -23.0000
 -22.0000
 -21.0000
 -20.0000
 -19.0000
 -18.0000
 -17.0000
 -16.0000
 -15.0000
 -14.0000
 -13.0000
 -12.0000
 -11.0000
 -10.0000
  -9.0000
  -8.0000
  -7.0000
  -6.0000
  -5.0000
  -4.0000
  -3.0000
  -2.0000
  -1.0000
        0
   1.0000
   2.0000
   3.0000
   4.0000
   5.0000
   6.0000
   7.0000
   8.0000
   9.0000
  10.0000
  11.0000
  12.0000
  13.0000
  14.0000
  15.0000
  16.0000
  17.0000
  18.0000
  19.0000
  20.0000
  21.0000
  22.0000
  23.0000
  24.0000
  25.0000
  26.0000
  27.0000
  28.0000
  29.0000
  30.0000
  31.0000
  32.0000
  33.0000
  34.0000
  35.0000
  36.0000
  37.0000
  38.0000
  39.0000
  40.0000
  41.0000
  42.0000
  43.0000
  44.0000
  45.0000
  46.0000
  47.0000
  48.0000
  49.0000
  50.0000
  51.0000
  52.0000
  53.0000
  54.0000
  55.0000
  56.0000
  57.0000
  58.0000
  59.0000
  60.0000
  61.0000
  62.0000
  63.0000
  64.0000
  65.0000
  66.0000
  67.0000
  68.0000
  69.0000
  70.0000
  71.0000
  72.0000
  73.0000
  74.0000
  75.0000
  76.0000
  77.0000
  78.0000
  79.0000
  80.0000
  81.0000
  82.0000
  83.0000
  84.0000
  85.0000
  86.0000
  87.0000
  88.0000
  89.0000
  90.0000
  91.0000
  92.0000
  93.0000
  94.0000
  95.0000
  96.0000
  97.0000
  98.0000
  99.0000
 100.0000
 101.0000
 102.7783
 105.9406
 111.5640
 121.5640
 132.0912
 142.6185
 153.1457
 163.6730
 174.2002
 184.7275
 195.2547
 205.7820
 216.3092
 226.8365
 237.3637
 247.8910
 258.4182
 268.9455
 279.4727
 290.0000
 300.0000
 310.0000
 320.0000
 330.0000
 340.0000
 350.0000
 360.0000
 370.0000
 380.0000
 390.0000
 400.0000
]
yprim = copy(xprim)
zprim = copy(xprim)

M = (length(xprim)-1, length(yprim)-1, length(zprim)-1)
println("M = $M")

lprim = (xprim, yprim, zprim)

isbloch = [true, true, true]

g3 = Grid(lprim, isbloch)

# Create materials.
εvac = 1.0
vac = Material("Vacuum", ε=εvac)

εdiel = 12.0
diel = Material("Dielectric", ε=εdiel)

# Create objects.
dom_vac = Object(Box(g3.bounds), vac)
obj_diel = Object(Sphere([0,0,0], 75), diel)
# obj_diel = Object(Box([0,0,0], [75,75,75]), diel)
obj_xn_diel = Object(Sphere([-150,0,0], 75), diel)
obj_xp_diel = Object(Sphere([150,0,0], 75), diel)
obj_yn_diel = Object(Sphere([0,-150,0], 75), diel)
obj_yp_diel = Object(Sphere([0,150,0], 75), diel)
obj_zn_diel = Object(Sphere([0,0,-150], 75), diel)
obj_zp_diel = Object(Sphere([0,0,150], 75), diel)

# Add objects.
ovec = Object3[]
paramset = (SMat3Complex[], SMat3Complex[])
# add!(ovec, paramset, dom_vac)
# add!(ovec, paramset, dom_vac, obj_diel)
add!(ovec, paramset, dom_vac, obj_diel, obj_xn_diel, obj_xp_diel, obj_yn_diel, obj_yp_diel, obj_zn_diel, obj_zp_diel)

N = g3.N
ε3d = create_param3d(N)
εobj3d = create_n3d(Object3, N)
εind3d = create_n3d(ParamInd, N)
εoind3d = create_n3d(ObjInd, N)

μ3d = create_param3d(N)
μobj3d = create_n3d(Object3, N)
μind3d = create_n3d(ParamInd, N)
μoind3d = create_n3d(ObjInd, N)

τl = g3.ghosted.τl

gt = PRIM
nw = 1
ngt = Int(gt)
ngt′ = alter(ngt)

M = length.(τl[nPR])  # N+1
ind = (map((m,b)->circshift(1:m, b), M, g3.isbloch.data),
       map((m,b)->circshift(1:m, -b), M, g3.isbloch.data))

gt_cmp = SVector(gt, gt, gt)
gt_cmp = map((k,g)->(k==nw ? alter(g) : g), nXYZ, gt_cmp)  # no change if nw = 0

# Choose the circularly shifted indices to use.
ind_cmp = MaxwellFDM.t_ind(ind, gt_cmp)

# Prepare the circularly shifted locations of the field components.
τlcmp = view.(MaxwellFDM.t_ind(τl,gt_cmp), ind_cmp)  # Tuple3{VecFloat}: locations of Fw = Uw or Vw

# Prepare the circularly shifted viewes of various arrays to match the sorted
# τl.  Even though all arrays are for same locations, param3d_cmp contains gt
# material, whereas obj3d_cmp, pind3d_cmp, oind3d_cmp contain alter(gt)
# material, so use ngt′ instead of ngt for them.
ε3d_cmp = view(ε3d, ind_cmp..., nXYZ, nXYZ)
μobj3d_cmp = view(μobj3d[nw], ind_cmp...)
μind3d_cmp = view(μind3d[nw], ind_cmp...)
μoind3d_cmp = view(μoind3d[nw], ind_cmp...)

# o = ovec[2]  # ovec[1]: Box, ovec[2]: Sphere
# shape = o.shape
#
# gt′ = alter(gt)
# param = matparam(o,gt)
# pind′ = paramind(o,gt′)
# oind = objind(o)
#
# arrays = (pind3d_cmp, oind3d_cmp, obj3d_cmp)
# vals = (pind′, oind, o)
#
# println("# of objects = $(length(ovec))")
# # @time assign_val_shape!((arrays..., @view(param3d_cmp[:,:,:,nw,nw])), (vals..., param[nw,nw]), shape, τlcmp)
# # @code_warntype assign_val_shape!((arrays..., param3d_cmp), (vals, param), shape, τlcmp)

# assign_val_shape!(pind3d_cmp, pind′, shape, τlcmp)
# assign_val_shape!(oind3d_cmp, oind, shape, τlcmp)
# assign_val_shape!(obj3d_cmp, o, shape, τlcmp)
# if nw == 4
#     assign_val_shape!(param3d_cmp, param, shape, τlcmp)
# else  # nw = 1, 2, 3
#     assign_val_shape!(@view(param3d_cmp[:,:,:,nw,nw]), param[nw,nw], shape, τlcmp)
# end


# @code_warntype MaxwellFDM.assign_param_cmp!(gt, nw, param3d_cmp, obj3d_cmp, pind3d_cmp, oind3d_cmp, ovec, τlcmp)
# @code_warntype assign_param!(param3d, obj3d, pind3d, oind3d, ovec, g3.ghosted.τl, g3.isbloch)

boundft = SVector(EE,EE,EE)
@time assign_param!((ε3d,μ3d), (εobj3d,μobj3d), (εind3d,μind3d), (εoind3d,μoind3d), boundft, ovec, g3.ghosted.τl, g3.isbloch)
ε3d_assigned = copy(ε3d)

# gt_cmp′ = alter.(gt_cmp)
# lcmp = MaxwellFDM.t_ind(g3.l, gt_cmp)
# σcmp = MaxwellFDM.t_ind(g3.σ, gt_cmp)
# lcmp′ = MaxwellFDM.t_ind(g3.ghosted.l, gt_cmp′)
# ∆τcmp′ = MaxwellFDM.t_ind(g3.ghosted.∆τ, gt_cmp′)
#
# param3d_gt = param3d[ngt]
# obj3d_cmp′ = obj3d[ngt][nw]
# pind3d_cmp′ = pind3d[ngt][nw]
# oind3d_cmp′ = oind3d[ngt][nw]

# @code_warntype MaxwellFDM.smooth_param_cmp!(gt, nw, param3d_gt, obj3d_cmp′, pind3d_cmp′, oind3d_cmp′, lcmp, lcmp′, σcmp, ∆τcmp′)

ft = EE
@time smooth_param!(ε3d, εobj3d, εind3d, εoind3d, ft, boundft, g3.l, g3.ghosted.l, g3.σ, g3.ghosted.∆τ)
ε3d_smoothed = copy(ε3d)
@save "benchmark/smoothing_result.jld2" ε3d_assigned ε3d_smoothed

# # Construct arguments and call assign_param!.
# kd = KDTree(ovec)
# param3d = create_default_param3d(g3.N)
# @time assign_param!(param3d, kd, g3.l, g3.lg, g3.N, g3.L, g3.isbloch)
# @btime assign_param!(param3d, kd, g3.l, g3.lg, g3.N, g3.L, g3.isbloch)
# @profile assign_param!(param3d, kd, g3.l, g3.lg, g3.N, g3.L, g3.isbloch)
# @code_warntype assign_param!(param3d, kd, g3.l, g3.lg, g3.N, g3.L, g3.isbloch)
