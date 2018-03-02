using MaxwellFDM, StaticArrays
using BenchmarkTools

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

L₀ = 1e-9
unit = PhysUnit(L₀)
Npml = ([10,10,10], [10,10,10])
ebc = [BLOCH, BLOCH, BLOCH]

g3 = Grid(unit, lprim, Npml, ebc)

# Create materials.
εvac = 1.0
vac = EncodedMaterial(PRIM, Material("Vacuum", ε=εvac))

εdiel = 12.0
diel = EncodedMaterial(PRIM, Material("Dielectric", ε=εdiel))

# Create objects.
dom_vac = Object(Box(g3.bounds), vac)
obj_diel = Object(Sphere([0,0,0], 75), diel)
# obj_diel = Object(Box([0,0,0], [75,75,75]), diel)

# Add objects.
ovec = Object3[]
paramset = (SMat3Complex[], SMat3Complex[])
add!(ovec, paramset, dom_vac, obj_diel)
# add!(ovec, paramset, dom_vac)

param3d = create_param3d(g3.N)
obj3d = create_n3d(Object3, g3.N)
pind3d = create_n3d(ParamInd, g3.N)
oind3d = create_n3d(ObjInd, g3.N)

τl = g3.ghosted.τl

gt = PRIM
nw = 4
ngt = Int(gt)
ngt′ = alter(ngt)

M = length.(τl[nPR])  # N+1
ind = (map((m,e)->circshift(1:m, e==BLOCH), M, g3.ebc.data),
       map((m,e)->circshift(1:m, -(e==BLOCH)), M, g3.ebc.data))

gt_cmp = SVector(gt, gt, gt)
gt_cmp = map((k,g)->(k==nw ? alter(g) : g), nXYZ, gt_cmp)  # no change if nw = 0

# Choose the circularly shifted indices to use.
ind_cmp = t_ind(ind, gt_cmp)

# Prepare the circularly shifted locations of the field components.
τlcmp = view.(t_ind(τl,gt_cmp), ind_cmp)  # Tuple3{Vector{Float}}: locations of Fw = Uw or Vw

# Prepare the circularly shifted viewes of various arrays to match the sorted
# τl.  Even though all arrays are for same locations, param3d_cmp contains gt
# material, whereas obj3d_cmp, pind3d_cmp, oind3d_cmp contain alter(gt)
# material, so use ngt′ instead of ngt for them.
param3d_cmp = view(param3d[ngt], ind_cmp..., 1:3, 1:3)
obj3d_cmp = view(obj3d[ngt′][nw], ind_cmp...)
pind3d_cmp = view(pind3d[ngt′][nw], ind_cmp...)
oind3d_cmp = view(oind3d[ngt′][nw], ind_cmp...)

# @code_warntype MaxwellFDM.assign_param_obj!(gt, nw, param3d_cmp, obj3d_cmp, pind3d_cmp, oind3d_cmp, ovec[1], τlcmp)
# @code_warntype MaxwellFDM.assign_param_cmp!(gt, nw, param3d_cmp, obj3d_cmp, pind3d_cmp, oind3d_cmp, ovec, τlcmp)
# @code_warntype assign_param!(param3d, obj3d, pind3d, oind3d, ovec, g3.ghosted.τl, g3.ebc)

@time assign_param!(param3d, obj3d, pind3d, oind3d, ovec, g3.ghosted.τl, g3.ebc)

gt_cmp′ = alter.(gt_cmp)
lcmp = t_ind(g3.l, gt_cmp)
σcmp = t_ind(g3.σ, gt_cmp)
lcmp′ = t_ind(g3.ghosted.l, gt_cmp′)
∆τcmp′ = t_ind(g3.ghosted.∆τ, gt_cmp′)

param3d_gt = param3d[ngt]
obj3d_cmp′ = obj3d[ngt][nw]
pind3d_cmp′ = pind3d[ngt][nw]
oind3d_cmp′ = oind3d[ngt][nw]

# @code_warntype MaxwellFDM.smooth_param_cmp!(gt, nw, param3d_gt, obj3d_cmp′, pind3d_cmp′, oind3d_cmp′, lcmp, lcmp′, σcmp, ∆τcmp′)

@time smooth_param!(param3d, obj3d, pind3d, oind3d, g3.l, g3.ghosted.l, g3.σ, g3.ghosted.∆τ)


# # Construct arguments and call assign_param!.
# kd = KDTree(ovec)
# param3d = create_default_param3d(g3.N)
# @time assign_param!(param3d, kd, g3.l, g3.lg, g3.N, g3.L, g3.ebc)
# @btime assign_param!(param3d, kd, g3.l, g3.lg, g3.N, g3.L, g3.ebc)
# @profile assign_param!(param3d, kd, g3.l, g3.lg, g3.N, g3.L, g3.ebc)
# @code_warntype assign_param!(param3d, kd, g3.l, g3.lg, g3.N, g3.L, g3.ebc)
