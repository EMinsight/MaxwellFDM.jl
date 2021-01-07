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

M = (length(xprim)-1, length(yprim)-1)
println("M = $M")

lprim = (xprim, yprim)

isbloch = [true, true]

g2 = Grid(lprim, isbloch)

# Create materials.
εvac = 1.0
vac = Material{2,1}("Vacuum", ε=εvac)  # Ke = 2, Km = 1

εdiel = 12.0
diel = Material{2,1}("Dielectric", ε=εdiel)  # Ke = 2, Km = 1

# Create objects.
dom_vac = Object(Box(g2.bounds), vac)
obj_diel = Object(Sphere([0,0], 75), diel)
# obj_diel = Object(Box([0,0,0], [75,75,75]), diel)
obj_xn_diel = Object(Sphere([-150,0], 75), diel)
obj_xp_diel = Object(Sphere([150,0], 75), diel)
obj_yn_diel = Object(Sphere([0,-150], 75), diel)
obj_yp_diel = Object(Sphere([0,150], 75), diel)

# Add objects.
oind2shp = Shape2[]
oind2εind = ParamInd[]
oind2μind = ParamInd[]
εind2ε = SSComplex2[]
μind2μ = SSComplex1[]

add!(oind2shp, (oind2εind,oind2μind), (εind2ε,μind2μ), dom_vac, obj_diel, obj_xn_diel, obj_xp_diel, obj_yn_diel, obj_yp_diel)


N = g2.N
ε2d = create_param_array(N, ncmp=2)
εxx_oind2d = create_oind_array(N)
εyy_oind2d = create_oind_array(N)
εoo_oind2d = create_oind_array(N)

μ2d = create_param_array(N, ncmp=1)
μzz_oind2d = create_oind_array(N)

boundft = SVector(EE,EE)
@time begin
    # Set the diagonal entries of ε2d.
    assign_param!(ε2d, (εyy_oind2d,εxx_oind2d), ft2gt.(EE,boundft), oind2shp, oind2εind, εind2ε, g2.ghosted.τl, g2.isbloch)

    # Set the off-diagonal entries of ε2d.
    assign_param!(ε2d, tuple(μzz_oind2d), ft2gt.(EE,boundft), oind2shp, oind2εind, εind2ε, g2.ghosted.τl, g2.isbloch)

    # Set μ2d (2D array of scalars).
    assign_param!(μ2d, tuple(εoo_oind2d), ft2gt.(HH,boundft), oind2shp, oind2μind, μind2μ, g2.ghosted.τl, g2.isbloch)
end

@time begin
    # Smooth the diagonal entries of ε2d.
    smooth_param!(ε2d, (εxx_oind2d,εyy_oind2d), oind2shp, oind2εind, εind2ε, ft2gt.(EE,boundft), g2.l, g2.ghosted.l, g2.σ, g2.ghosted.∆τ)

    # Smooth the off-diagonal entries of ε2d.
    smooth_param!(ε2d, tuple(εoo_oind2d), oind2shp, oind2εind, εind2ε, ft2gt.(EE,boundft), g2.l, g2.ghosted.l, g2.σ, g2.ghosted.∆τ)
end


# # Construct arguments and call assign_param!.
# kd = KDTree(oind2obj)
# param2d = create_default_param2d(g2.N)
# @time assign_param!(param2d, kd, g2.l, g2.lg, g2.N, g2.L, g2.isbloch)
# @btime assign_param!(param2d, kd, g2.l, g2.lg, g2.N, g2.L, g2.isbloch)
# @profile assign_param!(param2d, kd, g2.l, g2.lg, g2.N, g2.L, g2.isbloch)
# @code_warntype assign_param!(param2d, kd, g2.l, g2.lg, g2.N, g2.L, g2.isbloch)
