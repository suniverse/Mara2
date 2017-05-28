local mara = require 'mara'

local fieldMagnitude = 1.0

local function background(x, y, z)
	return {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, fieldMagnitude}
end

local function abc_driving_floor_ceil(x, y, z)
	local sgnz = z > 0 and 1 or -1
	local vx = 0.2 * math.sin(4 * math.pi * y) * sgnz
	local vy = 0.2 * math.cos(4 * math.pi * x) * sgnz
	return {vx, vy}
end

local function abc_driving_floor_only(x, y, z)
	if (z < 0) then
		local vx = 0.2 * math.sin(4 * math.pi * y)
		local vy = 0.2 * math.cos(4 * math.pi * x)
		return {vx, vy}
	else
		return {0, 0}
	end
end

local function single_vortex_pair(x, y, z)
	local sgnz = z > 0 and 1 or -1
	local R = (x^2 + y^2)^0.5
	local k = 5
	local v = 0.2 / k
	local vf = v * 3 * k * (k * R)^2 * math.exp (-(k * R)^3) -- v-phi
	local vx = vf * (-y / R) * sgnz
	local vy = vf * ( x / R) * sgnz
	return {vx, vy}
end

local setup = {
	run_name = 'BoundaryDrivenMHD',
	output_directory = 'data/bdrv-vortex-pair-164z',
	final_time = 64.,
	checkpoint_interval = 0.0,
	vtk_output_interval = 0.1,
	cfl_parameter = 0.3,
	initial_data = background,
	grid_geometry = 'cartesian',
	resolution = {16, 16, 64},
	domain_lower = {-0.5, -0.5, -1.0},
	domain_upper = { 0.5,  0.5,  1.0},
	conservation_law = {'newtonian_mhd'},
	riemann_solver = 'hlle',
	flux_scheme = {'method_of_lines_plm', plm_theta=1.5},
	runge_kutta_order = 2,
	boundary_condition = 'driven_mhd',
	boundary_velocity_function = single_vortex_pair,
	disable_ct = false,
	pressure_floor = 1e-2
}

mara.run(setup)