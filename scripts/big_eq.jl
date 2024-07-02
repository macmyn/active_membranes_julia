# open("shapef.m","r") do f
#   global a = read(f,String)
#   end
  
file_path = "shapef.m"
a = []

open(file_path) do file
    for line in eachline(file)
        push!(a, line)
    end
end


# Change 2., 3., 4. to 2.0., 3.0, 4.0 (except for the variable Ga2)
replace_2_regex = r"(?<!Ga)2\.\*"  # 2.* --> 2.0.* except when it's part of Ga2.*
a = replace.(a, replace_2_regex => "2.0*")
a = replace.(a, r"2\.\+" => "2.0.+")
a = replace.(a, r"2\.\/" => "2.0./")

a = replace.(a, r"3\.\*" => "3.0*")
a = replace.(a, r"4\.\*" => "4.0*")

# Replace y(N,:) with y[N,:]
a = replace.(a, r"y\((\d+),:\)" => s"y[\1,:]")

# Replace dy(N,:) with dy[N,:]
a = replace.(a, r"dy\((\d+),:\)" => s"dy[\1,:]")

a = replace.(a, r"\.\.\.\s*\n\s*" => "\n")
a = replace.(a, r"\.\.\." => "")

a = replace.(a, r"%" => "#")

a = replace.(a, "~=" => "!=")

# Replace spline function with Julia equivalent
old_var_pattern = r"spline\(([^,]+),([^,]+),([^)]+)\)"
matches = match(old_var_pattern, matlab_line)
uc_var, r0c_var, u_var = matches.captures
new_var_pattern = "itp = interpolate(($uc_var,), $r0c_var, BSpline(Cubic(Line(OnGrid()))))\n" * "r0 = itp.($u_var)  # Interpolate at the points in $u_var"


# Replace MATLAB-style spline function call with Julia-style code
a = replace.(a, old_var_pattern => new_var_pattern)


# Example usage
matlab_line = "r0 = spline(uc,r0c,u);"



a[1] = "function shapef(dy, y, u, t)"

push!(a, "\nend")



# open("big_eq.txt", "w") do file
#   write(file, a);
# end


# println(a)

open("big_eq2.jl", "w") do file
  for line in a
    println(file, line)
  end
end