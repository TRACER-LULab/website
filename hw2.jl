### A Pluto.jl notebook ###
# v0.17.7

using Markdown
using InteractiveUtils

# ╔═╡ 1ee74d1b-f96a-41e7-9c0b-d552cfc3f2c6
begin
	using Plots
	using LinearSolve
	using Symbolics
	using Tullio
	using Latexify
	using Unitful
	using PlutoUI
	using LinearAlgebra
	using LabelledArrays
	using ModelingToolkit
	using NonlinearSolve
	@variables x y z F_x F_y F_z
end;

# ╔═╡ d3a088a2-6495-4750-85f8-a01f9092bf6e
md"""
# ENGR 727 
## Homework 2
## Carson Farmer
## February 1, 2022
"""

# ╔═╡ d018bb86-7e27-11ec-3db2-cd579d41ad35
md"""
# Imports
"""

# ╔═╡ 79a774e2-14b4-4079-a673-404d7402a512
md"""
# Problem 1 - Problem 1.14 in Textbook
"""

# ╔═╡ a7f73c0c-8af7-4679-a42a-49652a554ff1
begin
	τ_1 = [-2x^2+3y^2-5z z+4x*y-7 -3x+y+1
		 	z+4x*y-7 -2y^2 0
		 	-3x+y+1 0 3x+y+3z-5
			]
end

# ╔═╡ bef16de0-e39b-425c-b0e4-37434730b48e
begin
	cart = [x, y, z]
	cart′ = Differential.(cart)
	F = [F_x, F_y, F_z]
	@tullio τ[i] := cart′[j](τ_1[i, j])
	eqs_1 = F .~ -τ
	latexify(eqs_1)
end

# ╔═╡ 65d33acc-6b38-4f33-8da9-085817f34bea
md"""
!!! note "Answer"
	$F_x = 0$
	$F_x = 0$
	$F_z = 0$
"""

# ╔═╡ 4cafc5a2-672f-47c3-b11a-445ec7a757e6
md"""
# Problem 2 - Problem 1.65 in Textbook
"""

# ╔═╡ 5393419e-026d-425a-8f93-5ea593202120
md"""
```math
\tau_{ij} = 
\begin{bmatrix}
	60 & 40 & -40\\
	40 &  0 & -20\\
	-40&-20 & 20
\end{bmatrix}
```
"""

# ╔═╡ 3fd8a785-96a3-451f-b383-6c87555eb1dc
begin
	@variables σ[1:3,1:3] σ_x σ_y σ_z τ_xy τ_xz τ_yz
	sigma = [σ_x τ_xy τ_xz; τ_xy σ_y τ_yz; τ_xz τ_yz σ_z]
	σ_2 = [60 40 -40;40 0 -20; -40 -20 20]
end;

# ╔═╡ 7de1b874-df71-43c9-849d-e421b3fec528
md"""
```math
R = 
\begin{bmatrix}
	\cos{\theta} & \cos{\frac{\pi}{2}-\theta} & 0\\
	\cos{\frac{3\pi}{2}-\theta} & \cos{\theta} & 0\\
	0 & 0 & 1
\end{bmatrix}
=
\begin{bmatrix}
	0.866 & 0.5 & 0.5\\
	-0.5  & 0.866 & 0\\
	0 & 0 & 1
\end{bmatrix}
```
"""

# ╔═╡ 80a4c591-97d5-4f0e-918a-c9520e9d0f3b
begin
	θ = 30*π/180
	r = round.([cos(θ) cos(π/2-θ) 0; cos(-(π/2+θ)) cos(θ) 0; 0 0 1], sigdigits=10)
end;

# ╔═╡ 7718d209-08bd-455b-a2df-b4deea805b61
md"""
```math
\tau_{ij}^\prime = R \sigma R^T = 
\begin{bmatrix}
	79.64 & -5.981 & -44.64 \\
	-5.981 & -19.64 & 2.679 \\
	-44.64 & 2.679 & 20.0
\end{bmatrix}
```
```math

{\tau_{ij}^\prime}^2 = R \sigma^2 R^T = 
\begin{bmatrix}
	8371.28 & -478.461 & -4464.1 \\
	-478.461 & 428.719 & 267.949 \\
	-4464.1 & 267.949 & 2400.0\\
\end{bmatrix}
```
"""

# ╔═╡ 86464ad1-ef40-4566-b1cb-5c4d7834107f
sol_2 = r*σ_2*r';

# ╔═╡ f4c74e07-a30b-4c3e-9f6f-30aa05e1e69e
md"""
$I_1 = tr(\tau_{ij}) = tr(\tau_{ij}^\prime) = 60+0+20 = 79.64-19.64-20 = 80$
$I_2 = \frac{1}{2}(tr(\tau_{ij})^2-tr(\tau_{ij}^2)) = \frac{1}{2}(tr(\tau_{ij}^\prime)^2-tr({\tau_{ij}^\prime}^2)) = 0.5(80^2-(6800+2000+2400)) = 0.5(80^2-(8371.28+428.719+2400)) = -2400$
$I_3 = |\tau_{ij}| = |\tau_{ij}^\prime| = (-40\cdot(40\cdot20+20\cdot40)+20\cdot(60\cdot-20)+40\cdot40) = 79.64\cdot(-19.64\cdot20-2.679^2)+5.981(-5.981\cdot2.679+2.679\cdot44.64)-44.64\cdot(-5.981\cdot2.679-44.64\cdot19.64) = 8000$
"""

# ╔═╡ 6eed76e7-24bf-4ca9-bf17-a83b9bc755de
begin
	I₁ = tr(sol_2)
	I₁ = tr(σ_2)
	I₂ = 1/2*(tr(sol_2)^2-tr(sol_2^2))
	# I₂ = 1/2*(tr(σ_2)^2-tr(σ_2^2))
	I₃ = det(σ_2)
	# I₃ = det(sol_2)
end;

# ╔═╡ 130ecbe0-1cfe-4547-84f5-094bf8044126
md"""
!!! note "Answer"
	$τ_{ij} = \begin{bmatrix}
		79.64 & -5.981 & -44.64\\
		-5.984& -19.64 & 2.679 \\
		-44.64& 2.679  & 20.0 
		\end{bmatrix}\text{MPa}$
	$I_1 = 80\text{MPa}$
	$I_2 = -2400\text{MPa}^2$
	$I_3 = 8000\text{MPa}^3$
"""

# ╔═╡ dbce3d98-51e3-456c-a9b2-1457ba62c700
md"""
# Problem 3 - Problem 1.69 in Textbook
"""

# ╔═╡ 1e3cd5d8-3075-4746-8637-c912d616ece2
md"""
```math
\tau_{ij} = 
\begin{bmatrix}
	12 & 6 & 9\\
	6  & 10& 3\\
	9  & 3 &14
\end{bmatrix}
```
"""

# ╔═╡ 93f83718-2ba3-4988-b5e1-4ef0b6d92039
md"""
$\sigma_{1,2,3} = |\tau_{ij} - \sigma_p I|$
"""

# ╔═╡ b97ee3f5-a198-4ba2-a087-70ef315db537
begin
	@variables s σ_p σ_123
	σ_3 = [12 6 9;6 10 3; 9 3 14]
	latexify(σ_123 ~ PolyForm(expand(det(σ_3-σ_p*I(3)))))
end

# ╔═╡ db3b3cbb-39a6-42e0-921a-bfa53e805eea
md"""
Solving for the roots of the equation yields: 

$\sigma_1 = 24.75\text{MPa}$

$\sigma_2 = 8.478\text{MPa}$

$\sigma_3 = 2.773\text{MPa}$
"""

# ╔═╡ 54520742-17a2-4a66-b7f5-feda3852745e
md"""
To find the direction cosines, the following equations are solved for each of the principle stresses: 

$(\sigma_x-\sigma_p)l+\tau_{xy}m+\tau_{xz}n = (12-\sigma_p)l+6m+9n = 0$

$\tau_{xy}l + (\sigma_y-\sigma_p)m + \tau_{yz}n = 6l + (10-\sigma_p)m + 3n = 0$

$\tau_{xz}l+ \tau_{yz}m + (\sigma_z-\sigma_p)n = 9l+ 3m + (14-\sigma_p)n = 0$

$l^2+m^2+n^2 = 0$
"""

# ╔═╡ b76a4896-f87e-4176-a28f-d9fbda2b7f19
begin
	@variables l m n
	direction_eqns = [
		0 ~ (12-σ_p)*l + 6*m + 9*n,
		0 ~ 6*l + (10-σ_p)*m + 3*n,
		0 ~ 9*l + 3*m + (14-σ_p)*n,
	]
end;

# ╔═╡ f88ce048-e1c1-47c0-900d-bd81739dfd67
md"""
For σₚ = σ₁ = 24.75MP:
"""

# ╔═╡ 3ccbb53c-80ab-45b2-8244-e96c225d4494
md"""
$-12.75l + 6m + 9n = 0$
$6l - 14.75m + 3n = 0$
$9l + 3m - 10.75n = 0$
"""

# ╔═╡ fc9dfb09-fb8e-4683-b0dc-4c1375744126
md"""
Solving for $l$, $m$, and $n$, the direction vector is $\langle -0.6466, -0.3957, -0.6522 \rangle$
"""

# ╔═╡ 6e09578e-e23f-4cfb-a813-10c5aa8be359
md"""
For σₚ = σ₂ = 8.478 MPa:
"""

# ╔═╡ bc07861d-b324-4a1c-bf77-c1290508baf4
md"""
$3.522l + 6m + 9n = 0$
$6l + 1.522m + 3n = 0$
$9l + 3m + 5.22n = 0$
"""

# ╔═╡ 8de4abab-f5e1-432e-902f-162169468d73
md"""
Solving for $l$, $m$, and $n$, the direction vector is $\langle -0.0809, -0.8144, 0.5746 \rangle$
"""

# ╔═╡ c2fa4f8f-3689-4a5d-87ee-d182ef250f4a
md"""
For σₚ = σ₃ = 2.773MPa:
"""

# ╔═╡ 0651538b-cf6e-4f92-87b2-e4ec940e1aa3
md"""
$9.227l + 6m + 9n = 0$
$6l + 7.227m + 3n = 0$
$9l + 3m + 11.227n = 0$
"""

# ╔═╡ 583a0ab9-db59-43be-8583-d97d39a8c486
md"""
Solving for $l$, $m$, and $n$, the direction vector is $\langle -0.7585, 0.4243, 0.4946 \rangle$
"""

# ╔═╡ 51ab8601-c9c5-44b9-b041-0055731e797a
md"""
!!! note "Answer" 
	$\sigma_1 = 24.75\text{MPa},\ \ \ \vec{n} = \langle -0.6466, -0.3957, -0.6522 \rangle$
	$\sigma_1 = 8.478\text{MPa},\ \ \ \vec{n} = \langle -0.0809, -0.8114,  0.5746 \rangle$
	$\sigma_1 = 2.773\text{MPa},\ \ \ \vec{n} = \langle -0.7585,  0.4243,  0.4946 \rangle$
"""

# ╔═╡ 617ac012-e2d0-4708-af65-8cd21c9d45ee
md"""
# Problem 4 - Problem 1.71 in Textbook
The given stresses, $\tau_{ij}$, are:
"""

# ╔═╡ 3ca25a1c-6753-4015-976f-77ccc806eba8
τ_4_sym = let
	[x^2+y 0 0;
	 0 y^2-5 0;
	 0 0 -x+6y+z]
end

# ╔═╡ 08674b25-909b-42ac-a487-0f5c56e1616b
md"""
in (x, y, z) coordinates. At point (3, 1, 5) the stress is:
```math
\tau_{ij} = \begin{bmatrix}
	10 & 0 & 0\\
	0 & -4 & 0\\
	0 & 0 & 8 
\end{bmatrix}\text{MPa}
```
"""

# ╔═╡ bff27eec-f0eb-4cce-b4d6-70fd7b4ec549
τ_4 = let
	substitute.(τ_4_sym, (Dict(x=>3, y=>1, z=>5), )).|>Symbolics.value
end;

# ╔═╡ 62ba9781-a790-4149-a399-0e06ad9e52d5
md"""
## Subproblem a
The rotation matrix is: 
```math
R = \begin{bmatrix}
	1 & 0 & 0 \\
	0 & \frac{1}{2} & \frac{\sqrt{3}}{2}\\
	0 & -\frac{\sqrt{3}}{2} & \frac{1}{2}\\
\end{bmatrix}
```
"""

# ╔═╡ 493608ce-70ee-468f-a2d4-f8f443adf7df
R_4a = [1 0 0 ; 0 1/2 sqrt(3)/2; 0 -sqrt(3)/2 1/2];

# ╔═╡ decad7b8-75d7-4a01-884b-bfafab73cf32
round.(R_4a*τ_4*R_4a', sigdigits = 3);

# ╔═╡ 3adc9eab-9642-494a-9ea2-0c0edd79935e
md"""
!!! note "Answer"
	```math
	\tau_{ij}^\prime = R^T\ \tau_{ij} R =
	\begin{bmatrix}
		10 & 0 & 0 \\
		0 & 5 & 5.2 \\
		0 & 5.2 & -1
	\end{bmatrix}\text{MPa}
	```
"""

# ╔═╡ 66d9df69-7996-4493-a516-317ac6f613d9
md"""
## Subproblem b
The given rotation matrix is:
```math
R = \begin{bmatrix}
	\frac{2}{\sqrt{5}} & -\frac{1}{\sqrt{5}} & 0\\
	\frac{1}{\sqrt{5}} & \frac{2}{\sqrt{5}} & 0\\
	0 & 0 & 1
\end{bmatrix}
```
"""

# ╔═╡ be145b51-6cb6-4628-8447-b8d17a703d96
R_4b = let
	[2/sqrt(5) -1/sqrt(5) 0; 1/sqrt(5) 2/sqrt(5) 0; 0 0 1]
end;

# ╔═╡ b0c491d5-30aa-4c85-ab41-9016e8cc64de
τ_4′′ = round.(Symbolics.value.(R_4b*τ_4*R_4b'), sigdigits = 2);

# ╔═╡ a59e79c6-a86f-4d60-b38a-090e99aae5e6
md"""
By applying the same matrix rotation operation as in Subproblem a, the resulting stress state is:
"""

# ╔═╡ 211ed5cf-6338-4b6c-ae9e-ffb3fd1d3a50
md"""
!!! note "Answer"
	```math
	\tau_{ij}^{\prime\prime} = 
	\begin{bmatrix}
		7.2 & 5.6 & 0\\
		5.6 & -1.2 & 0\\
		0 & 0 & 8
	\end{bmatrix}\text{MPa}
	```
"""

# ╔═╡ 230b2b41-802a-46d2-9d92-62ac9922c43d
md"""
## Subproblem c
"""

# ╔═╡ ee92bd0b-29cc-4fe7-a132-78568743fcb5
md"""
!!! note "Answer"
	For I₁:

	$I_1(\sigma_4) = 10-4+8 = 14\text{MPa}$
	$I_1(\sigma_4^\prime) = 10+5-1.0 = 14\text{MPa}$
	$I_1(\sigma_4^{\prime\prime}) = 7.2-1.2+8.0 = 1\text{MPa}$

	For I₂: 

	$I_2(\sigma_4) = 1/2(14^2-(100+16+64)) = 0.03125\text{MPa}^2$
	$I_2(\sigma_4^\prime) = 1/2(14^2-(100+52.04+28.04)) = 0.03125\text{MPa}^2$
	$I_2(\sigma_4^{\prime\prime} = 1/2(14^2-(83.2+32.8+64.0)) = 0.03125\text{MPa}^2$

	For I₃: 

	$I_3(\sigma_4) = \det(\sigma_4) = 10\cdot -4 \cdot 8 =  -320\text{MPa}^3$
	$I_3(\sigma_4^\prime) = \det(\sigma_4^\prime) = 10\cdot(-5-5.2^2) = -320\text{MPa}^3$
	$I_3(\sigma_4^{\prime\prime}) = \det(\sigma_4^{\prime\prime}) = 8\cdot(7.2\cdot-1.2-5.6^2) = -320\text{MPa}^3$
"""

# ╔═╡ fe20af7d-05a4-450d-b7ab-7003055586df
md"""
# Problem 5 - Problem 1.75 in Textbok
"""

# ╔═╡ 9729b10d-6c6d-465e-aff1-38eaf7fca34c
md"""
The general equation of a plane through 3 points (r₁, r₂, r₃) is: 

$f(x, y, z) = (\langle x, y, z \rangle - \vec{r}_1)\cdot(\vec{r}_2-\vec{r}_1)\times(\vec{r}_3-\vec{r}_1) = 0$

From the equation of a plane, $ax+by+cz+d = 0$, the cosines vector, \vec{n} is:

$\vec{n} = \langle \frac{a}{\sqrt{a^2+b^2+c^2}}, \frac{b}{\sqrt{a^2+b^2+c^2}}, \frac{c}{\sqrt{a^2+b^2+c^2}}\rangle$
The normal and shear stress on the plane, $\sigma^\prime$ and $\tau^\prime$ respectively, are defined as:

$\sigma^\prime = \vec{n}^T\sigma\vec{n}$
$\tau^\prime = (\vec{n}^T\sigma^T\sigma\vec{n}-{\sigma^\prime}^2)^{\frac{1}{2}}$

where σ is the stress tensor.

For this problem: 
```math
\sigma = \begin{bmatrix}
	100 & 40 &  0\\
	 40 & 60 & 80\\
	  0 & 80 & 20 
\end{bmatrix}
```
"""

# ╔═╡ cfb974e3-e281-478d-907b-15ec708af696
σ_5 = [100 40 0;40 60 80;0 80 20]u"MPa";

# ╔═╡ 23ecb4ec-92a5-465c-b425-058753c9b6f2
md"""
The following function is used to solve for the normal and shear stress on the oblique plane:
"""

# ╔═╡ 6c99da8e-c83a-4552-8e3e-1710e4932bba
function oblique_plane(r1, r2, r3, σ)
	@variables x y z 
	xs = [x, y, z]
	plane = 0~(xs-r1)⋅((r2-r1)×(r3-r1))
	coeffs, constant = Symbolics.polynomial_coeffs(Symbolics.PolyForm(plane.rhs), xs)
	for c in xs
		if !haskey(coeffs, c)
			coeffs[c] = 0.0
		end
	end
	dir = SLVector((x = coeffs[x], y = coeffs[y], z = coeffs[z]))
	cosines = [dir.x/norm(dir), dir.y/norm(dir), dir.z/norm(dir)]
	σ′ = cosines'*σ*cosines
	τ′ = sqrt(cosines'*σ'*σ*cosines-σ′^2)
	(σ = σ′, τ = τ′)
end

# ╔═╡ 224242da-8ff5-4485-a8aa-2b15820f8841
md"""
## Subproblem a
$\vec{r}_1 = \langle 0, 4, 2 \rangle$
$\vec{r}_2 = \langle 0, 4, 0 \rangle$
$\vec{r}_3 = \langle 3, 0, 2 \rangle$
"""

# ╔═╡ 93417d49-684d-4bb5-b220-c95dec571035
res_5a = oblique_plane([0, 4, 2], [0, 4, 0], [3, 0, 2], σ_5)

# ╔═╡ 6533e757-73f7-4f3c-a20e-a9e90b9943b7
md"""
!!! note "Answer"
	$\sigma^\prime = 124\text{MPa}$
	$\tau^\prime = 48.6621\text{MPa}$
"""

# ╔═╡ 3375015f-a720-4eb2-9c67-f133da64da54
md"""
## Subproblem b
$\vec{r}_1 = \langle 0, 0, 2 \rangle$
$\vec{r}_2 = \langle 3, 0, 2 \rangle$
$\vec{r}_3 = \langle 0, 4, 0 \rangle$
"""

# ╔═╡ 34c79484-8cf6-48ef-a08e-c0a0e54f32e9
res_5b = oblique_plane([0, 0, 2], [3, 0, 2], [0, 4, 0], σ_5)

# ╔═╡ fee104ad-6f7d-40d7-81e9-dd665e6532b5
md"""
!!! note "Answer"
	$\sigma^\prime = 92.0\text{MPa}$
	$\tau^\prime = 66.453\text{MPa}$
"""

# ╔═╡ 9df061be-bda1-45e3-a8c0-062b7c0b7e4b
md"""
## Subproblem c
$\vec{r}_1 = \langle 0, 0, 2 \rangle$
$\vec{r}_2 = \langle 0, 4, 0 \rangle$
$\vec{r}_3 = \langle 3, 0, 0 \rangle$
"""

# ╔═╡ 913f2bf9-0d2a-49ec-b6c9-3b87684f595c
res_5c = oblique_plane([0, 0, 2], [0, 4, 0], [3, 0, 0], σ_5)

# ╔═╡ 93a8de29-e4c5-4f2b-acb8-8b66712a45d0
md"""
!!! note "Answer"
	$\sigma^\prime = 109.836\text{MPa}$
	$\tau^\prime = 74.272\text{MPa}$
"""

# ╔═╡ 22592932-d897-4ad4-8549-a322ea2eb1c8
md"""
# Problem 6
The stress state is given by:

```math
[\tau_{ij}] = \begin{bmatrix}
4 & 1 & 2 \\
1 & 6 & 0 \\
2 & 0 & 8
\end{bmatrix}
```
and a reference system of: 

|  | x | y | z |
| :----: | :----: | :----: | :----: |
|x′ |0 | $\frac{\sqrt{2}}{2}$| $-\frac{\sqrt{2}}{2}$ |
|y′ | $-\frac{\sqrt{3}}{3}$ | $\frac{\sqrt{3}}{3}$| $\frac{\sqrt{3}}{3}$ |
|z′ | $\frac{\sqrt{6}}{3}$ | $\frac{\sqrt{6}}{6}$| $\frac{\sqrt{6}}{6}$ |
"""

# ╔═╡ 4f3391ca-1e2d-4f78-85b2-45d07ab6db6d
md"""
## Subproblem A - Do the stress transformation with the equations provided
"""

# ╔═╡ 7b8d5f1f-5b28-40da-9e30-0e060e1a2d8b
let
	R = [0 sqrt(2)/2 -sqrt(2)/2;-sqrt(3)/3 sqrt(3)/3 sqrt(3)/3; sqrt(6)/3 sqrt(6)/6 sqrt(6)/6]
	τij = [4 1 2; 1 6 0; 2 0 8]
	
	l = R[:, 1] # [l₁, l₂, l₃]
	m = R[:, 2] # [m₁, m₂, m₃]
	n = R[:, 3] # [n₁, n₂, n₃]
	
	σx′ = τij[1, 1]*l[1]^2+τij[2, 2]*m[1]^2+τij[3,3]*n[1]^2+2*(τij[1, 2]*l[1]*m[1]+τij[2,3]*m[1]*n[1]+τij[1,3]*l[1]*n[1])
	σy′ = τij[1, 1]*l[2]^2+τij[2, 2]*m[2]^2+τij[3,3]*n[2]^2+2*(τij[1, 2]*l[2]*m[2]+τij[2,3]*m[2]*n[2]+τij[1,3]*l[2]*n[2])
	σz′ = τij[1, 1]*l[3]^2+τij[2, 2]*m[3]^2+τij[3,3]*n[3]^2+2*(τij[1, 2]*l[3]*m[3]+τij[2,3]*m[3]*n[3]+τij[1,3]*l[3]*n[3])

	τx′y′ = τij[1,1]*l[1]*l[2]+τij[2,2]*m[1]*m[2]+τij[3,3]*n[1]*n[2]+τij[1,2]*(l[1]*m[2]+m[1]*l[2])+τij[2,3]*(m[1]*n[2]+m[2]*n[1])+τij[1,3]*(n[1]*l[2]+n[2]*l[1])
	τx′z′ = τij[1,1]*l[1]*l[3]+τij[2,2]*m[1]*m[3]+τij[3,3]*n[1]*n[3]+τij[1,2]*(l[1]*m[3]+m[1]*l[3])+τij[2,3]*(m[1]*n[3]+m[3]*n[1])+τij[1,3]*(n[1]*l[3]+n[3]*l[1])
	τy′z′ = τij[1,1]*l[3]*l[2]+τij[2,2]*m[3]*m[2]+τij[3,3]*n[3]*n[2]+τij[1,2]*(l[3]*m[2]+m[3]*l[2])+τij[2,3]*(m[3]*n[2]+m[2]*n[3])+τij[1,3]*(n[3]*l[2]+n[2]*l[3])
end;

# ╔═╡ 1917d5f7-b58f-4f21-ac77-4e0f4624d2d4
md"""
!!! note "Answer"
	```math
	\tau_{ij} = 
	\begin{bmatrix}
		7 & -0.408 & -1.15 \\
		-0.408 & 4.0 & 2.12\\
		-1.15 & 2.12 & 7.0
	\end{bmatrix}\text{MPa}
	```
"""

# ╔═╡ 766a0456-c582-4fc2-967c-09e60248c1dd
md"""
## Subproblem B - Use Matrices to do the Stress Transformation
"""

# ╔═╡ c18ec93d-46b3-4922-9293-f36bd5f91637
let
	R = [0 sqrt(2)/2 -sqrt(2)/2;-sqrt(3)/3 sqrt(3)/3 sqrt(3)/3; sqrt(6)/3 sqrt(6)/6 sqrt(6)/6]
	τij = [4 1 2; 1 6 0; 2 0 8]
	τij′ = R*τij*R'
	latexify(round.(τij′, sigdigits = 3)u"MPa")
end

# ╔═╡ b9a3af42-aaa5-496d-9412-8660bce7d27b
md"""
!!! note "Answer"
	```math
	\tau_{ij} = \begin{bmatrix}
		7 & -0.408 & -1.15 \\
		-0.408 & 4 & 2.12  \\
		-1.15 & 2.12 & 7 
	\end{bmatrix}\text{MPa}
	```
"""

# ╔═╡ 25fb6e02-5cf3-403f-af0b-80a2891cf757
md"""
## Subproblem C - Determine the normal and shear stress normal to the x′ axis
"""

# ╔═╡ 41c2bc2e-16b9-420b-86d6-b88db161fd0f
let
	n⃗ = [0 sqrt(2)/2 -sqrt(2)/2]
	τij = [4 1 2; 1 6 0; 2 0 8]
	σ′ = round(((n⃗*τij*transpose(n⃗)))[1], sigdigits = 3)*1.0u"MPa"
	τ′ = round((n⃗*τij*τij'*transpose(n⃗))[1]-((n⃗*τij*transpose(n⃗)))[1], sigdigits = 3)*1.0u"MPa"
end;

# ╔═╡ ffbf1717-c3ba-46f3-9660-cbfaf1928646
md"""
!!! note "Answer"
	$\sigma^\prime = 7\text{MPa}$
	$\tau^\prime = 43.5\text{MPa}$
"""

# ╔═╡ 82fdb709-c327-4e39-a758-f00b17af7500
md"""
# Problem 7
"""

# ╔═╡ 565759e5-93ea-43be-bed5-cdb81506fa4f
md"""
Given:
```math
\tau_{ij} = \begin{bmatrix}
	150 & -45 & 0\\
	-45 & 70 & 0\\
	0 & 0 & -80
\end{bmatrix}\text{MPa}
```
"""

# ╔═╡ 19c9b13a-f509-472a-9583-5e14ddcc4c54
τ_7 = [150 -45 0; -45 70 0; 0 0 -80];

# ╔═╡ 2d97d47c-ce56-4578-88f0-b20c2404b901
md"""
1. The principal stresses and directions are found by taking the eigenvalues and eigenvectors, respectively, of the stress tensor, $\tau_{ij}$:
"""

# ╔═╡ e795f3b6-1f7b-462c-825f-ff3a212ac053
eigen(τ_7);

# ╔═╡ 63d453a2-b6f2-4df7-b3b3-7b7acfe58b0e
md"""
$\sigma_1 = 170.208\text{MPa},\ \ \ \vec{n} = \langle 0.91224, -0.4096, 0.0 \rangle$
$\sigma_2 = 49.792\text{MPa},\ \ \ \vec{n} = \langle 0.409656, 0.91224, 0.0 \rangle$
$\sigma_3 = -80\text{MPa},\ \ \ \vec{n} = \langle 0, 0, 1 \rangle$
"""

# ╔═╡ beb0b4d5-360e-4dc0-97dc-62075b771feb
md"""
2. The direction between the x and x′ axis can be found with:
```math
	\theta = \arccos{(\langle 1, 0, 0 \rangle \cdot \langle 0.91224, -0.409656, 0.0 \rangle)} = 24.18^\circ
```
"""

# ╔═╡ edf4b202-d55b-48db-ac0f-9b50126ec994
let
	σp_w_vecs = eigen(τ_7)
	max_dir = σp_w_vecs.vectors[:,3]
	θ = acosd(dot([1,0,0], max_dir))
end;

# ╔═╡ fed7dd98-ded9-4266-b4a1-f8dcdf8f1acd
md"""
!!! note "Answer"
	$\sigma_1 = 170.208\text{MPa},\ \ \ \vec{n} = \vec{n} = \langle 0.91224, -0.4096, 0.0 \rangle$
	$\sigma_2 = 49.792\text{MPa},\ \ \ \vec{n} = \langle 0.409656, 0.91224, 0.0 \rangle$
	$\sigma_3 = -80\text{MPa},\ \ \ \vec{n} = \langle 0, 0, 1 \rangle$
	$\theta = 24.18^\circ \text{ between the x-axis and x}^\prime\text{-axis}$
"""

# ╔═╡ a1bfd01f-5b52-4154-a0f8-245d01fd8967
md"""
# Problem 8
"""

# ╔═╡ 9f11b130-82cf-4ba3-b118-fe5d0289b2d2
md"""
Given the stress tensor, $\tau_{ij}$:
```math
\begin{bmatrix}
	3 & 0 & -2\\
	0 & 1 & 2\\
	-2 & 2 & -1
\end{bmatrix}\text{MPa}
```
and the plane:

$x+2y-2z = -4$
The normal and shear stresses are found with:

$\sigma = \vec{n}^T\ \tau_{ij}\ \vec{n}$
$\tau = \sqrt{\vec{n}^T\ \tau{ij}^T\ \tau_{ij} \vec{n} - \sigma}$
where $\vec{n}$ is:

$\langle \frac{1}{\sqrt{1+4+4}}, \frac{2}{\sqrt{1+4+4}}, \frac{-2}{\sqrt{1+4+4}} \rangle = \langle \frac{1}{3}, \frac{2}{3}, -\frac{2}{3}\rangle$ which is determined from the coefficients of the variables in the plane equations.
"""

# ╔═╡ 7bea4c10-0719-432b-9abb-042dbb6bf9cc
let
	τij = [3 0 -2;0 1 2; -2 2 -1]
	coeffs = (x = 1, y = 2, z = -2)
	n = norm(coeffs)
	dir = [coeffs.x/n, coeffs.y/n, coeffs.z/n]
	σobq = round(dir'*τij*dir, sigdigits = 3)
	τobq = round(sqrt(dir'*τij'*τij*dir-σobq^2), sigdigits = 3)
end;

# ╔═╡ a3ade92d-5d57-4823-8448-f7f9749b64c3
md"""
!!! note "Answer"
	$\sigma = -0.556\text{MPa}$
	$\tau = 2.71\text{MPa}$
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
LabelledArrays = "2ee39098-c373-598a-b85f-a56591580800"
Latexify = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
LinearSolve = "7ed4a6bd-45f5-4d41-b270-4a48e9bafcae"
ModelingToolkit = "961ee093-0014-501f-94e3-6117800e7a78"
NonlinearSolve = "8913a72c-1f9b-4ce2-8d82-65094dcecaec"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Symbolics = "0c5d862f-8b57-4792-8d23-62f2024744c7"
Tullio = "bc48ee85-29a4-5162-ae0b-a64e1601d4bc"
Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[compat]
LabelledArrays = "~1.7.0"
Latexify = "~0.15.9"
LinearSolve = "~1.11.0"
ModelingToolkit = "~8.3.2"
NonlinearSolve = "~0.3.14"
Plots = "~1.25.6"
PlutoUI = "~0.7.32"
Symbolics = "~4.3.0"
Tullio = "~0.3.2"
Unitful = "~1.10.1"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.0"
manifest_format = "2.0"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "8eaf9f1b4921132a4cff3f36a1d9ba923b14a481"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.4"

[[deps.AbstractTrees]]
git-tree-sha1 = "03e0550477d86222521d254b741d470ba17ea0b5"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.3.4"

[[deps.Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "af92965fb30777147966f58acb05da51c5616b5f"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.3.3"

[[deps.ArgCheck]]
git-tree-sha1 = "dedbbb2ddb876f899585c4ec4433265e3017215a"
uuid = "dce04be8-c92d-5529-be00-80e4d2c0e197"
version = "2.1.0"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[deps.ArnoldiMethod]]
deps = ["LinearAlgebra", "Random", "StaticArrays"]
git-tree-sha1 = "62e51b39331de8911e4a7ff6f5aaf38a5f4cc0ae"
uuid = "ec485272-7323-5ecc-a04f-4719b315124d"
version = "0.2.0"

[[deps.ArrayInterface]]
deps = ["Compat", "IfElse", "LinearAlgebra", "Requires", "SparseArrays", "Static"]
git-tree-sha1 = "1ee88c4c76caa995a885dc2f22a5d548dfbbc0ba"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "3.2.2"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.AutoHashEquals]]
git-tree-sha1 = "45bb6705d93be619b81451bb2006b7ee5d4e4453"
uuid = "15f4f7f2-30c1-5605-9d31-71845cf9641f"
version = "0.2.0"

[[deps.BangBang]]
deps = ["Compat", "ConstructionBase", "Future", "InitialValues", "LinearAlgebra", "Requires", "Setfield", "Tables", "ZygoteRules"]
git-tree-sha1 = "a33794b483965bf49deaeec110378640609062b1"
uuid = "198e06fe-97b7-11e9-32a5-e1d131e6ad66"
version = "0.3.34"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.Baselet]]
git-tree-sha1 = "aebf55e6d7795e02ca500a689d326ac979aaf89e"
uuid = "9718e550-a3fa-408a-8086-8db961cd8217"
version = "0.1.1"

[[deps.Bijections]]
git-tree-sha1 = "705e7822597b432ebe152baa844b49f8026df090"
uuid = "e2ed5e7c-b2de-5872-ae92-c73ca462fb04"
version = "0.1.3"

[[deps.BitTwiddlingConvenienceFunctions]]
deps = ["Static"]
git-tree-sha1 = "bc1317f71de8dce26ea67fcdf7eccc0d0693b75b"
uuid = "62783981-4cbd-42fc-bca8-16325de8dc4b"
version = "0.1.1"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[deps.CPUSummary]]
deps = ["Hwloc", "IfElse", "Static"]
git-tree-sha1 = "87b0c9c6ee0124d6c1f4ce8cb035dcaf9f90b803"
uuid = "2a0fbf3d-bb9c-48f3-b0a9-814d99fd7ab9"
version = "0.1.6"

[[deps.CSTParser]]
deps = ["Tokenize"]
git-tree-sha1 = "f9a6389348207faf5e5c62cbc7e89d19688d338a"
uuid = "00ebfdb7-1f24-5e51-bd34-a7502290713f"
version = "3.3.0"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "4b859a208b2397a7a623a03449e4636bdb17bcf2"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.1+1"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "54fc4400de6e5c3e27be6047da2ef6ba355511f8"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.11.6"

[[deps.ChangesOfVariables]]
deps = ["ChainRulesCore", "LinearAlgebra", "Test"]
git-tree-sha1 = "bf98fa45a0a4cee295de98d4c1462be26345b9a1"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.2"

[[deps.CloseOpenIntervals]]
deps = ["ArrayInterface", "Static"]
git-tree-sha1 = "7b8f09d58294dc8aa13d91a8544b37c8a1dcbc06"
uuid = "fb6a15b2-703c-40df-9091-08a04967cfa9"
version = "0.1.4"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "Colors", "FixedPointNumbers", "Random"]
git-tree-sha1 = "6b6f04f93710c71550ec7e16b650c1b9a612d0b6"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.16.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "024fe24d83e4a5bf5fc80501a314ce0d1aa35597"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.0"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "417b0ed7b8b838aa6ca0a87aadf1bb9eb111ce40"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.8"

[[deps.Combinatorics]]
git-tree-sha1 = "08c8b6831dc00bfea825826be0bc8336fc369860"
uuid = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
version = "1.0.2"

[[deps.CommonMark]]
deps = ["Crayons", "JSON", "URIs"]
git-tree-sha1 = "4aff51293dbdbd268df314827b7f409ea57f5b70"
uuid = "a80b9123-70ca-4bc0-993e-6e3bcb318db6"
version = "0.8.5"

[[deps.CommonSolve]]
git-tree-sha1 = "68a0743f578349ada8bc911a5cbd5a2ef6ed6d1f"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.0"

[[deps.CommonSubexpressions]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "7b8a93dba8af7e3b42fecabf646260105ac373f7"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.0"

[[deps.Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "44c37b4636bc54afac5c574d2d02b625349d6582"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.41.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[deps.CompositeTypes]]
git-tree-sha1 = "d5b014b216dc891e81fea299638e4c10c657b582"
uuid = "b152e2b5-7a66-4b01-a709-34e65c35f657"
version = "0.1.2"

[[deps.CompositionsBase]]
git-tree-sha1 = "455419f7e328a1a2493cabc6428d79e951349769"
uuid = "a33af91c-f02d-484b-be07-31d278c5ca2b"
version = "0.1.1"

[[deps.ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f74e9d5388b8620b4cee35d4c5a618dd4dc547f4"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.3.0"

[[deps.Contour]]
deps = ["StaticArrays"]
git-tree-sha1 = "9f02045d934dc030edad45944ea80dbd1f0ebea7"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.5.7"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.DEDataArrays]]
deps = ["ArrayInterface", "DocStringExtensions", "LinearAlgebra", "RecursiveArrayTools", "SciMLBase", "StaticArrays"]
git-tree-sha1 = "31186e61936fbbccb41d809ad4338c9f7addf7ae"
uuid = "754358af-613d-5f8d-9788-280bf1605d4c"
version = "0.2.0"

[[deps.DataAPI]]
git-tree-sha1 = "cc70b17275652eb47bc9e5f81635981f13cea5c8"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.9.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "3daef5523dd2e769dad2365274f760ff5f282c7d"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.11"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DefineSingletons]]
git-tree-sha1 = "0fba8b706d0178b4dc7fd44a96a92382c9065c2c"
uuid = "244e2a9f-e319-4986-a169-4d1fe445cd52"
version = "0.1.2"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[deps.DensityInterface]]
deps = ["InverseFunctions", "Test"]
git-tree-sha1 = "80c3e8639e3353e5d2912fb3a1916b8455e2494b"
uuid = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
version = "0.4.0"

[[deps.DiffEqBase]]
deps = ["ArrayInterface", "ChainRulesCore", "DEDataArrays", "DataStructures", "Distributions", "DocStringExtensions", "FastBroadcast", "ForwardDiff", "FunctionWrappers", "IterativeSolvers", "LabelledArrays", "LinearAlgebra", "Logging", "MuladdMacro", "NonlinearSolve", "Parameters", "PreallocationTools", "Printf", "RecursiveArrayTools", "RecursiveFactorization", "Reexport", "Requires", "SciMLBase", "Setfield", "SparseArrays", "StaticArrays", "Statistics", "SuiteSparse", "ZygoteRules"]
git-tree-sha1 = "d75333ab19d6d01c53bb350a9aabb074ba768a9d"
uuid = "2b5f629d-d688-5b77-993f-72d75c75574e"
version = "6.81.3"

[[deps.DiffEqCallbacks]]
deps = ["DataStructures", "DiffEqBase", "ForwardDiff", "LinearAlgebra", "NLsolve", "OrdinaryDiffEq", "Parameters", "RecipesBase", "RecursiveArrayTools", "SciMLBase", "StaticArrays"]
git-tree-sha1 = "e57ecaf9f7875714c164ccca3c802711589127cf"
uuid = "459566f4-90b8-5000-8ac3-15dfb0a30def"
version = "2.20.1"

[[deps.DiffEqJump]]
deps = ["ArrayInterface", "Compat", "DataStructures", "DiffEqBase", "FunctionWrappers", "Graphs", "LinearAlgebra", "PoissonRandom", "Random", "RandomNumbers", "RecursiveArrayTools", "Reexport", "StaticArrays", "TreeViews", "UnPack"]
git-tree-sha1 = "628ddc7e2b44e214232e747b22f1a1d9a8f14467"
uuid = "c894b116-72e5-5b58-be3c-e6d8d4ac2b12"
version = "8.1.0"

[[deps.DiffResults]]
deps = ["StaticArrays"]
git-tree-sha1 = "c18e98cba888c6c25d1c3b048e4b3380ca956805"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.0.3"

[[deps.DiffRules]]
deps = ["LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "9bc5dac3c8b6706b58ad5ce24cffd9861f07c94f"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.9.0"

[[deps.Distances]]
deps = ["LinearAlgebra", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "3258d0659f812acde79e8a74b11f17ac06d0ca04"
uuid = "b4f34e82-e78d-54a5-968a-f98e89d6e8f7"
version = "0.10.7"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.Distributions]]
deps = ["ChainRulesCore", "DensityInterface", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SparseArrays", "SpecialFunctions", "Statistics", "StatsBase", "StatsFuns", "Test"]
git-tree-sha1 = "5863b0b10512ed4add2b5ec07e335dc6121065a5"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.41"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "b19534d1895d702889b219c382a6e18010797f0b"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.8.6"

[[deps.DomainSets]]
deps = ["CompositeTypes", "IntervalSets", "LinearAlgebra", "StaticArrays", "Statistics"]
git-tree-sha1 = "5f5f0b750ac576bcf2ab1d7782959894b304923e"
uuid = "5b8099bc-c8ec-5219-889f-1d9e522a28bf"
version = "0.5.9"

[[deps.Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[deps.DynamicPolynomials]]
deps = ["DataStructures", "Future", "LinearAlgebra", "MultivariatePolynomials", "MutableArithmetics", "Pkg", "Reexport", "Test"]
git-tree-sha1 = "74e63cbb0fda19eb0e69fbe622447f1100cd8690"
uuid = "7c1d4256-1411-5781-91ec-d7bc3513ac07"
version = "0.4.3"

[[deps.EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3f3a2501fa7236e9b911e0f7a588c657e822bb6d"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.3+0"

[[deps.EllipsisNotation]]
deps = ["ArrayInterface"]
git-tree-sha1 = "d7ab55febfd0907b285fbf8dc0c73c0825d9d6aa"
uuid = "da5c29d0-fa7d-589e-88eb-ea29b0a81949"
version = "1.3.0"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b3bfd02e98aedfa5cf885665493c5598c350cd2f"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.2.10+0"

[[deps.ExponentialUtilities]]
deps = ["ArrayInterface", "LinearAlgebra", "Printf", "Requires", "SparseArrays"]
git-tree-sha1 = "3e1289d9a6a54791c1ee60da0850f4fd71188da6"
uuid = "d4d017d3-3776-5f7e-afef-a10c40355c18"
version = "1.11.0"

[[deps.ExprTools]]
git-tree-sha1 = "56559bbef6ca5ea0c0818fa5c90320398a6fbf8d"
uuid = "e2ba6199-217a-4e67-a87a-7c52f15ade04"
version = "0.1.8"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "Pkg", "Zlib_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "d8a578692e3077ac998b50c0217dfd67f21d1e5f"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.0+0"

[[deps.FastBroadcast]]
deps = ["LinearAlgebra", "Polyester", "Static"]
git-tree-sha1 = "0f8ef5dcb040dbb9edd98b1763ac10882ee1ff03"
uuid = "7034ab61-46d4-4ed7-9d0f-46aef9175898"
version = "0.1.12"

[[deps.FastClosures]]
git-tree-sha1 = "acebe244d53ee1b461970f8910c235b259e772ef"
uuid = "9aa1b823-49e4-5ca5-8b0f-3971ec8bab6a"
version = "0.3.2"

[[deps.FillArrays]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "Statistics"]
git-tree-sha1 = "8756f9935b7ccc9064c6eef0bff0ad643df733a3"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "0.12.7"

[[deps.FiniteDiff]]
deps = ["ArrayInterface", "LinearAlgebra", "Requires", "SparseArrays", "StaticArrays"]
git-tree-sha1 = "6eae72e9943d8992d14359c32aed5f892bda1569"
uuid = "6a86dc24-6348-571c-b903-95158fe2bd41"
version = "2.10.0"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "21efd19106a55620a188615da6d3d06cd7f6ee03"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.93+0"

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[deps.ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions", "StaticArrays"]
git-tree-sha1 = "1bd6fc0c344fc0cbee1f42f8d2e7ec8253dda2d2"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.25"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "87eb71354d8ec1a96d4a7636bd57a7347dde3ef9"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.10.4+0"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

[[deps.FunctionWrappers]]
git-tree-sha1 = "241552bc2209f0fa068b6415b1942cc0aa486bcc"
uuid = "069b7b12-0de2-55c6-9aab-29f3d0a68a2e"
version = "1.1.2"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pkg", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll"]
git-tree-sha1 = "0c603255764a1fa0b61752d2bec14cfbd18f7fe8"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.3.5+1"

[[deps.GR]]
deps = ["Base64", "DelimitedFiles", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Pkg", "Printf", "Random", "RelocatableFolders", "Serialization", "Sockets", "Test", "UUIDs"]
git-tree-sha1 = "4a740db447aae0fbeb3ee730de1afbb14ac798a1"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.63.1"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Pkg", "Qt5Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "aa22e1ee9e722f1da183eb33370df4c1aeb6c2cd"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.63.1+0"

[[deps.GeometryBasics]]
deps = ["EarCut_jll", "IterTools", "LinearAlgebra", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "58bcdf5ebc057b085e58d95c138725628dd7453c"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.4.1"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "a32d672ac2c967f3deb8a81d828afc739c838a06"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.68.3+2"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[deps.Graphs]]
deps = ["ArnoldiMethod", "Compat", "DataStructures", "Distributed", "Inflate", "LinearAlgebra", "Random", "SharedArrays", "SimpleTraits", "SparseArrays", "Statistics"]
git-tree-sha1 = "d727758173afef0af878b29ac364a0eca299fc6b"
uuid = "86223c79-3864-5bf0-83f7-82e725a168b6"
version = "1.5.1"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HTTP]]
deps = ["Base64", "Dates", "IniFile", "Logging", "MbedTLS", "NetworkOptions", "Sockets", "URIs"]
git-tree-sha1 = "0fa77022fe4b511826b39c894c90daf5fce3334a"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "0.9.17"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[deps.HostCPUFeatures]]
deps = ["BitTwiddlingConvenienceFunctions", "IfElse", "Libdl", "Static"]
git-tree-sha1 = "8f0dc80088981ab55702b04bba38097a44a1a3a9"
uuid = "3e5b6fbb-0976-4d2c-9146-d79de83f2fb0"
version = "0.1.5"

[[deps.Hwloc]]
deps = ["Hwloc_jll"]
git-tree-sha1 = "92d99146066c5c6888d5a3abc871e6a214388b91"
uuid = "0e44f5e4-bd66-52a0-8798-143a42290a1d"
version = "2.0.0"

[[deps.Hwloc_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "d8bccde6fc8300703673ef9e1383b11403ac1313"
uuid = "e33a78d0-f292-5ffc-b300-72abe9b543c8"
version = "2.7.0+0"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[deps.HypertextLiteral]]
git-tree-sha1 = "2b078b5a615c6c0396c77810d92ee8c6f470d238"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.3"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "f7be53659ab06ddc986428d3a9dcc95f6fa6705a"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.2"

[[deps.IfElse]]
git-tree-sha1 = "debdd00ffef04665ccbb3e150747a77560e8fad1"
uuid = "615f187c-cbe4-4ef1-ba3b-2fcf58d6d173"
version = "0.1.1"

[[deps.Inflate]]
git-tree-sha1 = "f5fc07d4e706b84f72d54eedcc1c13d92fb0871c"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.2"

[[deps.IniFile]]
deps = ["Test"]
git-tree-sha1 = "098e4d2c533924c921f9f9847274f2ad89e018b8"
uuid = "83e8ac13-25f8-5344-8a64-a9f2b223428f"
version = "0.5.0"

[[deps.InitialValues]]
git-tree-sha1 = "4da0f88e9a39111c2fa3add390ab15f3a44f3ca3"
uuid = "22cec73e-a1b8-11e9-2c92-598750a2cf9c"
version = "0.3.1"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.IntervalSets]]
deps = ["Dates", "EllipsisNotation", "Statistics"]
git-tree-sha1 = "3cc368af3f110a767ac786560045dceddfc16758"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.5.3"

[[deps.InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "a7254c0acd8e62f1ac75ad24d5db43f5f19f3c65"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.2"

[[deps.IrrationalConstants]]
git-tree-sha1 = "7fd44fd4ff43fc60815f8e764c0f352b83c49151"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.1"

[[deps.IterTools]]
git-tree-sha1 = "fa6287a4469f5e048d763df38279ee729fbd44e5"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.4.0"

[[deps.IterativeSolvers]]
deps = ["LinearAlgebra", "Printf", "Random", "RecipesBase", "SparseArrays"]
git-tree-sha1 = "1169632f425f79429f245113b775a0e3d121457c"
uuid = "42fd0dbc-a981-5370-80f2-aaf504508153"
version = "0.9.2"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "22df5b96feef82434b07327e2d3c770a9b21e023"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.0"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "8076680b162ada2a031f707ac7b4953e30667a37"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.2"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "d735490ac75c5cb9f1b00d8b5509c11984dc6943"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.0+0"

[[deps.JuliaFormatter]]
deps = ["CSTParser", "CommonMark", "DataStructures", "Pkg", "Tokenize"]
git-tree-sha1 = "da0c8830cebe2337093bb46fc117498517a9df80"
uuid = "98e50ef6-434e-11e9-1051-2b60c6c9e899"
version = "0.21.2"

[[deps.KLU]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse_jll"]
git-tree-sha1 = "1ed18ccf6292d89abf85beba35b9399aeddff9b2"
uuid = "ef3ab10e-7fda-4108-b977-705223b18434"
version = "0.2.3"

[[deps.Krylov]]
deps = ["LinearAlgebra", "Printf", "SparseArrays"]
git-tree-sha1 = "e60270d7871e7ffe66b3a90b477ecb5df037aa0c"
uuid = "ba0b0d4f-ebba-5204-a429-3ac8c609bfb7"
version = "0.7.11"

[[deps.KrylovKit]]
deps = ["LinearAlgebra", "Printf"]
git-tree-sha1 = "0328ad9966ae29ccefb4e1b9bfd8c8867e4360df"
uuid = "0b1a1467-8014-51b9-945f-bf0ae24f4b77"
version = "0.5.3"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6250b16881adf048549549fba48b1161acdac8c"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.1+0"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[deps.LabelledArrays]]
deps = ["ArrayInterface", "ChainRulesCore", "LinearAlgebra", "MacroTools", "StaticArrays"]
git-tree-sha1 = "41158dee1d434944570b02547d404e075da15690"
uuid = "2ee39098-c373-598a-b85f-a56591580800"
version = "1.7.0"

[[deps.Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "Printf", "Requires"]
git-tree-sha1 = "a8f4f279b6fa3c3c4f1adadd78a621b13a506bce"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.15.9"

[[deps.LayoutPointers]]
deps = ["ArrayInterface", "LinearAlgebra", "ManualMemory", "SIMDTypes", "Static"]
git-tree-sha1 = "83b56449c39342a47f3fcdb3bc782bd6d66e1d97"
uuid = "10f19ff3-798f-405d-979b-55457f8fc047"
version = "0.1.4"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[deps.Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "64613c82a59c120435c067c2b809fc61cf5166ae"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.7+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "7739f837d6447403596a75d19ed01fd08d6f56bf"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.3.0+3"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c333716e46366857753e273ce6a69ee0945a6db9"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.42.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "42b62845d70a619f063a7da093d995ec8e15e778"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.16.1+1"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c30530bf0effd46e15e0fdcf2b8636e78cbbd73"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.35.0+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "340e257aada13f95f98ee352d316c3bed37c8ab9"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.3.0+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7f3efec06033682db852f8b3bc3c1d2b0a0ab066"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.36.0+0"

[[deps.LineSearches]]
deps = ["LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "Printf"]
git-tree-sha1 = "f27132e551e959b3667d8c93eae90973225032dd"
uuid = "d3d80556-e9d4-5f37-9878-2ab0fcc64255"
version = "7.1.1"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LinearSolve]]
deps = ["ArrayInterface", "DocStringExtensions", "IterativeSolvers", "KLU", "Krylov", "KrylovKit", "LinearAlgebra", "RecursiveFactorization", "Reexport", "Requires", "SciMLBase", "Setfield", "SparseArrays", "SuiteSparse", "UnPack"]
git-tree-sha1 = "c954090c0a7327a52beccf984610cd505b18d6ce"
uuid = "7ed4a6bd-45f5-4d41-b270-4a48e9bafcae"
version = "1.11.0"

[[deps.LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "e5718a00af0ab9756305a0392832c8952c7426c1"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.6"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.LoopVectorization]]
deps = ["ArrayInterface", "CPUSummary", "ChainRulesCore", "CloseOpenIntervals", "DocStringExtensions", "ForwardDiff", "HostCPUFeatures", "IfElse", "LayoutPointers", "LinearAlgebra", "OffsetArrays", "PolyesterWeave", "SIMDDualNumbers", "SLEEFPirates", "SpecialFunctions", "Static", "ThreadingUtilities", "UnPack", "VectorizationBase"]
git-tree-sha1 = "67c0dfeae307972b50009ce220aae5684ea852d1"
uuid = "bdcacae8-1622-11e9-2a5c-532679323890"
version = "0.12.101"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "3d3e902b31198a27340d0bf00d6ac452866021cf"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.9"

[[deps.ManualMemory]]
git-tree-sha1 = "9cb207b18148b2199db259adfa923b45593fe08e"
uuid = "d125e4d3-2237-4719-b19c-fa641b8a4667"
version = "0.1.6"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "Random", "Sockets"]
git-tree-sha1 = "1c38e51c3d08ef2278062ebceade0e46cefc96fe"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.0.3"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[deps.Measures]]
git-tree-sha1 = "e498ddeee6f9fdb4551ce855a46f54dbd900245f"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.1"

[[deps.Metatheory]]
deps = ["AutoHashEquals", "DataStructures", "Dates", "DocStringExtensions", "Parameters", "Reexport", "TermInterface", "ThreadsX", "TimerOutputs"]
git-tree-sha1 = "0886d229caaa09e9f56bcf1991470bd49758a69f"
uuid = "e9d8d322-4543-424a-9be4-0cc815abe26c"
version = "1.3.3"

[[deps.MicroCollections]]
deps = ["BangBang", "InitialValues", "Setfield"]
git-tree-sha1 = "6bb7786e4f24d44b4e29df03c69add1b63d88f01"
uuid = "128add7d-3638-4c79-886c-908ea0c25c34"
version = "0.1.2"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "bf210ce90b6c9eed32d25dbcae1ebc565df2687f"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.0.2"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.ModelingToolkit]]
deps = ["AbstractTrees", "ArrayInterface", "ConstructionBase", "DataStructures", "DiffEqBase", "DiffEqCallbacks", "DiffEqJump", "DiffRules", "Distributed", "Distributions", "DocStringExtensions", "DomainSets", "Graphs", "IfElse", "InteractiveUtils", "JuliaFormatter", "LabelledArrays", "Latexify", "Libdl", "LinearAlgebra", "MacroTools", "NaNMath", "NonlinearSolve", "RecursiveArrayTools", "Reexport", "Requires", "RuntimeGeneratedFunctions", "SafeTestsets", "SciMLBase", "Serialization", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArrays", "SymbolicUtils", "Symbolics", "UnPack", "Unitful"]
git-tree-sha1 = "0aca0097c38d75d8114844a6c2c048273bbdcfae"
uuid = "961ee093-0014-501f-94e3-6117800e7a78"
version = "8.3.2"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[deps.MuladdMacro]]
git-tree-sha1 = "c6190f9a7fc5d9d5915ab29f2134421b12d24a68"
uuid = "46d2c3a1-f734-5fdb-9937-b9b9aeba4221"
version = "0.2.2"

[[deps.MultivariatePolynomials]]
deps = ["DataStructures", "LinearAlgebra", "MutableArithmetics"]
git-tree-sha1 = "fa6ce8c91445e7cd54de662064090b14b1089a6d"
uuid = "102ac46a-7ee4-5c85-9060-abc95bfdeaa3"
version = "0.4.2"

[[deps.MutableArithmetics]]
deps = ["LinearAlgebra", "SparseArrays", "Test"]
git-tree-sha1 = "73deac2cbae0820f43971fad6c08f6c4f2784ff2"
uuid = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"
version = "0.3.2"

[[deps.NLSolversBase]]
deps = ["DiffResults", "Distributed", "FiniteDiff", "ForwardDiff"]
git-tree-sha1 = "50310f934e55e5ca3912fb941dec199b49ca9b68"
uuid = "d41bc354-129a-5804-8e4c-c37616107c6c"
version = "7.8.2"

[[deps.NLsolve]]
deps = ["Distances", "LineSearches", "LinearAlgebra", "NLSolversBase", "Printf", "Reexport"]
git-tree-sha1 = "019f12e9a1a7880459d0173c182e6a99365d7ac1"
uuid = "2774e3e8-f4cf-5e23-947b-6d7e65073b56"
version = "4.5.1"

[[deps.NaNMath]]
git-tree-sha1 = "f755f36b19a5116bb580de457cda0c140153f283"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "0.3.6"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[deps.NonlinearSolve]]
deps = ["ArrayInterface", "FiniteDiff", "ForwardDiff", "IterativeSolvers", "LinearAlgebra", "RecursiveArrayTools", "RecursiveFactorization", "Reexport", "SciMLBase", "Setfield", "StaticArrays", "UnPack"]
git-tree-sha1 = "b61c51cd5b9d8b197dfcbbf2077a0a4e1505278d"
uuid = "8913a72c-1f9b-4ce2-8d82-65094dcecaec"
version = "0.3.14"

[[deps.OffsetArrays]]
deps = ["Adapt"]
git-tree-sha1 = "043017e0bdeff61cfbb7afeb558ab29536bbb5ed"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.10.8"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "648107615c15d4e09f7eca16307bc821c1f718d8"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.13+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[deps.OrdinaryDiffEq]]
deps = ["Adapt", "ArrayInterface", "DataStructures", "DiffEqBase", "DocStringExtensions", "ExponentialUtilities", "FastClosures", "FiniteDiff", "ForwardDiff", "LinearAlgebra", "LinearSolve", "Logging", "LoopVectorization", "MacroTools", "MuladdMacro", "NLsolve", "NonlinearSolve", "Polyester", "PreallocationTools", "RecursiveArrayTools", "Reexport", "SparseArrays", "SparseDiffTools", "StaticArrays", "UnPack"]
git-tree-sha1 = "1475c25a9dc4de848a5234543f2cb8601ada67d4"
uuid = "1dea7af3-3e70-54e6-95c3-0bf5283fa5ed"
version = "6.6.5"

[[deps.PCRE_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b2a7af664e098055a7529ad1a900ded962bca488"
uuid = "2f80f16e-611a-54ab-bc61-aa92de5b98fc"
version = "8.44.0+0"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "ee26b350276c51697c9c2d88a072b339f9f03d73"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.5"

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Parsers]]
deps = ["Dates"]
git-tree-sha1 = "92f91ba9e5941fc781fecf5494ac1da87bdac775"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.2.0"

[[deps.Pixman_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b4f5d02549a10e20780a24fce72bea96b6329e29"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.40.1+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Requires", "Statistics"]
git-tree-sha1 = "a3a964ce9dc7898193536002a6dd892b1b5a6f1d"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "2.0.1"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "6f1b25e8ea06279b5689263cc538f51331d7ca17"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.1.3"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "GeometryBasics", "JSON", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "PlotThemes", "PlotUtils", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "UUIDs", "UnicodeFun", "Unzip"]
git-tree-sha1 = "db7393a80d0e5bef70f2b518990835541917a544"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.25.6"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "UUIDs"]
git-tree-sha1 = "ae6145ca68947569058866e443df69587acc1806"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.32"

[[deps.PoissonRandom]]
deps = ["Random", "Statistics", "Test"]
git-tree-sha1 = "44d018211a56626288b5d3f8c6497d28c26dc850"
uuid = "e409e4f3-bfea-5376-8464-e040bb5c01ab"
version = "0.4.0"

[[deps.Polyester]]
deps = ["ArrayInterface", "BitTwiddlingConvenienceFunctions", "CPUSummary", "IfElse", "ManualMemory", "PolyesterWeave", "Requires", "Static", "StrideArraysCore", "ThreadingUtilities"]
git-tree-sha1 = "55f5db122f19d8b5b26fe9576edc1ff819e499bb"
uuid = "f517fe37-dbe3-4b94-8317-1923a5111588"
version = "0.6.3"

[[deps.PolyesterWeave]]
deps = ["BitTwiddlingConvenienceFunctions", "CPUSummary", "IfElse", "Static", "ThreadingUtilities"]
git-tree-sha1 = "a3ff99bf561183ee20386aec98ab8f4a12dc724a"
uuid = "1d0040c9-8b98-4ee7-8388-3f51789ca0ad"
version = "0.1.2"

[[deps.PreallocationTools]]
deps = ["Adapt", "ArrayInterface", "ForwardDiff", "LabelledArrays"]
git-tree-sha1 = "e4cb8d4a2edf9b3804c1fb2c2de57d634ff3f36e"
uuid = "d236fae5-4411-538c-8e31-a6e3d9e00b46"
version = "0.2.3"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "2cf929d64681236a2e074ffafb8d568733d2e6af"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.2.3"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.Qt5Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "xkbcommon_jll"]
git-tree-sha1 = "ad368663a5e20dbb8d6dc2fddeefe4dae0781ae8"
uuid = "ea2cea3b-5b76-57ae-a6ef-0a8af62496e1"
version = "5.15.3+0"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "78aadffb3efd2155af139781b8a8df1ef279ea39"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.4.2"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.RandomNumbers]]
deps = ["Random", "Requires"]
git-tree-sha1 = "043da614cc7e95c703498a491e2c21f58a2b8111"
uuid = "e6cf234a-135c-5ec9-84dd-332b85af5143"
version = "1.5.3"

[[deps.RecipesBase]]
git-tree-sha1 = "6bf3f380ff52ce0832ddd3a2a7b9538ed1bcca7d"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.2.1"

[[deps.RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "RecipesBase"]
git-tree-sha1 = "37c1631cb3cc36a535105e6d5557864c82cd8c2b"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.5.0"

[[deps.RecursiveArrayTools]]
deps = ["Adapt", "ArrayInterface", "ChainRulesCore", "DocStringExtensions", "FillArrays", "LinearAlgebra", "RecipesBase", "Requires", "StaticArrays", "Statistics", "ZygoteRules"]
git-tree-sha1 = "5144e1eafb2ecc75765888a4bdcd3a30a6a08b14"
uuid = "731186ca-8d62-57ce-b412-fbd966d074cd"
version = "2.24.1"

[[deps.RecursiveFactorization]]
deps = ["LinearAlgebra", "LoopVectorization", "Polyester", "StrideArraysCore", "TriangularSolve"]
git-tree-sha1 = "832379c5df67f4bab32ed0253ac299cf1e9c36e6"
uuid = "f2c3362d-daeb-58d1-803e-2bc74f2840b4"
version = "0.2.8"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.Referenceables]]
deps = ["Adapt"]
git-tree-sha1 = "e681d3bfa49cd46c3c161505caddf20f0e62aaa9"
uuid = "42d2dcc6-99eb-4e98-b66c-637b7d73030e"
version = "0.1.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "cdbd3b1338c72ce29d9584fdbe9e9b70eeb5adca"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "0.1.3"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "bf3188feca147ce108c76ad82c2792c57abe7b1f"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.7.0"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "68db32dff12bb6127bac73c209881191bf0efbb7"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.3.0+0"

[[deps.RuntimeGeneratedFunctions]]
deps = ["ExprTools", "SHA", "Serialization"]
git-tree-sha1 = "cdc1e4278e91a6ad530770ebb327f9ed83cf10c4"
uuid = "7e49a35a-f44a-4d26-94aa-eba1b4ca6b47"
version = "0.5.3"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[deps.SIMDDualNumbers]]
deps = ["ForwardDiff", "IfElse", "SLEEFPirates", "VectorizationBase"]
git-tree-sha1 = "62c2da6eb66de8bb88081d20528647140d4daa0e"
uuid = "3cdde19b-5bb0-4aaf-8931-af3e248e098b"
version = "0.1.0"

[[deps.SIMDTypes]]
git-tree-sha1 = "330289636fb8107c5f32088d2741e9fd7a061a5c"
uuid = "94e857df-77ce-4151-89e5-788b33177be4"
version = "0.1.0"

[[deps.SLEEFPirates]]
deps = ["IfElse", "Static", "VectorizationBase"]
git-tree-sha1 = "1410aad1c6b35862573c01b96cd1f6dbe3979994"
uuid = "476501e8-09a2-5ece-8869-fb82de89a1fa"
version = "0.6.28"

[[deps.SafeTestsets]]
deps = ["Test"]
git-tree-sha1 = "36ebc5622c82eb9324005cc75e7e2cc51181d181"
uuid = "1bc83da4-3b8d-516f-aca4-4fe02f6d838f"
version = "0.0.1"

[[deps.SciMLBase]]
deps = ["ArrayInterface", "CommonSolve", "ConstructionBase", "Distributed", "DocStringExtensions", "IteratorInterfaceExtensions", "LinearAlgebra", "Logging", "RecipesBase", "RecursiveArrayTools", "StaticArrays", "Statistics", "Tables", "TreeViews"]
git-tree-sha1 = "40c1c606543c0130cd3673f0dd9e11f2b5d76cd0"
uuid = "0bca4576-84f4-4d90-8ffe-ffa030f20462"
version = "1.26.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "0b4b7f1393cff97c33891da2a0bf69c6ed241fda"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.1.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "Requires"]
git-tree-sha1 = "0afd9e6c623e379f593da01f20590bacc26d1d14"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "0.8.1"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "5d7e3f4e11935503d3ecaf7186eac40602e7d231"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.4"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "b3363d7460f7d098ca0912c69b082f75625d7508"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.0.1"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SparseDiffTools]]
deps = ["Adapt", "ArrayInterface", "Compat", "DataStructures", "FiniteDiff", "ForwardDiff", "Graphs", "LinearAlgebra", "Requires", "SparseArrays", "StaticArrays", "VertexSafeGraphs"]
git-tree-sha1 = "75c89362201983c500dd34923b015dbecdae7a90"
uuid = "47a9eef4-7e08-11e9-0b38-333d64bd3804"
version = "1.20.0"

[[deps.SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "e08890d19787ec25029113e88c34ec20cac1c91e"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.0.0"

[[deps.SplittablesBase]]
deps = ["Setfield", "Test"]
git-tree-sha1 = "39c9f91521de844bad65049efd4f9223e7ed43f9"
uuid = "171d559e-b47b-412a-8079-5efa626c420e"
version = "0.1.14"

[[deps.Static]]
deps = ["IfElse"]
git-tree-sha1 = "7f5a513baec6f122401abfc8e9c074fdac54f6c1"
uuid = "aedffcd0-7271-4cad-89d0-dc628f76c6d3"
version = "0.4.1"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "2884859916598f974858ff01df7dfc6c708dd895"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.3.3"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StatsAPI]]
git-tree-sha1 = "d88665adc9bcf45903013af0982e2fd05ae3d0a6"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.2.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "51383f2d367eb3b444c961d485c565e4c0cf4ba0"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.14"

[[deps.StatsFuns]]
deps = ["ChainRulesCore", "InverseFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "f35e1879a71cca95f4826a14cdbf0b9e253ed918"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "0.9.15"

[[deps.StrideArraysCore]]
deps = ["ArrayInterface", "CloseOpenIntervals", "IfElse", "LayoutPointers", "ManualMemory", "Requires", "SIMDTypes", "Static", "ThreadingUtilities"]
git-tree-sha1 = "12cf3253ebd8e2a3214ae171fbfe51e7e8d8ad28"
uuid = "7792a7ef-975c-4747-a70f-980b88e8d1da"
version = "0.2.9"

[[deps.StructArrays]]
deps = ["Adapt", "DataAPI", "StaticArrays", "Tables"]
git-tree-sha1 = "d21f2c564b21a202f4677c0fba5b5ee431058544"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.4"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "Pkg", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"

[[deps.SymbolicUtils]]
deps = ["AbstractTrees", "Bijections", "ChainRulesCore", "Combinatorics", "ConstructionBase", "DataStructures", "DocStringExtensions", "DynamicPolynomials", "IfElse", "LabelledArrays", "LinearAlgebra", "Metatheory", "MultivariatePolynomials", "NaNMath", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArrays", "TermInterface", "TimerOutputs"]
git-tree-sha1 = "bfa211c9543f8c062143f2a48e5bcbb226fd790b"
uuid = "d1185830-fcd6-423d-90d6-eec64667417b"
version = "0.19.7"

[[deps.Symbolics]]
deps = ["ArrayInterface", "ConstructionBase", "DataStructures", "DiffRules", "Distributions", "DocStringExtensions", "DomainSets", "IfElse", "Latexify", "Libdl", "LinearAlgebra", "MacroTools", "Metatheory", "NaNMath", "RecipesBase", "Reexport", "Requires", "RuntimeGeneratedFunctions", "SciMLBase", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArrays", "SymbolicUtils", "TermInterface", "TreeViews"]
git-tree-sha1 = "074e08aea1c745664da5c4b266f50b840e528b1c"
uuid = "0c5d862f-8b57-4792-8d23-62f2024744c7"
version = "4.3.0"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "TableTraits", "Test"]
git-tree-sha1 = "bb1064c9a84c52e277f1096cf41434b675cd368b"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.6.1"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[deps.TermInterface]]
git-tree-sha1 = "7aa601f12708243987b88d1b453541a75e3d8c7a"
uuid = "8ea1fca8-c5ef-4a55-8b96-4e9afe9c9a3c"
version = "0.2.3"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.ThreadingUtilities]]
deps = ["ManualMemory"]
git-tree-sha1 = "884539ba8c4584a3a8173cb4ee7b61049955b79c"
uuid = "8290d209-cae3-49c0-8002-c8c24d57dab5"
version = "0.4.7"

[[deps.ThreadsX]]
deps = ["ArgCheck", "BangBang", "ConstructionBase", "InitialValues", "MicroCollections", "Referenceables", "Setfield", "SplittablesBase", "Transducers"]
git-tree-sha1 = "6dad289fe5fc1d8e907fa855135f85fb03c8fa7a"
uuid = "ac1d9e8a-700a-412c-b207-f0111f4b6c0d"
version = "0.1.9"

[[deps.TimerOutputs]]
deps = ["ExprTools", "Printf"]
git-tree-sha1 = "97e999be94a7147d0609d0b9fc9feca4bf24d76b"
uuid = "a759f4b9-e2f1-59dc-863e-4aeb61b1ea8f"
version = "0.5.15"

[[deps.Tokenize]]
git-tree-sha1 = "0952c9cee34988092d73a5708780b3917166a0dd"
uuid = "0796e94c-ce3b-5d07-9a54-7f471281c624"
version = "0.5.21"

[[deps.Transducers]]
deps = ["Adapt", "ArgCheck", "BangBang", "Baselet", "CompositionsBase", "DefineSingletons", "Distributed", "InitialValues", "Logging", "Markdown", "MicroCollections", "Requires", "Setfield", "SplittablesBase", "Tables"]
git-tree-sha1 = "3f0945b47207a41946baee6d1385e4ca738c25f7"
uuid = "28d57a85-8fef-5791-bfe6-a80928e7c999"
version = "0.4.68"

[[deps.TreeViews]]
deps = ["Test"]
git-tree-sha1 = "8d0d7a3fe2f30d6a7f833a5f19f7c7a5b396eae6"
uuid = "a2a6695c-b41b-5b7d-aed9-dbfdeacea5d7"
version = "0.3.0"

[[deps.TriangularSolve]]
deps = ["CloseOpenIntervals", "IfElse", "LayoutPointers", "LinearAlgebra", "LoopVectorization", "Polyester", "Static", "VectorizationBase"]
git-tree-sha1 = "c3ab8b77b82fd92e2b6eea8a275a794d5a6e4011"
uuid = "d5829a12-d9aa-46ab-831f-fb7c9ab06edf"
version = "0.1.9"

[[deps.Tullio]]
deps = ["ChainRulesCore", "DiffRules", "LinearAlgebra", "Requires"]
git-tree-sha1 = "0288b7a395fc412952baf756fac94e4f28bfec65"
uuid = "bc48ee85-29a4-5162-ae0b-a64e1601d4bc"
version = "0.3.2"

[[deps.URIs]]
git-tree-sha1 = "97bbe755a53fe859669cd907f2d96aee8d2c1355"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.3.0"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unitful]]
deps = ["ConstructionBase", "Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "b95e0b8a8d1b6a6c3e0b3ca393a7a285af47c264"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.10.1"

[[deps.Unzip]]
git-tree-sha1 = "34db80951901073501137bdbc3d5a8e7bbd06670"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.1.2"

[[deps.VectorizationBase]]
deps = ["ArrayInterface", "CPUSummary", "HostCPUFeatures", "Hwloc", "IfElse", "LayoutPointers", "Libdl", "LinearAlgebra", "SIMDTypes", "Static"]
git-tree-sha1 = "6e261bff5c9f2537776165dea3067df9de4440cf"
uuid = "3d5dd08c-fd9d-11e8-17fa-ed2836048c2f"
version = "0.21.23"

[[deps.VertexSafeGraphs]]
deps = ["Graphs"]
git-tree-sha1 = "8351f8d73d7e880bfc042a8b6922684ebeafb35c"
uuid = "19fa3120-7c27-5ec5-8db8-b0b0aa330d6f"
version = "0.2.0"

[[deps.Wayland_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "3e61f0b86f90dacb0bc0e73a0c5a83f6a8636e23"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.19.0+0"

[[deps.Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "66d72dc6fcc86352f01676e8f0f698562e60510f"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.23.0+0"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "1acf5bdf07aa0907e0a37d3718bb88d4b687b74a"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.9.12+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "5be649d550f3f4b95308bf0183b82e2582876527"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.6.9+4"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4e490d5c960c314f33885790ed410ff3a94ce67e"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.9+4"

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fe47bd2247248125c428978740e18a681372dd4"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.3+4"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "b7c0aa8c376b31e4852b360222848637f481f8c3"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.4+4"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[deps.Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[deps.Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "19560f30fd49f4d4efbe7002a1037f8c43d43b96"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.10+4"

[[deps.Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6783737e45d3c59a4a4c4091f5f88cdcf0908cbb"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.0+3"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "daf17f441228e7a3833846cd048892861cff16d6"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.13.0+3"

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "926af861744212db0eb001d9e40b5d16292080b2"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.0+4"

[[deps.Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[deps.Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[deps.Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "4bcbf660f6c2e714f87e960a171b119d06ee163b"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.2+4"

[[deps.Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "5c8424f8a67c3f2209646d4425f3d415fee5931d"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.27.0+4"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "79c31e7844f6ecf779705fbc12146eb190b7d845"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.4.0+3"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "cc4bf3fdde8b7e3e9fa0351bdeedba1cf3b7f6e6"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.0+0"

[[deps.ZygoteRules]]
deps = ["MacroTools"]
git-tree-sha1 = "8c1a8e4dfacb1fd631745552c8db35d0deb09ea0"
uuid = "700de1a5-db45-46bc-99cf-38207098b444"
version = "0.2.2"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "94d180a6d2b5e55e447e2d27a29ed04fe79eb30c"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.38+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "b910cb81ef3fe6e78bf6acee440bda86fd6ae00c"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+1"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"

[[deps.xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "ece2350174195bb31de1a63bea3a41ae1aa593b6"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "0.9.1+5"
"""

# ╔═╡ Cell order:
# ╟─d3a088a2-6495-4750-85f8-a01f9092bf6e
# ╟─d018bb86-7e27-11ec-3db2-cd579d41ad35
# ╠═1ee74d1b-f96a-41e7-9c0b-d552cfc3f2c6
# ╟─79a774e2-14b4-4079-a673-404d7402a512
# ╟─a7f73c0c-8af7-4679-a42a-49652a554ff1
# ╟─bef16de0-e39b-425c-b0e4-37434730b48e
# ╟─65d33acc-6b38-4f33-8da9-085817f34bea
# ╟─4cafc5a2-672f-47c3-b11a-445ec7a757e6
# ╟─5393419e-026d-425a-8f93-5ea593202120
# ╟─3fd8a785-96a3-451f-b383-6c87555eb1dc
# ╟─7de1b874-df71-43c9-849d-e421b3fec528
# ╟─80a4c591-97d5-4f0e-918a-c9520e9d0f3b
# ╟─7718d209-08bd-455b-a2df-b4deea805b61
# ╟─86464ad1-ef40-4566-b1cb-5c4d7834107f
# ╟─f4c74e07-a30b-4c3e-9f6f-30aa05e1e69e
# ╟─6eed76e7-24bf-4ca9-bf17-a83b9bc755de
# ╟─130ecbe0-1cfe-4547-84f5-094bf8044126
# ╟─dbce3d98-51e3-456c-a9b2-1457ba62c700
# ╟─1e3cd5d8-3075-4746-8637-c912d616ece2
# ╟─93f83718-2ba3-4988-b5e1-4ef0b6d92039
# ╟─b97ee3f5-a198-4ba2-a087-70ef315db537
# ╟─db3b3cbb-39a6-42e0-921a-bfa53e805eea
# ╟─54520742-17a2-4a66-b7f5-feda3852745e
# ╟─b76a4896-f87e-4176-a28f-d9fbda2b7f19
# ╟─f88ce048-e1c1-47c0-900d-bd81739dfd67
# ╟─3ccbb53c-80ab-45b2-8244-e96c225d4494
# ╟─fc9dfb09-fb8e-4683-b0dc-4c1375744126
# ╟─6e09578e-e23f-4cfb-a813-10c5aa8be359
# ╟─bc07861d-b324-4a1c-bf77-c1290508baf4
# ╟─8de4abab-f5e1-432e-902f-162169468d73
# ╟─c2fa4f8f-3689-4a5d-87ee-d182ef250f4a
# ╟─0651538b-cf6e-4f92-87b2-e4ec940e1aa3
# ╟─583a0ab9-db59-43be-8583-d97d39a8c486
# ╟─51ab8601-c9c5-44b9-b041-0055731e797a
# ╟─617ac012-e2d0-4708-af65-8cd21c9d45ee
# ╟─3ca25a1c-6753-4015-976f-77ccc806eba8
# ╟─08674b25-909b-42ac-a487-0f5c56e1616b
# ╟─bff27eec-f0eb-4cce-b4d6-70fd7b4ec549
# ╟─62ba9781-a790-4149-a399-0e06ad9e52d5
# ╠═493608ce-70ee-468f-a2d4-f8f443adf7df
# ╠═decad7b8-75d7-4a01-884b-bfafab73cf32
# ╟─3adc9eab-9642-494a-9ea2-0c0edd79935e
# ╟─66d9df69-7996-4493-a516-317ac6f613d9
# ╟─be145b51-6cb6-4628-8447-b8d17a703d96
# ╟─b0c491d5-30aa-4c85-ab41-9016e8cc64de
# ╟─a59e79c6-a86f-4d60-b38a-090e99aae5e6
# ╟─211ed5cf-6338-4b6c-ae9e-ffb3fd1d3a50
# ╟─230b2b41-802a-46d2-9d92-62ac9922c43d
# ╟─ee92bd0b-29cc-4fe7-a132-78568743fcb5
# ╟─fe20af7d-05a4-450d-b7ab-7003055586df
# ╟─9729b10d-6c6d-465e-aff1-38eaf7fca34c
# ╟─cfb974e3-e281-478d-907b-15ec708af696
# ╟─23ecb4ec-92a5-465c-b425-058753c9b6f2
# ╠═6c99da8e-c83a-4552-8e3e-1710e4932bba
# ╟─224242da-8ff5-4485-a8aa-2b15820f8841
# ╠═93417d49-684d-4bb5-b220-c95dec571035
# ╟─6533e757-73f7-4f3c-a20e-a9e90b9943b7
# ╟─3375015f-a720-4eb2-9c67-f133da64da54
# ╠═34c79484-8cf6-48ef-a08e-c0a0e54f32e9
# ╟─fee104ad-6f7d-40d7-81e9-dd665e6532b5
# ╟─9df061be-bda1-45e3-a8c0-062b7c0b7e4b
# ╠═913f2bf9-0d2a-49ec-b6c9-3b87684f595c
# ╟─93a8de29-e4c5-4f2b-acb8-8b66712a45d0
# ╟─22592932-d897-4ad4-8549-a322ea2eb1c8
# ╟─4f3391ca-1e2d-4f78-85b2-45d07ab6db6d
# ╠═7b8d5f1f-5b28-40da-9e30-0e060e1a2d8b
# ╟─1917d5f7-b58f-4f21-ac77-4e0f4624d2d4
# ╟─766a0456-c582-4fc2-967c-09e60248c1dd
# ╠═c18ec93d-46b3-4922-9293-f36bd5f91637
# ╟─b9a3af42-aaa5-496d-9412-8660bce7d27b
# ╟─25fb6e02-5cf3-403f-af0b-80a2891cf757
# ╠═41c2bc2e-16b9-420b-86d6-b88db161fd0f
# ╟─ffbf1717-c3ba-46f3-9660-cbfaf1928646
# ╟─82fdb709-c327-4e39-a758-f00b17af7500
# ╟─565759e5-93ea-43be-bed5-cdb81506fa4f
# ╟─19c9b13a-f509-472a-9583-5e14ddcc4c54
# ╟─2d97d47c-ce56-4578-88f0-b20c2404b901
# ╠═e795f3b6-1f7b-462c-825f-ff3a212ac053
# ╟─63d453a2-b6f2-4df7-b3b3-7b7acfe58b0e
# ╟─beb0b4d5-360e-4dc0-97dc-62075b771feb
# ╠═edf4b202-d55b-48db-ac0f-9b50126ec994
# ╟─fed7dd98-ded9-4266-b4a1-f8dcdf8f1acd
# ╟─a1bfd01f-5b52-4154-a0f8-245d01fd8967
# ╟─9f11b130-82cf-4ba3-b118-fe5d0289b2d2
# ╠═7bea4c10-0719-432b-9abb-042dbb6bf9cc
# ╟─a3ade92d-5d57-4823-8448-f7f9749b64c3
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
