### A Pluto.jl notebook ###
# v0.19.5

using Markdown
using InteractiveUtils

# ╔═╡ 86e4051d-ede6-4e1d-9676-98d92140cc78
using BenchmarkTools

# ╔═╡ 27793eba-f876-4e92-849d-945989e39216
using ForwardDiff

# ╔═╡ ac232fbf-827e-4363-9eac-04885334368a
using JET

# ╔═╡ 144e4f59-c0d0-4f75-aa49-d05cf781d83e
using LoopVectorization

# ╔═╡ 6b3592a6-43d0-4ebd-a06d-d325dd66a829
using PlutoUI

# ╔═╡ 06e6e4b8-31b2-4f9e-9da3-093a030021af
using ProfileCanvas

# ╔═╡ 22d3e9db-f4ec-4810-b2d5-94bdbf8d52ac
using ProgressLogging

# ╔═╡ 155a8c77-2837-44ee-9b7a-68d5a89a2fdc
md"""
!!! danger "High-performance Julia"

	_by Guillaume Dalle_ ([website](https://gdalle.github.io/))

	This tutorial introduces various tools to analyze and speed up Julia programs.
"""

# ╔═╡ c48f5b63-48b0-4bb6-997d-7d40745a46fa
md"""
To experiment with it, click on `Edit or run this notebook`. You will then have a choice between running the code remotely (on Binder) or locally (after installing Julia and Pluto).

If you're used to Jupyter, you may find Pluto surprising: the cells are not executed in the order in which they appear, but in the _topological order of variable dependencies_. Besides, the output (last expression) of a cell is displayed _above it_ and not below. Read the [FAQ](https://github.com/fonsp/Pluto.jl/wiki) to know more.
"""

# ╔═╡ 5de2a556-f3af-4a64-a5c6-32d30f758be3
TableOfContents()

# ╔═╡ 1ac5ba38-0eef-41bb-8f9c-3bbf057cae21
VERSION

# ╔═╡ e1852c8d-4028-409e-8e1a-8253bbd6e6a5
html"<button onclick='present()'>Toggle presentation mode</button>"

# ╔═╡ 9331fad2-f29e-11eb-0349-477bd2e7e412
md"""
# Analysis toolbox

Before trying to improve the efficiency of our code, it is essential to analyze it and locate potential improvements.
"""

# ╔═╡ 3d98e7db-c643-4500-987d-4a225e55b2a5
md"""
## Tracking loops

In long-running code, the best way to track loops is not a periodic `println(i)`. There are packages designed for this purpose, such as [ProgressMeter.jl](https://github.com/timholy/ProgressMeter.jl).
However, since the REPL doesn't work well in Pluto notebooks, we can use the `@progress` macro of [ProgressLogging.jl](https://github.com/JuliaLogging/ProgressLogging.jl) instead.
"""

# ╔═╡ b4f2a99e-de45-49d2-be86-9f2d03357462
@progress for i in 1:10
	sleep(0.2)
end

# ╔═╡ 068b3e45-5105-48aa-a547-536470f6abda
md"""
Julia also has a built-in [logging system](https://docs.julialang.org/en/v1/stdlib/Logging/) which Pluto natively understands.
"""

# ╔═╡ 7a75c990-46b4-484e-acc4-65e34f41a9f2
for i in 1:9
	sleep(0.2)
	if i % 3 == 1
		@info "Good integer" i i%3
	elseif i % 3 == 2
		@warn "Okay integer" i i%3
	else
		@error "Bad integer" i i%3
	end
end

# ╔═╡ f7b1b44f-2aa6-4c5c-97a2-ac7037fb48ce
md"""
## Benchmarking

To evaluate the efficiency of a function, we need to know how long it takes and how much memory it uses. Julia provides built-in macros for these tasks:
- `@elapsed` returns the computation time (in seconds)
- `@allocated` returns the allocated memory (in bytes)
- `@time` prints both (in the REPL!) and returns the function result
"""

# ╔═╡ 1fb43343-083b-4b1a-b622-d88c9aa0808c
@elapsed exp(rand(100, 100))

# ╔═╡ a28f7911-3dbb-45fb-a82d-2834d3c8502c
@allocated exp(rand(100, 100))

# ╔═╡ 4da8a7ca-3cea-4629-a66d-44f3b907af09
@time exp(rand(100, 100));

# ╔═╡ c0a7c1fe-457f-4e52-b0ea-2821e40817ea
md"""
However, the built-in macros have shortcomings: they only run the function once, and their measurements may be biased by the presence of global variables.
We can get a more accurate evaluation thanks to [BenchmarkTools.jl](https://github.com/JuliaCI/BenchmarkTools.jl):

- `@belapsed` is similar to `@elapsed`
- `@ballocated` is similar to `@allocated`
- `@btime` is similar to `@time`
- `@benchmark` prints a pretty graphical summary of `@btime`
"""

# ╔═╡ 4e73e893-31c5-4c34-82da-ddabaa7316a4
@benchmark exp(rand(100, 100))

# ╔═╡ ddbe8924-756a-46a8-9af1-6df0a21e9057
md"""
When using BenchmarkTools.jl, it is important to [interpolate](https://juliaci.github.io/BenchmarkTools.jl/stable/manual/#Interpolating-values-into-benchmark-expressions) any external (especially global) variables with a dollar sign, so that they're evaluated at definition time (see below).
"""

# ╔═╡ 94c78148-c651-4a59-9e62-5c7e9576d1e8
md"""
## Profiling

Profiling is more precise than benchmarking: it tells you how much time you spend _in each nested function call_.
Julia provides a sampling-based [profiler](https://docs.julialang.org/en/v1/manual/profile/), but it is hard to use without a good visualization.
We recommend using the [`@profview` macro](https://www.julia-vscode.org/docs/stable/userguide/profiler/) from the VSCode extension (or from [ProfileCanvas.jl](https://github.com/pfitzseb/ProfileCanvas.jl) when you're working on a notebook).
However, there are many other options: see the [FlameGraphs.jl](https://github.com/timholy/FlameGraphs.jl) README for a list.

A profiling "flame graph" represents the entire call stack (without C routines), with each layer corresponding to a call depth.
The width of a tile is proportional to its execution time, but you can click on it to make it fill the whole window.
Note that the first few layers are usually boilerplate code, and we need to scroll down to reach user-defined functions, usually below a tile called `eval`.
"""

# ╔═╡ c44f3dc9-ff19-4ba4-9388-73cfaf23f8e8
@profview exp(rand(1000, 1000))

# ╔═╡ a7de0ec9-6b01-4b42-8cce-bb2295da779f
md"""
The colors in the flame graph have special meanings:
- blue $\implies$ everything is fine
- gray $\implies$ compilation overhead from the first function call (just run the profiling step a second time)
- red $\implies$ "runtime dispatch" flag, a sign of bad type inference (except in the first layers where it is normal)
- yellow $\implies$ "garbage collection" (GC) flag, a sign of excessive allocations

The importance of these warnings suggests two basic principles for high-performance Julia:

!!! warning "How to write efficient Julia code"
	1. Facilitate type inference
	2. Reduce memory allocations

We will come back to them in the next section.
"""

# ╔═╡ 9ca598c1-dae8-40b9-a18d-c74f30524b35
md"""
## Diagnosing

**Type inference**

The built-in macro [`@code_warntype`](https://docs.julialang.org/en/v1/manual/performance-tips/#man-code-warntype) shows the result of type inference on a function call.
Non-concrete types are displayed in red: they are those for which inference failed.

Sometimes `@code_warntype` is not enough, because it only studies the outermost function and doesn't dive deeper into the call stack. This is what the macro `@report_opt` from [JET.jl](https://github.com/aviatesk/JET.jl) is for.

**Allocations**

Julia 1.8 (which is still a pre-release) introduced many novelties.
Among them is a built-in memory profiler, which mimics the behavior of the temporal profiler shown above. We will not demonstrate it here.
"""

# ╔═╡ 0fb6ed33-601c-4392-b7d9-32230c979d39
md"""
# Performance tips

Now that we know how to detect bad performance, we will discover how to achieve good performance! But always remember the golden rule: only optimize what needs optimizing! In other words,

> Premature optimization is the root of all evil. (Donald Knuth)

The primary source for this section is the [Julia manual page on performance tips](https://docs.julialang.org/en/v1/manual/performance-tips/) (read it!).
"""

# ╔═╡ a6e9da76-1ff0-4b54-9b55-4856ca32b251
md"""
## General advice

- Avoid [global variables](https://docs.julialang.org/en/v1/manual/performance-tips/#Avoid-global-variables), or turn them into constants with the keyword `const`
- Put critical code [inside functions](https://docs.julialang.org/en/v1/manual/performance-tips/#Performance-critical-code-should-be-inside-a-function)
- Vectorized operations (using the [dot syntax](https://docs.julialang.org/en/v1/manual/functions/#man-vectorized)) are not faster than loops, except linear algebra routines
- Beware of [closures](https://docs.julialang.org/en/v1/manual/performance-tips/#man-performance-captured) (i.e. functions that return functions)
"""

# ╔═╡ d3c1a86c-8c8f-4ad6-ac3c-2ba0f838d139
md"""
## Facilitate type inference

Julia is fastest when it can infer the type of each variable during just-in-time compilation: then it can decide ahead of runtime (statically) which method to dispatch where.
When this fails, types have to be inferred at runtime (dynamically), and "runtime dispatch" of methods is much slower.

!!! note "The key to successful type inference"
	In each function, the types of the arguments (*not their values*) should suffice to deduce the type of every other variable, especially the output.

Here are a few ways to make this happen.

- Always declare concrete or parametric types (no abstract types) in the following places:
  - [container initializations](https://docs.julialang.org/en/v1/manual/performance-tips/#man-performance-abstract-container)
  - [`struct` field values](https://docs.julialang.org/en/v1/manual/performance-tips/#Avoid-fields-with-abstract-type)
  - [`struct` field containers](https://docs.julialang.org/en/v1/manual/performance-tips/#Avoid-fields-with-abstract-containers)
- Never write `if typeof(x) == ...`: exploit [multiple dispatch](https://docs.julialang.org/en/v1/manual/performance-tips/#Break-functions-into-multiple-definitions) or [function barriers](https://docs.julialang.org/en/v1/manual/performance-tips/#kernel-functions) instead
- Define functions that [do not change the type of variables](https://docs.julialang.org/en/v1/manual/performance-tips/#Avoid-changing-the-type-of-a-variable) and [always output the same type](https://docs.julialang.org/en/v1/manual/performance-tips/#Write-%22type-stable%22-functions)
"""

# ╔═╡ fa483fea-bf9f-4764-8d4f-c6d33e3336fb
md"""
## Reduce memory allocations

Allocations and garbage collection are significant performance bottlenecks. Here are some ways to avoid them:

- Prefer in-place functions that reuse available containers (they name usually [ends with `!`](https://docs.julialang.org/en/v1/manual/style-guide/#bang-convention))
- [Pre-allocate](https://docs.julialang.org/en/v1/manual/performance-tips/#Pre-allocating-outputs) output memory
- Use [views instead of slices](https://docs.julialang.org/en/v1/manual/performance-tips/#man-performance-views) when you don't need copies: `view(A, :, 1)` instead of `A[:, 1]`
- [Combine vectorized operations](https://docs.julialang.org/en/v1/manual/performance-tips/#More-dots:-Fuse-vectorized-operations)
- Fix type inference bugs (they often lead to increased memory use)
"""

# ╔═╡ 00e96665-ccbb-42c1-b5de-f626d57bb0ff
md"""
# Examples
"""

# ╔═╡ 43bad028-9d16-426f-9cdb-a37b1ee1a623
md"""
## Matrix multiplication

Our goal here is to write a function that stores the product of two matrices $A \in \mathbb{R}^{m \times n}$ and $B \in \mathbb{R}^{n \times p}$ within the matrix $C \in \mathbb{R}^{m \times p}$, without using built-in linear algebra.
Starting from a simple implementation, we will enhance it step by step until we are satisfied.
"""

# ╔═╡ 857509a7-f07a-4bf0-9383-207984b95faa
A, B = rand(100, 300), rand(300, 200);

# ╔═╡ 7386749b-b2ab-48a7-a1d2-46e7f31e72e3
md"""
**Version 1**

Our first attempt aims at the simplest possible correct code. A programmer used to Python or R may try something like this.
"""

# ╔═╡ 6cd86e7a-9f82-4da1-a8f0-4ed2c1068ab9
function matmul1(A, B)
	@assert size(A, 2) == size(B, 1)
	m, p = size(A, 1), size(B, 2)
	C = zeros(m, p)
	for i = 1:m, j = 1:p
		C[i, j] = sum(A[i, :] .* B[:, j])
	end
	return C
end

# ╔═╡ 5d10a00b-bfa9-49c7-9f4b-503351fa2842
@btime matmul1($A, $B);

# ╔═╡ 38cc6383-c7d8-46b4-8531-251bd196d960
md"""
Is that a good running time? Can we do better? Hard to tell from the benchmark alone, which is why we need profiling.
"""

# ╔═╡ 62163e17-4265-4c97-95bb-29d608e80b07
@profview matmul1(A, B)

# ╔═╡ 9d8b7e25-c9c6-4aba-a33c-66fd18d804c0
md"""
Scrolling down the flame graph, we see that a lot of time is spent in the `Array` constructor, which triggers garbage collection.
This is caused by unnecessary allocations: when we do
```julia
C[i,j] = sum(A[i, :] .* B[:, j])
```
we create one copy of the row `A[i, :]`, one copy of the column `B[:, j]`, and a whole new vector to store their componentwise product before summing over it. Such a waste of memory!
"""

# ╔═╡ fe45168c-8cf1-435e-86fc-16cfffef3ec1
md"""
**Version 2**

Our second version remedies this problem by computing the dot product manually.
"""

# ╔═╡ 0400175c-5a3c-44a7-9a8a-c30a4756b88c
function matmul2(A, B)
	@assert size(A, 2) == size(B, 1)
	m, n, p = size(A, 1), size(A, 2), size(B, 2)
	C = zeros(m, p)
	for i = 1:m, j = 1:p, k = 1:n
		C[i, j] += A[i, k] * B[k, j]
	end
	return C
end

# ╔═╡ cd0cc22f-2d4d-4848-8f15-8f0127a4245b
@btime matmul2($A, $B);

# ╔═╡ 638de554-1bec-453d-9e30-796247aaa4cc
md"""
The running time has decreased a lot, and so has the allocated memory.
"""

# ╔═╡ fd4401cf-69e8-4444-92c3-478035301006
@profview matmul2(A, B)

# ╔═╡ 23053665-e058-43de-95d9-c688e3a80b0c
md"""
This time we see that the main bottlenecks are `setindex!` and `getindex` (which are used to access and modify components of an array), along with `+` (addition) and `iterate` (from the `for` loop). This is much more coherent, and if we want to do better, we will need to roll up our sleeves.
"""

# ╔═╡ 9a181530-02e7-47b0-9a86-c191baefac54
md"""
**Version 3**

Our third version uses some dark magic from [LoopVectorization.jl](https://github.com/JuliaSIMD/LoopVectorization.jl) to tell the compiler that operations can be vectorized (in a hardware sense, see the "Going further" section). We show this for the sake of completeness, but you should not use the `@turbo` macro unless you really know what you are doing (read the docs first!).
"""

# ╔═╡ c171555a-0166-476e-8ec6-1860745d84f2
function matmul3(A, B)
	@assert size(A, 2) == size(B, 1)
	m, n, p = size(A, 1), size(A, 2), size(B, 2)
	C = zeros(m, p)
	@turbo for i = 1:m, j = 1:p, k = 1:n
		C[i, j] += A[i, k] * B[k, j]
	end
	return C
end

# ╔═╡ 610f6d6f-9d37-4f3d-be78-ab9847162f4d
@btime matmul3($A, $B);

# ╔═╡ c7b551a0-8c2e-4785-b575-8d58e37c14ec
md"""
What a mind-blowing performance gain! This time we seem to have reached the speed of the built-in matrix product, so it is a good time to stop.
"""

# ╔═╡ 0d8af577-9275-490c-a689-65e7177c4d65
@btime $A * $B;

# ╔═╡ 69e8bf4e-d98d-4804-b6bf-f299c3452565
md"""
## Point storage

Our goal here is to compare different ways to store a point in 2D space.
We first define an abstract type (or interface).
"""

# ╔═╡ dacdb662-f46d-4032-a8b8-cdfbaf5317fc
abstract type AbstractPoint end

# ╔═╡ 253a9547-a2d4-4d17-b3b8-22194233bed3
md"""
To compute the norm of an `AbstractPoint`, we assume that all concrete subtypes must have fields called `x` and `y`. Of course, it is up to us to enforce this convention.
"""

# ╔═╡ 8178e06d-0632-4600-803a-09ed96816f61
mynorm(p::AbstractPoint) = sqrt(p.x^2 + p.y^2)

# ╔═╡ 3f9a432e-bab3-4357-b834-a2aaebe9fe31
md"""
**Version 1**

A struct written by a Julia beginner may look somewhat like this.
"""

# ╔═╡ 22b04135-f762-4331-8091-c8c3fa46655f
struct BeginnerPoint <: AbstractPoint
    x
    y
end

# ╔═╡ 3683d09a-7799-4bef-9d59-93f7fdb767a5
p_beginner = BeginnerPoint(3., 5.)

# ╔═╡ 9757e3ab-ecff-49e4-8fd9-44633e49b95c
@btime mynorm($p_beginner);

# ╔═╡ 19e2af3a-c409-4c9e-afa9-8874750ae909
md"""
Again, it is hard to judge performance in absolute terms, which is why we need profiling. Beware that, in order to profile such a short function, we need to run it many many times.
"""

# ╔═╡ 0ed838d3-32bc-4f40-82a7-066d50746f51
@profview for i = 1:1000000; mynorm(p_beginner); end

# ╔═╡ 76842b03-b2f1-482f-9982-e8903e35cb25
md"""
This time, the problem is not garbage collection, but runtime dispatch. Let's see what a type inference diagnosis can tell us.
"""

# ╔═╡ 9063e65e-15ef-420a-94a4-28a0b1f5335b
with_terminal() do
	@code_warntype mynorm(p_beginner)
end

# ╔═╡ d35a4f16-b5d4-4827-9b45-dbe28c9c4ff0
md"""
There are many red annotations in the `@code_warntype` output, meaning that the types of several intermediate variables cannot be inferred.
For instance, the line
```
%1  = Base.getproperty(p, :x)::Any
```
means that the field `p.x` has a type that cannot be determined during compilation alone.
More problematically, the line
```
Body::Any
```
states that the return type of the method is itself unknown. This means that type uncertainties may propagate if `mysqnorm(p_beginner)` is part of a larger code.
"""

# ╔═╡ 848fefa1-824b-4076-8149-b3a8869c172a
 @report_opt mynorm(p_beginner)

# ╔═╡ 23c83abe-0904-4faf-b5c7-e6f04b30da71
md"""
These problems are confirmed by the report of JET.jl, which detects several occurrences of runtime dispatch.
"""

# ╔═╡ 0a1dd5c2-d164-4b88-aa5d-a73ede91c56c
md"""
**Version 2**

The natural way to fix our first version would be to add type annotations to both fields. Indeed, without annotations, each field is considered of type `Any`, which is really bad.
"""

# ╔═╡ 40d777cc-7cf0-44f7-b179-fe3abbf4e030
struct CleverPoint <: AbstractPoint
    x::Float64
    y::Float64
end

# ╔═╡ 8b1b31e8-1f7c-427c-b69b-9fa5d4f654cc
p_clever = CleverPoint(3., 5.)

# ╔═╡ a13f3093-a2a7-441f-acaf-c4b9b099024c
@btime mynorm($p_clever);

# ╔═╡ 9f14261e-6bb9-4426-ae99-26fa35e531c1
md"""
As we can see, performance has greatly improved, because the method can now be fully inferred.
"""

# ╔═╡ 5970836a-5f14-446d-b05f-5beec9b05f8a
with_terminal() do
	@code_warntype mynorm(p_clever)
end

# ╔═╡ 7d487376-9651-45c6-bc8a-21117af8e745
 @report_opt mynorm(p_clever)

# ╔═╡ 5b5b3949-2a82-415a-8e2d-6b497c257a3f
md"""
But what if our points have other coordinate types? Maybe we don't want to convert them to `Float64` by default?
"""

# ╔═╡ bd06e581-1757-43f2-bdef-0fe4c8f9d238
md"""
**Version 3**

The most generic way to encode a point with real coordinates is to use a [parametric type](https://docs.julialang.org/en/v1/manual/types/#Parametric-Types).
"""

# ╔═╡ bb734c3b-d981-4473-aa04-9262206ee746
struct GeniusPoint{R <: Real} <: AbstractPoint
    x::R
    y::R
end

# ╔═╡ 1756c3bc-8662-4f76-bc6a-1b7448b36913
p_genius = GeniusPoint(3., 5.)

# ╔═╡ 008dcb2f-d32b-425d-bc7c-55b512d53b8a
@btime mynorm($p_genius)

# ╔═╡ f08f0153-e8c0-4ece-8cb6-5083539fb36c
md"""
`GeniusPoint` is just as fast as `CleverPoint` when applied to coordinates of type `Float64`. However, the former is generic: it can adapt to any coordinate type, which is often useful. A well-known example is [forward-mode automatic differentiation](https://en.wikipedia.org/wiki/Automatic_differentiation#Forward_accumulation), which uses dual numbers instead of standard floats.
"""

# ╔═╡ f34276f4-b267-4369-8563-1e1abe363a5f
ForwardDiff.gradient(a -> mynorm(GeniusPoint(a[1], a[2])), [3., 5.])

# ╔═╡ efd5cf6a-68e3-44b3-9b6f-eae396901e4e
md"""
We can compare this behavior with that of `CleverPoint`, for which the `Float64` conversion throws an error.
"""

# ╔═╡ 5c3eb0ba-dfef-4faa-87c5-009317b6faaa
try
	ForwardDiff.gradient(a -> mynorm(CleverPoint(a[1], a[2])), [3., 5.])
catch e
	@error e
end

# ╔═╡ fe04e854-1393-42fc-b6d7-6a4b3848e0ef
md"""
# Going further
"""

# ╔═╡ 41268586-52ac-47f3-8f53-b52072a9ae46
md"""
## Other performance tutorials

Here are some of my inspirations for this post:
- [Performance tips](https://www.juliafordatascience.com/performance-tips/) (Josh Day)
- [Optimizing Julia code](https://huijzer.xyz/posts/inference/) (Rik Huijzer)
- [7 Julia gotchas and how to handle them](https://www.stochasticlifestyle.com/7-julia-gotchas-handle/) (Chris Rackauckas)
"""

# ╔═╡ 107a2cec-d4bf-494b-b0c7-0d7e038369f6
md"""
## Good practices

Scientific programming needs fast code, but that's not all it needs. Adequate algorithms and sound coding habits are also very important. I recommend Antoine Levitt's notes on [Practical scientific computing](http://antoine.levitt.fr/calsci.pdf) as a starting point. 
"""

# ╔═╡ fdf97758-26c1-4157-a5d1-af89578f6277
md"""
## Generic programming

Multiple dispatch allows the right method to be chosen based on the type of every argument (not just the first one).
This is what makes it possible for multiple packages to work together without knowing about each other... as long as they remain generic.
In particular, it is not a good idea to overspecify input types: it usually doesn't improve performance, and can prevent unexpected uses of your code.

This is explained in great detail in the blog post [Type-Dispatch Design: Post Object-Oriented Programming for Julia](https://www.stochasticlifestyle.com/type-dispatch-design-post-object-oriented-programming-julia/) (Chris Rackauckas).
"""

# ╔═╡ 6437292a-2922-4219-a5e9-b7c8e2db20c7
md"""
## Hardware considerations

In order to optimize Julia code to the limit, it quickly becomes useful to know how a modern computer works.
The following blog post is an absolute masterpiece on this topic: [What scientists must know about hardware to write fast code](https://viralinstruction.com/posts/hardware/) (Jakob Nybo Nissen).
"""

# ╔═╡ ada6d5f4-f5fc-4c5f-9724-d29f4bb2a06a
md"""
## Package latency

A major source of frustration for Julia beginners is the time that elapses from the creation of the REPL until the first useful output, also called "time to first plot".
As a package developer, there are many resources available to help you address this problem:
- Tim Holy's great talk at JuliaCon 2021: [Package latency and what developers can do to reduce it](https://youtu.be/rVBgrWYKLHY)
- Several blog posts:
  - [Analyzing sources of compiler latency in Julia: method invalidations ](https://julialang.org/blog/2020/08/invalidations/) (Tim Holy, Jeff Bezanson, and Jameson Nash)
  - [Tutorial on precompilation](https://julialang.org/blog/2021/01/precompile_tutorial/) (Tim Holy)
  - [Profiling type inference](https://julialang.org/blog/2021/01/snoopi_deep/) (Tim Holy and Nathan Daly)
"""

# ╔═╡ 154d79cd-0925-4742-8a59-199e12014f0a
Multiple dispatch allows the right method to be chosen based on the type of every argument (not just the first one). This is what makes it possible for multiple packages to work together without knowing about each other... as long as they remain generic. In particular, it is not a good idea to overspecify input types: it usually doesn't improve performance, and can prevent unexpected uses of your code.

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
BenchmarkTools = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
JET = "c3a54625-cd67-489e-a8e7-0a5a0ff4e31b"
LoopVectorization = "bdcacae8-1622-11e9-2a5c-532679323890"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
ProfileCanvas = "efd6af41-a80b-495e-886c-e51b0c7d77a3"
ProgressLogging = "33c8b6b6-d38a-422a-b730-caa89a2f386c"

[compat]
BenchmarkTools = "~1.3.1"
ForwardDiff = "~0.10.30"
JET = "~0.5.16"
LoopVectorization = "~0.12.117"
PlutoUI = "~0.7.39"
ProfileCanvas = "~0.1.1"
ProgressLogging = "~0.1.4"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.3"
manifest_format = "2.0"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "8eaf9f1b4921132a4cff3f36a1d9ba923b14a481"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.4"

[[deps.Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "af92965fb30777147966f58acb05da51c5616b5f"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.3.3"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[deps.ArrayInterface]]
deps = ["ArrayInterfaceCore", "Compat", "IfElse", "LinearAlgebra", "Static"]
git-tree-sha1 = "8169fed5a73eff6978749e8be2b8185bb7a0fde3"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "6.0.11"

[[deps.ArrayInterfaceCore]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "d0f59ebfe8d3ea2799fb3fb88742d69978e5843e"
uuid = "30b0a656-2188-435a-8636-2ec0e6a096e2"
version = "0.1.10"

[[deps.ArrayInterfaceOffsetArrays]]
deps = ["ArrayInterface", "OffsetArrays", "Static"]
git-tree-sha1 = "3cbe45d8cc9cff51f302df1f87df64095423fd96"
uuid = "015c0d05-e682-4f19-8f0a-679ce4c54826"
version = "0.1.2"

[[deps.ArrayInterfaceStaticArrays]]
deps = ["Adapt", "ArrayInterface", "LinearAlgebra", "Static", "StaticArrays"]
git-tree-sha1 = "d7dc30474e73173a990eca86af76cae8790fa9f2"
uuid = "b0d46f97-bff5-4637-a19a-dd75974142cd"
version = "0.1.2"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.BenchmarkTools]]
deps = ["JSON", "Logging", "Printf", "Profile", "Statistics", "UUIDs"]
git-tree-sha1 = "4c10eee4af024676200bc7752e536f858c6b8f93"
uuid = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
version = "1.3.1"

[[deps.BitTwiddlingConvenienceFunctions]]
deps = ["Static"]
git-tree-sha1 = "28bbdbf0354959db89358d1d79d421ff31ef0b5e"
uuid = "62783981-4cbd-42fc-bca8-16325de8dc4b"
version = "0.1.3"

[[deps.CPUSummary]]
deps = ["CpuId", "IfElse", "Static"]
git-tree-sha1 = "0eaf4aedad5ccc3e39481db55d72973f856dc564"
uuid = "2a0fbf3d-bb9c-48f3-b0a9-814d99fd7ab9"
version = "0.1.22"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "9489214b993cd42d17f44c36e359bf6a7c919abf"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.15.0"

[[deps.ChangesOfVariables]]
deps = ["ChainRulesCore", "LinearAlgebra", "Test"]
git-tree-sha1 = "1e315e3f4b0b7ce40feded39c73049692126cf53"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.3"

[[deps.CloseOpenIntervals]]
deps = ["ArrayInterface", "Static"]
git-tree-sha1 = "eb61d6b97041496058245821e3bb7eba2b2cf4db"
uuid = "fb6a15b2-703c-40df-9091-08a04967cfa9"
version = "0.1.8"

[[deps.CodeTracking]]
deps = ["InteractiveUtils", "UUIDs"]
git-tree-sha1 = "6d4fa04343a7fc9f9cb9cff9558929f3d2752717"
uuid = "da1fd8a2-8d9e-5ec2-8556-3022fb5608a2"
version = "1.0.9"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "0f4e115f6f34bbe43c19751c90a38b2f380637b9"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.3"

[[deps.CommonSubexpressions]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "7b8a93dba8af7e3b42fecabf646260105ac373f7"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.0"

[[deps.Compat]]
deps = ["Dates", "LinearAlgebra", "UUIDs"]
git-tree-sha1 = "924cdca592bc16f14d2f7006754a621735280b74"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.1.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[deps.CpuId]]
deps = ["Markdown"]
git-tree-sha1 = "fcbb72b032692610bfbdb15018ac16a36cf2e406"
uuid = "adafc99b-e345-5852-983c-f28acb93d879"
version = "0.3.1"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DiffResults]]
deps = ["StaticArrays"]
git-tree-sha1 = "c18e98cba888c6c25d1c3b048e4b3380ca956805"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.0.3"

[[deps.DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "28d605d9a0ac17118fe2c5e9ce0fbb76c3ceb120"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.11.0"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "b19534d1895d702889b219c382a6e18010797f0b"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.8.6"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions", "StaticArrays"]
git-tree-sha1 = "2f18915445b248731ec5db4e4a17e451020bf21e"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.30"

[[deps.HostCPUFeatures]]
deps = ["BitTwiddlingConvenienceFunctions", "IfElse", "Libdl", "Static"]
git-tree-sha1 = "18be5268cf415b5e27f34980ed25a7d34261aa83"
uuid = "3e5b6fbb-0976-4d2c-9146-d79de83f2fb0"
version = "0.1.7"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "c47c5fa4c5308f27ccaac35504858d8914e102f9"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.4"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "f7be53659ab06ddc986428d3a9dcc95f6fa6705a"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.2"

[[deps.IfElse]]
git-tree-sha1 = "debdd00ffef04665ccbb3e150747a77560e8fad1"
uuid = "615f187c-cbe4-4ef1-ba3b-2fcf58d6d173"
version = "0.1.1"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "c6cf981474e7094ce044168d329274d797843467"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.6"

[[deps.IrrationalConstants]]
git-tree-sha1 = "7fd44fd4ff43fc60815f8e764c0f352b83c49151"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.1"

[[deps.JET]]
deps = ["InteractiveUtils", "JuliaInterpreter", "LoweredCodeUtils", "MacroTools", "Pkg", "Revise", "Test"]
git-tree-sha1 = "8e78b0c297cfa6cefd579f87232c89bd6ed7a081"
uuid = "c3a54625-cd67-489e-a8e7-0a5a0ff4e31b"
version = "0.5.16"

[[deps.JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "abc9885a7ca2052a736a600f7fa66209f96506e1"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "3c837543ddb02250ef42f4738347454f95079d4e"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.3"

[[deps.JuliaInterpreter]]
deps = ["CodeTracking", "InteractiveUtils", "Random", "UUIDs"]
git-tree-sha1 = "52617c41d2761cc05ed81fe779804d3b7f14fff7"
uuid = "aa1ae85d-cabe-5617-a682-6adf51b2e16a"
version = "0.9.13"

[[deps.LayoutPointers]]
deps = ["ArrayInterface", "ArrayInterfaceOffsetArrays", "ArrayInterfaceStaticArrays", "LinearAlgebra", "ManualMemory", "SIMDTypes", "Static"]
git-tree-sha1 = "a575de5a424a395217930fea6d0934ea853d0158"
uuid = "10f19ff3-798f-405d-979b-55457f8fc047"
version = "0.1.9"

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

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "09e4b894ce6a976c354a69041a04748180d43637"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.15"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.LoopVectorization]]
deps = ["ArrayInterface", "ArrayInterfaceCore", "ArrayInterfaceOffsetArrays", "ArrayInterfaceStaticArrays", "CPUSummary", "ChainRulesCore", "CloseOpenIntervals", "DocStringExtensions", "ForwardDiff", "HostCPUFeatures", "IfElse", "LayoutPointers", "LinearAlgebra", "OffsetArrays", "PolyesterWeave", "SIMDDualNumbers", "SIMDTypes", "SLEEFPirates", "SpecialFunctions", "Static", "ThreadingUtilities", "UnPack", "VectorizationBase"]
git-tree-sha1 = "ebfce259dfe319b404b9cddda0648e042f3192ac"
uuid = "bdcacae8-1622-11e9-2a5c-532679323890"
version = "0.12.117"

[[deps.LoweredCodeUtils]]
deps = ["JuliaInterpreter"]
git-tree-sha1 = "dedbebe234e06e1ddad435f5c6f4b85cd8ce55f7"
uuid = "6f1432cf-f94c-5a45-995e-cdbf5db27b0b"
version = "2.2.2"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "3d3e902b31198a27340d0bf00d6ac452866021cf"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.9"

[[deps.ManualMemory]]
git-tree-sha1 = "bcaef4fc7a0cfe2cba636d84cda54b5e4e4ca3cd"
uuid = "d125e4d3-2237-4719-b19c-fa641b8a4667"
version = "0.1.8"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[deps.NaNMath]]
git-tree-sha1 = "737a5957f387b17e74d4ad2f440eb330b39a62c5"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.0"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[deps.OffsetArrays]]
deps = ["Adapt"]
git-tree-sha1 = "b4975062de00106132d0b01b5962c09f7db7d880"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.12.5"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[deps.Parsers]]
deps = ["Dates"]
git-tree-sha1 = "1285416549ccfcdf0c50d4997a94331e88d68413"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.3.1"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "UUIDs"]
git-tree-sha1 = "8d1f54886b9037091edf146b517989fc4a09efec"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.39"

[[deps.PolyesterWeave]]
deps = ["BitTwiddlingConvenienceFunctions", "CPUSummary", "IfElse", "Static", "ThreadingUtilities"]
git-tree-sha1 = "7e597df97e46ffb1c8adbaddfa56908a7a20194b"
uuid = "1d0040c9-8b98-4ee7-8388-3f51789ca0ad"
version = "0.1.5"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "47e5f437cc0e7ef2ce8406ce1e7e24d44915f88d"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.3.0"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.Profile]]
deps = ["Printf"]
uuid = "9abbd945-dff8-562f-b5e8-e1ebf5ef1b79"

[[deps.ProfileCanvas]]
deps = ["JSON", "Pkg", "Profile", "REPL"]
git-tree-sha1 = "eb5b798f906389b255247a22c5eff8b554fbedca"
uuid = "efd6af41-a80b-495e-886c-e51b0c7d77a3"
version = "0.1.1"

[[deps.ProgressLogging]]
deps = ["Logging", "SHA", "UUIDs"]
git-tree-sha1 = "80d919dee55b9c50e8d9e2da5eeafff3fe58b539"
uuid = "33c8b6b6-d38a-422a-b730-caa89a2f386c"
version = "0.1.4"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.Revise]]
deps = ["CodeTracking", "Distributed", "FileWatching", "JuliaInterpreter", "LibGit2", "LoweredCodeUtils", "OrderedCollections", "Pkg", "REPL", "Requires", "UUIDs", "Unicode"]
git-tree-sha1 = "4d4239e93531ac3e7ca7e339f15978d0b5149d03"
uuid = "295af30f-e4ad-537b-8983-00126c2a3abe"
version = "3.3.3"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[deps.SIMDDualNumbers]]
deps = ["ForwardDiff", "IfElse", "SLEEFPirates", "VectorizationBase"]
git-tree-sha1 = "dd4195d308df24f33fb10dde7c22103ba88887fa"
uuid = "3cdde19b-5bb0-4aaf-8931-af3e248e098b"
version = "0.1.1"

[[deps.SIMDTypes]]
git-tree-sha1 = "330289636fb8107c5f32088d2741e9fd7a061a5c"
uuid = "94e857df-77ce-4151-89e5-788b33177be4"
version = "0.1.0"

[[deps.SLEEFPirates]]
deps = ["IfElse", "Static", "VectorizationBase"]
git-tree-sha1 = "ac399b5b163b9140f9c310dfe9e9aaa225617ff6"
uuid = "476501e8-09a2-5ece-8869-fb82de89a1fa"
version = "0.6.32"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "a9e798cae4867e3a41cae2dd9eb60c047f1212db"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.1.6"

[[deps.Static]]
deps = ["IfElse"]
git-tree-sha1 = "5d2c08cef80c7a3a8ba9ca023031a85c263012c5"
uuid = "aedffcd0-7271-4cad-89d0-dc628f76c6d3"
version = "0.6.6"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "383a578bdf6e6721f480e749d503ebc8405a0b22"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.4.6"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.ThreadingUtilities]]
deps = ["ManualMemory"]
git-tree-sha1 = "f8629df51cab659d70d2e5618a430b4d3f37f2c3"
uuid = "8290d209-cae3-49c0-8002-c8c24d57dab5"
version = "0.5.0"

[[deps.Tricks]]
git-tree-sha1 = "6bac775f2d42a611cdfcd1fb217ee719630c4175"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.6"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.VectorizationBase]]
deps = ["ArrayInterface", "CPUSummary", "HostCPUFeatures", "IfElse", "LayoutPointers", "Libdl", "LinearAlgebra", "SIMDTypes", "Static"]
git-tree-sha1 = "7d3de169cd221392082a5abc7f363726e1a30628"
uuid = "3d5dd08c-fd9d-11e8-17fa-ed2836048c2f"
version = "0.21.36"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
"""

# ╔═╡ Cell order:
# ╟─155a8c77-2837-44ee-9b7a-68d5a89a2fdc
# ╟─c48f5b63-48b0-4bb6-997d-7d40745a46fa
# ╠═86e4051d-ede6-4e1d-9676-98d92140cc78
# ╠═27793eba-f876-4e92-849d-945989e39216
# ╠═ac232fbf-827e-4363-9eac-04885334368a
# ╠═144e4f59-c0d0-4f75-aa49-d05cf781d83e
# ╠═6b3592a6-43d0-4ebd-a06d-d325dd66a829
# ╠═06e6e4b8-31b2-4f9e-9da3-093a030021af
# ╠═22d3e9db-f4ec-4810-b2d5-94bdbf8d52ac
# ╠═5de2a556-f3af-4a64-a5c6-32d30f758be3
# ╠═1ac5ba38-0eef-41bb-8f9c-3bbf057cae21
# ╟─e1852c8d-4028-409e-8e1a-8253bbd6e6a5
# ╟─9331fad2-f29e-11eb-0349-477bd2e7e412
# ╟─3d98e7db-c643-4500-987d-4a225e55b2a5
# ╠═b4f2a99e-de45-49d2-be86-9f2d03357462
# ╟─068b3e45-5105-48aa-a547-536470f6abda
# ╠═7a75c990-46b4-484e-acc4-65e34f41a9f2
# ╟─f7b1b44f-2aa6-4c5c-97a2-ac7037fb48ce
# ╠═1fb43343-083b-4b1a-b622-d88c9aa0808c
# ╠═a28f7911-3dbb-45fb-a82d-2834d3c8502c
# ╠═4da8a7ca-3cea-4629-a66d-44f3b907af09
# ╟─c0a7c1fe-457f-4e52-b0ea-2821e40817ea
# ╠═4e73e893-31c5-4c34-82da-ddabaa7316a4
# ╟─ddbe8924-756a-46a8-9af1-6df0a21e9057
# ╟─94c78148-c651-4a59-9e62-5c7e9576d1e8
# ╠═c44f3dc9-ff19-4ba4-9388-73cfaf23f8e8
# ╟─a7de0ec9-6b01-4b42-8cce-bb2295da779f
# ╟─9ca598c1-dae8-40b9-a18d-c74f30524b35
# ╟─0fb6ed33-601c-4392-b7d9-32230c979d39
# ╟─a6e9da76-1ff0-4b54-9b55-4856ca32b251
# ╟─d3c1a86c-8c8f-4ad6-ac3c-2ba0f838d139
# ╟─fa483fea-bf9f-4764-8d4f-c6d33e3336fb
# ╟─00e96665-ccbb-42c1-b5de-f626d57bb0ff
# ╟─43bad028-9d16-426f-9cdb-a37b1ee1a623
# ╠═857509a7-f07a-4bf0-9383-207984b95faa
# ╟─7386749b-b2ab-48a7-a1d2-46e7f31e72e3
# ╠═6cd86e7a-9f82-4da1-a8f0-4ed2c1068ab9
# ╠═5d10a00b-bfa9-49c7-9f4b-503351fa2842
# ╟─38cc6383-c7d8-46b4-8531-251bd196d960
# ╠═62163e17-4265-4c97-95bb-29d608e80b07
# ╟─9d8b7e25-c9c6-4aba-a33c-66fd18d804c0
# ╟─fe45168c-8cf1-435e-86fc-16cfffef3ec1
# ╠═0400175c-5a3c-44a7-9a8a-c30a4756b88c
# ╠═cd0cc22f-2d4d-4848-8f15-8f0127a4245b
# ╟─638de554-1bec-453d-9e30-796247aaa4cc
# ╠═fd4401cf-69e8-4444-92c3-478035301006
# ╟─23053665-e058-43de-95d9-c688e3a80b0c
# ╟─9a181530-02e7-47b0-9a86-c191baefac54
# ╠═c171555a-0166-476e-8ec6-1860745d84f2
# ╠═610f6d6f-9d37-4f3d-be78-ab9847162f4d
# ╟─c7b551a0-8c2e-4785-b575-8d58e37c14ec
# ╠═0d8af577-9275-490c-a689-65e7177c4d65
# ╟─69e8bf4e-d98d-4804-b6bf-f299c3452565
# ╠═dacdb662-f46d-4032-a8b8-cdfbaf5317fc
# ╟─253a9547-a2d4-4d17-b3b8-22194233bed3
# ╠═8178e06d-0632-4600-803a-09ed96816f61
# ╟─3f9a432e-bab3-4357-b834-a2aaebe9fe31
# ╠═22b04135-f762-4331-8091-c8c3fa46655f
# ╠═3683d09a-7799-4bef-9d59-93f7fdb767a5
# ╠═9757e3ab-ecff-49e4-8fd9-44633e49b95c
# ╟─19e2af3a-c409-4c9e-afa9-8874750ae909
# ╠═0ed838d3-32bc-4f40-82a7-066d50746f51
# ╟─76842b03-b2f1-482f-9982-e8903e35cb25
# ╠═9063e65e-15ef-420a-94a4-28a0b1f5335b
# ╟─d35a4f16-b5d4-4827-9b45-dbe28c9c4ff0
# ╠═848fefa1-824b-4076-8149-b3a8869c172a
# ╟─23c83abe-0904-4faf-b5c7-e6f04b30da71
# ╟─0a1dd5c2-d164-4b88-aa5d-a73ede91c56c
# ╠═40d777cc-7cf0-44f7-b179-fe3abbf4e030
# ╠═8b1b31e8-1f7c-427c-b69b-9fa5d4f654cc
# ╠═a13f3093-a2a7-441f-acaf-c4b9b099024c
# ╟─9f14261e-6bb9-4426-ae99-26fa35e531c1
# ╠═5970836a-5f14-446d-b05f-5beec9b05f8a
# ╠═7d487376-9651-45c6-bc8a-21117af8e745
# ╟─5b5b3949-2a82-415a-8e2d-6b497c257a3f
# ╟─bd06e581-1757-43f2-bdef-0fe4c8f9d238
# ╠═bb734c3b-d981-4473-aa04-9262206ee746
# ╠═1756c3bc-8662-4f76-bc6a-1b7448b36913
# ╠═008dcb2f-d32b-425d-bc7c-55b512d53b8a
# ╟─f08f0153-e8c0-4ece-8cb6-5083539fb36c
# ╠═f34276f4-b267-4369-8563-1e1abe363a5f
# ╟─efd5cf6a-68e3-44b3-9b6f-eae396901e4e
# ╠═5c3eb0ba-dfef-4faa-87c5-009317b6faaa
# ╟─fe04e854-1393-42fc-b6d7-6a4b3848e0ef
# ╟─41268586-52ac-47f3-8f53-b52072a9ae46
# ╟─107a2cec-d4bf-494b-b0c7-0d7e038369f6
# ╟─fdf97758-26c1-4157-a5d1-af89578f6277
# ╟─6437292a-2922-4219-a5e9-b7c8e2db20c7
# ╟─ada6d5f4-f5fc-4c5f-9724-d29f4bb2a06a
# ╠═154d79cd-0925-4742-8a59-199e12014f0a
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
